import sys
import os
from pymol import cmd
from pymol.Qt import QtWidgets, utils, QtCore
import mdtraj as md
import json
import copy
import webbrowser
import re

from fretraj import cloud
from fretraj import isosurf

# class App(QtWidgets.QWidget):


class App(QtWidgets.QMainWindow):

    def __init__(self, _pymol_running=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.uiIconPath = '../docs/source/_static/fretraj_logo'
        utils.loadUi('fretraj.ui', self)
        self.setWindowIcon(utils.QtGui.QIcon(self.uiIconPath))
        self._pymol_running = _pymol_running
        self.statusBar().showMessage("Ready", 2000)
        self.docsPath = '../docs/index.html'

        # initialize labels dictionary with defaults from GUI
        self.struct = None
        self.av = {}
        self.traj = {}
        self.labels = {'Position': {}, 'Distance': {}}
        self.labelName = self.comboBox_labelName.currentText()
        self.labelName_default = self.comboBox_labelName.currentText()
        self.distanceName = self.comboBox_distanceName.currentText()
        self.distanceName_default = self.comboBox_distanceName.currentText()
        self.update_labelDict()
        self.labels_default = copy.copy(self.labels)

        # set root path
        self.settings = {'root_path': None, 'browser': None}
        if os.path.isfile('fretraj_settings.conf'):
            with open('fretraj_settings.conf', 'r') as f:
                self.settings = json.load(f)
                self.lineEdit_rootDirectory.setText(self.settings['root_path'])
        else:
            self.setRootDirectory()
            if not self.settings['root_path']:
                sys.exit()

        # activate / deactivate GUI elements
        self.dyeRadius_spinBox_OnOff()
        self.push_computeACV.setEnabled(False)
        self.doubleSpinBox_isoValue.setEnabled(False)
        self.push_calculateFRET.setEnabled(False)
        self.spinBox_statePDB.setEnabled(False)

        # signals
        self.push_computeACV.clicked.connect(self.computeACV)
        self.push_loadPDB.clicked.connect(self.loadPDB)
        self.comboBox_simType.activated.connect(self.dyeRadius_spinBox_OnOff)
        self.push_calculateFRET.clicked.connect(self.calculateFRET)
        self.push_deleteLabel.clicked.connect(self.deleteLabel)
        self.push_deleteFRETparam.clicked.connect(self.deleteFRETparam)
        self.push_loadParameterFile.clicked.connect(self.loadParameterFile)
        self.doubleSpinBox_isoValue.valueChanged.connect(self.pymol_update_isosurface)
        self.lineEdit_labelName.returnPressed.connect(self.makeLabel)
        self.lineEdit_distanceName.returnPressed.connect(self.makeFRETparam)
        self.comboBox_labelName.currentIndexChanged.connect(self.update_GUIfields)
        self.actionAbout_FRETraj.triggered.connect(self.openAbout)
        self.comboBox_donorName.currentIndexChanged.connect(self.define_DA)
        self.comboBox_acceptorName.currentIndexChanged.connect(self.define_DA)
        self.comboBox_distanceName.currentIndexChanged.connect(self.update_GUIfields)
        self.spinBox_statePDB.valueChanged.connect(self.update_PDBstate)
        self.push_setRootDirectory.clicked.connect(self.setRootDirectory)
        self.actionDocumentation.triggered.connect(self.openDocumentation)

    def update_labelDict(self):
        """
        Update the label dictionary with the values from the GUI fields
        """

        pos = self.labelName
        dis = self.distanceName
        self.labels['Position'][pos] = {}
        self.labels['Distance'][dis] = {}
        self.labels['Position'][pos]['attach_id'] = self.spinBox_atomID.value()
        self.labels['Position'][pos]['linker_length'] = self.doubleSpinBox_linkerLength.value()
        self.labels['Position'][pos]['linker_width'] = self.doubleSpinBox_linkerWidth.value()
        self.labels['Position'][pos]['cv_thickness'] = self.doubleSpinBox_CVthickness.value()
        self.labels['Position'][pos]['cv_fraction'] = self.doubleSpinBox_CVfraction.value()
        self.labels['Position'][pos]['simulation_type'] = self.comboBox_simType.currentText()
        self.labels['Position'][pos]['grid_spacing'] = self.doubleSpinBox_gridSpacing.value()
        self.labels['Position'][pos]['dye_radius1'] = self.doubleSpinBox_dyeRadius1.value()
        self.labels['Position'][pos]['dye_radius2'] = self.doubleSpinBox_dyeRadius2.value()
        self.labels['Position'][pos]['dye_radius3'] = self.doubleSpinBox_dyeRadius3.value()
        self.labels['Position'][pos]['mol_selection'] = self.lineEdit_molSelection.text()
        self.labels['Distance'][dis]['R0'] = self.doubleSpinBox_R0.value()
        self.labels['Distance'][dis]['n_dist'] = self.spinBox_nDist.value()
        self.labels['Distance'][dis]['use_LabelLib'] = self.checkBox_useLabelLib.isChecked()
        self.frame = self.spinBox_statePDB.value()
        self.donorName = self.comboBox_donorName.currentText()
        self.acceptorName = self.comboBox_acceptorName.currentText()
        self._frame_pos = '{:d}-{}'.format(self.frame, self.labelName)

    def update_GUIfields(self):
        """
        Update the fields of the GUI upon changing the label in the dropdown
        """
        self.update_labelDict()
        pos = self.comboBox_labelName.currentText()
        dis = self.comboBox_distanceName.currentText()
        self.labelName = pos
        self.distanceName = dis
        self.spinBox_atomID.setValue(self.labels['Position'][pos]['attach_id'])
        self.doubleSpinBox_linkerLength.setValue(self.labels['Position'][pos]['linker_length'])
        self.doubleSpinBox_linkerWidth.setValue(self.labels['Position'][pos]['linker_width'])
        self.doubleSpinBox_CVthickness.setValue(self.labels['Position'][pos]['cv_thickness'])
        self.doubleSpinBox_CVfraction.setValue(self.labels['Position'][pos]['cv_fraction'])
        self.comboBox_simType.setCurrentText(self.labels['Position'][pos]['simulation_type'])
        self.dyeRadius_spinBox_OnOff()
        self.doubleSpinBox_gridSpacing.setValue(self.labels['Position'][pos]['grid_spacing'])
        self.doubleSpinBox_dyeRadius1.setValue(self.labels['Position'][pos]['dye_radius1'])
        self.doubleSpinBox_dyeRadius2.setValue(self.labels['Position'][pos]['dye_radius2'])
        self.doubleSpinBox_dyeRadius3.setValue(self.labels['Position'][pos]['dye_radius3'])
        self.lineEdit_molSelection.setText(self.labels['Position'][pos]['mol_selection'])
        self.doubleSpinBox_R0.setValue(self.labels['Distance'][dis]['R0'])
        self.spinBox_nDist.setValue(self.labels['Distance'][dis]['n_dist'])
        self.checkBox_useLabelLib.setChecked(self.labels['Distance'][dis]['use_LabelLib'])
        self.spinBox_statePDB.setValue(self.frame)

    def loadPDB(self):
        """
        Load PDB of CIF file. CIF files will be converted to PDB internally
        """
        self.fileNamePath_pdb, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load PDB / CIF', '', "PDB / CIF file (*.pdb *cif);;All Files (*)")
        if self.fileNamePath_pdb:
            self.fileName_pdb = self.fileNamePath_pdb.split("/")[-1]
            if self.fileName_pdb[-3:] == 'cif':
                self.fileNamePath_pdb, _ = cif2pdb(self.fileNamePath_pdb)
            try:
                self.struct = md.load_pdb(self.fileNamePath_pdb)
                self.spinBox_statePDB.setMaximum(self.struct.n_frames - 1)
            except IndexError:
                self.openErrorWin("PDB File Error", "The specified file \"{}.pdb\" can not be loaded".format(name))
            else:
                self.lineEdit_pdbFile.setText(self.fileName_pdb)
                self.push_computeACV.setEnabled(True)
                self.spinBox_statePDB.setEnabled(True)
                if self._pymol_running:
                    cmd.load(self.fileNamePath_pdb)

    def cif2pdb(self, pathtoPDBfile):
        """
        Convert .cif to .pdb with PyMOL

        Returns
        -------

        """
        name = pathtoPDBfile.split("/")[-1][:-4]
        try:
            cmd.load(pathtoPDBfile)
        except:
            self.openErrorWin("PDB File Error", "The specified file \"{}.pdb\" can not be loaded".format(name))
        else:
            fileNamePath_pdb = '{}/{}.pdb'.format(self.settings['root_path'], name)
            fileName_pdb = '{}.pdb'.format(fileNamePath_pdb.split("/")[-1])
            cmd.save(fileNamePath_pdb, name)
        return fileNamePath_pdb, fileName_pdb

    def update_PDBstate(self):
        self.update_labelDict()
        if self._pymol_running:
            cmd.set('state', self.frame)

    def loadParameterFile(self):
        """
        Load dye, linker and simulation settings from parameter file
        """
        self.fileNamePath_param, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load label parameter file', '', "JSON Files (*.json);;All Files (*)")
        if self.fileNamePath_param:
            self.fileName_param = self.fileNamePath_param.split("/")[-1]
            self.lineEdit_paramFile.setText(self.fileName_param)

            with open(self.fileNamePath_param) as f:
                try:
                    labels_json = json.load(f)
                except:
                    self.openErrorWin('Parameter File Error', "The specified file \"{}\" is not of type JSON".format(self.fileName_param))
                else:
                    error_msg = "The specified file \"{}\" has a wrong format".format(self.fileName_param)

                    for field in ['Position', 'Distance']:
                        if field == 'Position':
                            name_default = self.labelName_default
                        else:
                            name_default = self.distanceName_default
                            if 'Distance' not in labels_json:
                                labels_json['Distance'] = {self.distanceName_default: {}}
                        if 'Position' in labels_json:
                            if isinstance(labels_json['Position'], dict):
                                for pos in labels_json[field].keys():
                                    self.labels[field][pos] = {}

                                    for key in self.labels_default[field][name_default].keys():
                                        if isinstance(labels_json[field][pos], dict):
                                            if key in labels_json[field][pos].keys():
                                                self.labels[field][pos][key] = labels_json[field][pos][key]
                                            else:
                                                self.labels[field][pos][key] = copy.copy(self.labels_default[field][name_default][key])
                                        else:
                                            self.openErrorWin('Parameter File Error', error_msg)
                                            return
                            else:
                                self.openErrorWin('Parameter File Error', error_msg)
                                return
                        else:
                            self.openErrorWin('Parameter File Error', error_msg)
                            return

                    for newlabel in self.labels['Position'].keys():
                        self.addLabel(newlabel)
                    for newdistance in self.labels['Distance'].keys():
                        self.addFRETparam(newdistance)

    def makeLabel(self):
        """
        Create new label from edit box string
        """
        newlabel = self.lineEdit_labelName.text()
        self.labels['Position'][newlabel] = copy.copy(self.labels_default['Position'][self.labelName_default])
        self.addLabel(newlabel)

    def addLabel(self, newlabel):
        """
        Add a new label to Combobox
        """
        if self.comboBox_labelName.findText(newlabel) == -1:
            self.comboBox_labelName.addItem(newlabel)
            self.comboBox_labelName.setCurrentText(newlabel)

        self.lineEdit_labelName.clear()
        if self.struct:
            self.push_computeACV.setEnabled(True)
        if self.tableWidget_MeanPos.rowCount() > 1 and self.struct:
            self.push_calculateFRET.setEnabled(True)
        self.update_labelDict()

    def deleteLabel(self):
        """
        Remove the selected label from Combobox
        """
        self.deleteLabelFromList()
        self.deleteDistanceFromList()
        if self.labelName != self.labelName_default:
            self.comboBox_labelName.removeItem(self.comboBox_labelName.currentIndex())
            self.lineEdit_labelName.clear()

    def addLabelToList(self, av):
        """
        Add a label and the computed mean dye position (MP) to the list
        """
        rowCount = self.tableWidget_MeanPos.rowCount()
        for r in range(rowCount):
            if self._frame_pos == self.tableWidget_MeanPos.item(r, 0).text():
                break
        else:
            r = 0
            self.tableWidget_MeanPos.insertRow(r)
        self.tableWidget_MeanPos.setItem(r, 0, QtWidgets.QTableWidgetItem(self._frame_pos))
        self.tableWidget_MeanPos.setItem(r, 1, QtWidgets.QTableWidgetItem(', '.join('{:.2f}'.format(p) for p in av.acv.mp)))
        if self.comboBox_donorName.findText(self._frame_pos) == -1:
            self.comboBox_donorName.addItem(self._frame_pos)
        if self.comboBox_acceptorName.findText(self._frame_pos) == -1:
            self.comboBox_acceptorName.addItem(self._frame_pos)

    def deleteLabelFromList(self):
        """
        Remove a label and the computed mean dye position (MP) from the list
        """
        r = 0
        rowCount = self.tableWidget_MeanPos.rowCount()
        while r < rowCount:
            if self.labelName in self.tableWidget_MeanPos.item(r, 0).text():
                self.tableWidget_MeanPos.removeRow(r)
                rowCount -= 1
            else:
                r += 1

        i = 0
        itemCount = self.comboBox_donorName.count()
        while i < itemCount:
            if self.labelName in self.comboBox_donorName.itemText(i):
                self.comboBox_donorName.removeItem(i)
                self.comboBox_acceptorName.removeItem(i)
                itemCount -= 1
            else:
                i += 1

    def makeFRETparam(self):
        """
        Create new FRET parameters from edit box string
        """
        newdistance = self.lineEdit_distanceName.text()
        self.labels['Distance'][newdistance] = copy.copy(self.labels_default['Distance'][self.distanceName_default])
        self.addFRETparam(newdistance)

    def addFRETparam(self, newdistance):
        """
        Add FRET parameters to the Combobox
        """
        if self.comboBox_distanceName.findText(newdistance) == -1:
            self.comboBox_distanceName.addItem(newdistance)
            self.comboBox_distanceName.setCurrentText(newdistance)
        self.lineEdit_distanceName.clear()

    def deleteFRETparam(self):
        """
        Remove FRET parameters from the Combobox
        """
        dis = self.comboBox_distanceName.currentText()
        if dis != self.distanceName_default:
            self.comboBox_distanceName.removeItem(self.comboBox_distanceName.currentIndex())

    def addDistanceToList(self, traj):
        """
        Add a FRET set to the list
        """
        DA = '{} -> {}'.format(self.donorName, self.acceptorName)
        rowCount = self.tableWidget_FRET.rowCount()
        for r in range(rowCount):
            if DA == self.tableWidget_FRET.item(r, 0).text():
                break
        else:
            r = 0
            self.tableWidget_FRET.insertRow(r)
        self.tableWidget_FRET.setItem(r, 0, QtWidgets.QTableWidgetItem('{} -> {}'.format(self.donorName, self.acceptorName)))
        self.tableWidget_FRET.setItem(r, 1, QtWidgets.QTableWidgetItem('{:.1f}'.format(traj.mean_R_DA)))
        self.tableWidget_FRET.setItem(r, 2, QtWidgets.QTableWidgetItem('{:.1f}'.format(traj.sigma_R_DA)))
        self.tableWidget_FRET.setItem(r, 3, QtWidgets.QTableWidgetItem('{:.2f}'.format(traj.mean_E_DA)))
        self.tableWidget_FRET.setItem(r, 4, QtWidgets.QTableWidgetItem('{:.1f}'.format(traj.mean_R_DA_E)))
        self.tableWidget_FRET.setItem(r, 5, QtWidgets.QTableWidgetItem('{:.1f}'.format(traj.R_mp)))
        self.tableWidget_FRET.setItem(r, 6, QtWidgets.QTableWidgetItem('{:.1f}'.format(traj.R_attach)))

    def deleteDistanceFromList(self):
        """
        Remove a FRET set from the list
        """
        r = 0
        rowCount = self.tableWidget_FRET.rowCount()
        while r < rowCount:
            if self.labelName in self.tableWidget_FRET.item(r, 0).text():
                self.tableWidget_FRET.removeRow(r)
                rowCount -= 1
            else:
                r += 1

    def define_DA(self):
        """
        Define the donor and acceptor labels and color them in the mean position table
        """
        self.update_labelDict()
        for i in range(self.tableWidget_MeanPos.rowCount()):
            for j in range(self.tableWidget_MeanPos.columnCount()):
                self.tableWidget_MeanPos.item(i, j).setBackground(utils.QtGui.QColor(255, 255, 255))
        if self.tableWidget_MeanPos.rowCount() > 1:
            if self.donorName == self.acceptorName:
                self.push_calculateFRET.setEnabled(False)
            else:
                self.push_calculateFRET.setEnabled(True)
                row_don = self.tableWidget_MeanPos.findItems(self.donorName, QtCore.Qt.MatchExactly)[0].row()
                row_acc = self.tableWidget_MeanPos.findItems(self.acceptorName, QtCore.Qt.MatchExactly)[0].row()
                for j in range(self.tableWidget_MeanPos.columnCount()):
                    self.tableWidget_MeanPos.item(row_don, j).setBackground(utils.QtGui.QColor(102, 184, 99))
                    self.tableWidget_MeanPos.item(row_acc, j).setBackground(utils.QtGui.QColor(212, 76, 81))
        else:
            self.push_calculateFRET.setEnabled(False)

    def dyeRadius_spinBox_OnOff(self):
        """
        Activate or deactivate dyeRadius spinboxes depending on the type of simulation
        """
        if self.comboBox_simType.currentText() == 'AV1':
            self.doubleSpinBox_dyeRadius2.setEnabled(False)
            self.doubleSpinBox_dyeRadius3.setEnabled(False)
        else:
            self.doubleSpinBox_dyeRadius2.setEnabled(True)
            self.doubleSpinBox_dyeRadius3.setEnabled(True)

    def computeACV(self):
        """
        Create a cloud.Volume object at the specified label position.
        If the GUI is run as a PyMOL plugin, the ACV will be rendered in the current window
        """
        self.update_labelDict()
        self.av[self._frame_pos] = cloud.Volume(self.struct, self.frame, self.labelName, self.labels)
        if self.av[self._frame_pos].acv is None:
            msg = 'ACV could not be calculated!'
            self.statusBar().showMessage(msg, 3000)
            print(msg)
        else:
            av_filename = '{}/{}.xyz'.format(self.settings['root_path'], self._frame_pos)
            self.av[self._frame_pos].save_acv(av_filename, format='xyz')
            self.addLabelToList(self.av[self._frame_pos])
            self.define_DA()
            msg = 'ACV successfully calculated!'
            self.statusBar().showMessage(msg, 3000)
            if self._pymol_running:
                cmd.load(av_filename)
                self.doubleSpinBox_isoValue.setEnabled(False)
                isoval = self.doubleSpinBox_isoValue.value()
                isosurf.smooth_map_from_xyz(self.labelName, self.labelName, isoval)

    def calculateFRET(self):
        self.update_labelDict()
        self.traj[(self.donorName, self.acceptorName)] = cloud.FRET_Trajectory(self.av[self.donorName], self.av[self.acceptorName], R0=self.labels['Distance'][self.distanceName]['R0'], n_dist=self.labels['Distance'][self.distanceName]['n_dist'], use_LabelLib=self.labels['Distance'][self.distanceName]['use_LabelLib'])
        self.addDistanceToList(self.traj[(self.donorName, self.acceptorName)])

    def pymol_update_isosurface(self):
        """
        Update the isosurface when spin box value is changed.
        The spin box is only active when PYMOL is running.
        """
        isoval = self.doubleSpinBox_isoValue.value()
        isosurf.smooth_map_from_xyz(self.labelName, self.labelName, isoval)

    def setRootDirectory(self):
        rootDir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Set root directory')
        if rootDir:
            self.lineEdit_rootDirectory.setText(rootDir)
            self.settings['root_path'] = rootDir
            with open('fretraj_settings.conf', 'w') as f:
                json.dump(self.settings, f)

    def openDocumentation(self):
        if not self.settings['browser']:
            browser_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Search for default browser')
            self.settings['browser'] = re.escape(browser_path)
            with open('fretraj_settings.conf', 'w') as f:
                json.dump(self.settings, f)
        if self.settings['browser']:
            try:
                browser = webbrowser.get('{} %s'.format(self.settings['browser']))
                browser.open(self.docsPath)
            except:
                self.settings['browser'] = None
                with open('fretraj_settings.conf', 'w') as f:
                    json.dump(self.settings, f)

    def openAbout(self):
        """
        Open About Window
        """
        msg = QtWidgets.QMessageBox()
        pixmap = utils.QtGui.QPixmap(self.uiIconPath)
        msg.setIconPixmap(pixmap.scaledToWidth(64))
        msg.setWindowTitle("About FRETraj")
        msg.setText('FRETraj v.0.1\n\nFRETraj is program to calculate accessible contact volumes\n\n(C) Fabio Steffen\nUiversity of Zurich\n2018-2019')
        msg.exec_()

    def openErrorWin(self, title, message):
        """
        Open Error Window
        """
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Critical)
        msg.setWindowTitle(title)
        msg.setText(message)
        msg.exec_()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = App()
    window.show()
    app.exec_()
