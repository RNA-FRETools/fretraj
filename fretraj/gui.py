import sys
import os
from pymol import cmd
from pymol.Qt import QtWidgets, utils, QtCore
import mdtraj as md
import json
import copy
import webbrowser
import re
import functools

try: 
    from fretraj import cloud
    from fretraj import isosurf
    from fretraj import export
except ModuleNotFoundError: # for PyMOL Plugin
    from . import cloud
    from . import isosurf
    from . import export

try:
    import LabelLib as ll
except ModuleNotFoundError:
    _LabelLib_found = False
else:
    _LabelLib_found = True

# class App(QtWidgets.QWidget):
package_directory = os.path.dirname(os.path.abspath(cloud.__file__))

about = {}
with open(os.path.join(package_directory, '__about__.py')) as a:
    exec(a.read(), about)


class App(QtWidgets.QMainWindow):

    def __init__(self, _pymol_running=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.uiIconPath = '{}/icon'.format(package_directory)
        self.exampleDataPath = '{}/data/'.format(package_directory)
        fretrajUI = os.path.join(os.path.dirname(__file__), 'fretraj.ui')
        self.settingsUI = os.path.join(os.path.dirname(__file__), 'settings.ui')
        self.textUI = os.path.join(os.path.dirname(__file__), 'textpad.ui')
        utils.loadUi(fretrajUI, self)
        self.setWindowTitle("FRETraj")
        self.setWindowIcon(utils.QtGui.QIcon(self.uiIconPath))
        self._pymol_running = _pymol_running
        self.statusBar().showMessage("Ready", 2000)
        self.readTheDocsURL = None
        self.settingsWindow = QtWidgets.QDialog(self)
        utils.loadUi(self.settingsUI, self.settingsWindow)
        self.settingsWindow.setWindowTitle("FRETraj - Settings")
        self.textWindow = QtWidgets.QDialog(self)
        utils.loadUi(self.textUI, self.textWindow)
        

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
        self.labels_default = copy.deepcopy(self.labels)

        # set root path
        root_err = 'No root directory specified. FRETraj is not initialized.'
        self.settings = {'root_path': None, 'browser': None, 'local_docs': None}
        if os.path.isfile('{}/.fretraj_settings.conf'.format(package_directory)):
            with open('{}/.fretraj_settings.conf'.format(package_directory), 'r') as f:
                self.settings = json.load(f)
                if os.path.isdir(self.settings['root_path']):
                    self.lineEdit_rootDirectory.setText(self.settings['root_path'])
                    self.settingsWindow.lineEdit_rootDirectory.setText(self.settings['root_path'])
                else:
                    self.settings['root_path'] = None
                    self.setRootDirectory()
                    if not self.settings['root_path']:
                        raise ValueError(root_err)
        else:
            self.setRootDirectory()
            if not self.settings['root_path']:
                raise ValueError(root_err)

        # examples
        n_examples = 0
        fileformat = '.pdb'
        self.exampleQAction = {}
        files = [f for f in os.listdir(self.exampleDataPath) if f.endswith(fileformat)]
        for example in files:
            self.exampleQAction[example] = self.menuLoad_Example.addAction(example.replace(fileformat, ''))
            self.exampleQAction[example].triggered.connect(functools.partial(self.openExample, self.exampleQAction[example].text(), fileformat))
            n_examples += 1
        if n_examples > 0:
            self.menuLoad_Example.removeAction(self.example_placeholder)

        # activate / deactivate GUI elements
        self.dyeRadius_spinBox_OnOff()
        self.push_computeACV.setEnabled(False)
        self.doubleSpinBox_contourValue.setEnabled(False)
        self.doubleSpinBox_contourValue_CV.setEnabled(False)
        self.spinBox_bfactor.setEnabled(False)
        self.spinBox_gaussRes.setEnabled(False)
        self.doubleSpinBox_gridBuffer.setEnabled(False)
        self.push_calculateFRET.setEnabled(False)
        self.spinBox_statePDB.setEnabled(False)
        self.spinBox_atomID.setEnabled(False)
        self.push_transfer.setEnabled(False)
        self.push_loadParameterFile.setEnabled(False)
        self.push_showText.setEnabled(False)
        if not _LabelLib_found:
            self.checkBox_useLabelLib.setEnabled(False)
            self.checkBox_useLabelLib.setChecked(False)

        # signals
        self.push_computeACV.clicked.connect(self.computeACV)
        self.push_loadPDB.clicked.connect(self.loadPDB)
        self.comboBox_simType.activated.connect(self.dyeRadius_spinBox_OnOff)
        self.push_calculateFRET.clicked.connect(self.calculateFRET)
        self.push_deleteLabel.clicked.connect(self.deleteLabel)
        self.push_deleteFRETparam.clicked.connect(self.deleteFRETparam)
        self.push_loadParameterFile.clicked.connect(self.loadParameterFile)
        self.doubleSpinBox_contourValue.valueChanged.connect(self.pymol_update_isosurface)
        self.doubleSpinBox_contourValue_CV.valueChanged.connect(self.pymol_update_isosurface)
        self.spinBox_bfactor.valueChanged.connect(self.pymol_update_isosurface)
        self.spinBox_gaussRes.valueChanged.connect(self.pymol_update_isosurface)
        self.doubleSpinBox_gridBuffer.valueChanged.connect(self.pymol_update_isosurface)
        self.lineEdit_labelName.returnPressed.connect(self.makeLabel)
        self.lineEdit_distanceName.returnPressed.connect(self.makeFRETparam)
        self.comboBox_labelName.currentIndexChanged.connect(self.update_comboBox)
        self.actionAbout_FRETraj.triggered.connect(self.openAbout)
        self.comboBox_donorName.currentIndexChanged.connect(self.define_DA)
        self.comboBox_acceptorName.currentIndexChanged.connect(self.define_DA)
        self.comboBox_distanceName.currentIndexChanged.connect(self.update_comboBox)
        self.spinBox_statePDB.valueChanged.connect(self.update_PDBstate)
        self.push_setRootDirectory.clicked.connect(self.setRootDirectory)
        self.actionDocumentation.triggered.connect(self.openDocumentation)
        self.spinBox_atomID.valueChanged.connect(self.update_atom)
        self.push_transfer.clicked.connect(self.transferToLabel)
        self.push_deleteFRET.clicked.connect(self.deleteFRET)
        self.actionSettings.triggered.connect(self.openSettings)
        self.push_showText.clicked.connect(self.openPDBFile)


    def update_labelDict(self, pos=None):
        """
        Update the label dictionary with the values from the GUI fields
        """
        if pos is None:
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
        self.labels['Position'][pos]['state'] = self.spinBox_statePDB.value()
        self.labels['Position'][pos]['frame_mdtraj'] = self.spinBox_statePDB.value() - 1
        self.labels['Position'][pos]['use_LabelLib'] = self.checkBox_useLabelLib.isChecked()
        self.labels['Distance'][dis]['R0'] = self.doubleSpinBox_R0.value()
        self.labels['Distance'][dis]['n_dist'] = self.spinBox_nDist.value()
        self.donorName = self.comboBox_donorName.currentText()
        self.acceptorName = self.comboBox_acceptorName.currentText()
        self.labels['Position'][pos]['contour_level_AV'] = self.doubleSpinBox_contourValue.value()
        self.labels['Position'][pos]['contour_level_CV'] = self.doubleSpinBox_contourValue_CV.value()
        self.labels['Position'][pos]['b_factor'] = self.spinBox_bfactor.value()
        self.labels['Position'][pos]['gaussian_resolution'] = self.spinBox_gaussRes.value()
        self.labels['Position'][pos]['grid_buffer'] = self.doubleSpinBox_gridBuffer.value()


    def update_comboBox(self):
        self.update_labelDict()
        self.update_GUIfields()

    def update_GUIfields(self):
        """
        Update the fields of the GUI upon changing the label in the dropdown
        """
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
        self.spinBox_statePDB.setValue(self.labels['Position'][pos]['state'])
        self.checkBox_useLabelLib.setChecked(self.labels['Position'][pos]['use_LabelLib'])
        self.doubleSpinBox_R0.setValue(self.labels['Distance'][dis]['R0'])
        self.spinBox_nDist.setValue(self.labels['Distance'][dis]['n_dist'])
        # self.update_atom()
        self.doubleSpinBox_contourValue.setValue(self.labels['Position'][pos]['contour_level_AV'])
        self.doubleSpinBox_contourValue_CV.setValue(self.labels['Position'][pos]['contour_level_CV'])
        self.spinBox_bfactor.setValue(self.labels['Position'][pos]['b_factor'])
        self.spinBox_gaussRes.setValue(self.labels['Position'][pos]['gaussian_resolution'])
        self.doubleSpinBox_gridBuffer.setValue(self.labels['Position'][pos]['grid_buffer'])

        
    def loadPDB(self, fileNamePath_pdb=False):
        """
        Load PDB of CIF file. CIF files will be converted to PDB internally
        """
        if not fileNamePath_pdb:
            self.fileNamePath_pdb, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load PDB / CIF', '', "PDB / CIF file (*.pdb *cif);;All Files (*)")
        else:
            self.fileNamePath_pdb = fileNamePath_pdb
        if self.fileNamePath_pdb:
            self.fileName_pdb = self.fileNamePath_pdb.split("/")[-1]
            if self.fileName_pdb[-3:] == 'cif':
                self.fileNamePath_pdb, _ = self.cif2pdb(self.fileNamePath_pdb)
            try:
                self.struct = md.load_pdb(self.fileNamePath_pdb)
            except IndexError:
                self.openErrorWin("PDB File Error", "The specified file \"{}.pdb\" cannot be loaded".format(name))
                return 0
            else:
                self.lineEdit_pdbFile.setText(self.fileName_pdb)
                self.spinBox_statePDB.setMaximum(self.struct.n_frames)
                self.spinBox_atomID.setMaximum(self.struct.n_atoms)
                self.update_atom()
                self.push_computeACV.setEnabled(True)
                self.spinBox_statePDB.setEnabled(True)
                self.spinBox_atomID.setEnabled(True)
                self.push_transfer.setEnabled(True)
                self.push_loadParameterFile.setEnabled(True)
                with open(self.fileNamePath_pdb, 'r') as f:
                    self.pdbText = f.read()
                self.push_showText.setEnabled(True)
                if self._pymol_running:
                    cmd.reinitialize()
                    cmd.load(self.fileNamePath_pdb)
                return 1

    def cif2pdb(self, pathtoPDBfile):
        """
        Convert .cif to .pdb with PyMOL

        Returns
        -------
        fileNamePath_pdb : string
        fileName_pdb : string
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
        if self._pymol_running:
            self.update_labelDict()
            cmd.set('state', self.labels['Position'][self.labelName]['state'])

    def update_atom(self):
        self.lineEdit_pdbAtom.setText(str(self.struct.top.atom(self.spinBox_atomID.value() - 1)))

    def loadParameterFile(self, fileNamePath_param=False):
        """
        Load dye, linker and simulation settings from parameter file
        """
        if not fileNamePath_param:
            self.fileNamePath_param, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load label parameter file', '', "JSON Files (*.json);;All Files (*)")
        else:
            self.fileNamePath_param = fileNamePath_param
        if self.fileNamePath_param:
            self.fileName_param = self.fileNamePath_param.split("/")[-1]
            try:
                with open(self.fileNamePath_param) as f:
                    try:
                        labels_json = json.load(f)
                    except json.decoder.JSONDecodeError:
                        self.openErrorWin('Parameter File Error', "The specified file \"{}\" is not of type JSON".format(self.fileName_param))
                        return 0
                    else:
                        self.lineEdit_paramFile.setText(self.fileName_param)
                        error_msg = "The specified file \"{}\" has a wrong format".format(self.fileName_param)
                        for field in labels_json.keys():
                            if field == 'Distance':
                                name_default = self.distanceName_default
                            else:
                                name_default = self.labelName_default

                            if not isinstance(labels_json[field], dict):
                                self.openErrorWin('Parameter File Error', error_msg)
                                return 0

                            #print(self.labels_default)
                            for pos_dis in labels_json[field].keys():
                                
                                #if pos_dis != name_default:

                                self.labels[field][pos_dis] = {}
                                #print(self.labels_default)
                                #print(self.labels_default[field][name_default].keys())
                                
                                for key in self.labels_default[field][name_default].keys():
                                    if isinstance(labels_json[field][pos_dis], dict):
                                        if key in labels_json[field][pos_dis].keys():
                                            self.labels[field][pos_dis][key] = labels_json[field][pos_dis][key]

                                        else:
                                            self.labels[field][pos_dis][key] = copy.deepcopy(self.labels_default[field][name_default][key])
                                    else:
                                        self.openErrorWin('Parameter File Error', error_msg)
                                        return 0
                    for newlabel in self.labels['Position'].keys():
                        #self.update_GUIfields()
                        self.addLabel(newlabel)
                    for newdistance in self.labels['Distance'].keys():
                        self.addFRETparam(newdistance)
                    return 1
            except FileNotFoundError:
                self.openErrorWin('Parameter File Error', "The specified file \"{}\" could not be found".format(self.fileName_param))
                return 0

    def transferToLabel(self):
        newlabel = '{:d}-{}'.format(self.spinBox_statePDB.value(), self.lineEdit_pdbAtom.text())
        self.labels['Position'][newlabel] = copy.deepcopy(self.labels['Position'][self.labelName])
        self.update_labelDict(newlabel)
        self.update_GUIfields()
        self.addLabel(newlabel)

    def makeLabel(self):
        """
        Create new label from edit box string
        """
        newlabel = self.lineEdit_labelName.text()
        self.labels['Position'][newlabel] = copy.deepcopy(self.labels['Position'][self.labelName])
        self.update_labelDict(newlabel)
        self.update_GUIfields()
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
        self.update_GUIfields()
        #self.update_labelDict()

    def deleteLabel(self):
        """
        Remove the selected label from Combobox
        """
        self.deleteLabelFromList()
        self.deleteDistanceFromList()
        self.deleteIsosurface()
        if self.labelName != self.labelName_default:
            self.comboBox_labelName.removeItem(self.comboBox_labelName.currentIndex())
            self.lineEdit_labelName.clear()

    def addLabelToList(self, av):
        """
        Add a label and the computed mean dye position (MP) to the list
        """
        rowCount = self.tableWidget_MeanPos.rowCount()
        for r in range(rowCount):
            if self.labelName == self.tableWidget_MeanPos.item(r, 0).text():
                break
        else:
            r = 0
            self.tableWidget_MeanPos.insertRow(r)
        self.tableWidget_MeanPos.setItem(r, 0, QtWidgets.QTableWidgetItem(self.labelName))
        self.tableWidget_MeanPos.setItem(r, 1, QtWidgets.QTableWidgetItem(', '.join('{:.2f}'.format(p) for p in av.acv.mp)))
        if self.comboBox_donorName.findText(self.labelName) == -1:
            self.comboBox_donorName.addItem(self.labelName)
        if self.comboBox_acceptorName.findText(self.labelName) == -1:
            self.comboBox_acceptorName.addItem(self.labelName)

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

    def deleteIsosurface(self):
        av_name = self.labelName.replace('\'', 'p')
        cmd.delete(av_name)
        cmd.delete(av_name + '_map')
        cmd.delete(av_name + '_isosurf')
        cmd.delete(av_name + '_CV_map')
        cmd.delete(av_name + '_CV_isosurf')

    def makeFRETparam(self):
        """
        Create new FRET parameters from edit box string
        """
        newdistance = self.lineEdit_distanceName.text()
        self.labels['Distance'][newdistance] = copy.deepcopy(self.labels_default['Distance'][self.distanceName_default])
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
        Define the donor and acceptor labels and color them in the mean position and in the FRET table
        """
        self.update_labelDict()
        for i in range(self.tableWidget_MeanPos.rowCount()):
            for j in range(self.tableWidget_MeanPos.columnCount()):
                self.tableWidget_MeanPos.item(i, j).setBackground(utils.QtGui.QColor(255, 255, 255))

        for i in range(self.tableWidget_FRET.rowCount()):
            for j in range(self.tableWidget_FRET.columnCount()):
                self.tableWidget_FRET.item(i, j).setBackground(utils.QtGui.QColor(255, 255, 255))

        if self.tableWidget_MeanPos.rowCount() > 1:
            if self.donorName == self.acceptorName:
                self.push_calculateFRET.setEnabled(False)
                if self._pymol_running:
                    cmd.color('white', '*_isosurf')
            else:
                self.push_calculateFRET.setEnabled(True)
                row_don = self.tableWidget_MeanPos.findItems(self.donorName, QtCore.Qt.MatchExactly)[0].row()
                row_acc = self.tableWidget_MeanPos.findItems(self.acceptorName, QtCore.Qt.MatchExactly)[0].row()
                for j in range(self.tableWidget_MeanPos.columnCount()):
                    self.tableWidget_MeanPos.item(row_don, j).setBackground(utils.QtGui.QColor(108, 179, 129))
                    self.tableWidget_MeanPos.item(row_acc, j).setBackground(utils.QtGui.QColor(194, 84, 73))

                DA = '{} -> {}'.format(self.donorName, self.acceptorName)
                for r in range(self.tableWidget_FRET.rowCount()):
                    if DA in self.tableWidget_FRET.item(r, 0).text():
                        for j in range(self.tableWidget_FRET.columnCount()):
                            self.tableWidget_FRET.item(r, j).setBackground(utils.QtGui.QColor(200, 200, 200))
                        break

                if self._pymol_running:
                    don_name = self.donorName.replace('\'', 'p')
                    acc_name = self.acceptorName.replace('\'', 'p')
                    cmd.set_color('don_green', [108, 179, 129])
                    cmd.set_color('acc_red', [194, 84, 73])
                    cmd.color('don_green', don_name + '_isosurf')
                    cmd.color('don_green', don_name + '_CV_isosurf')
                    cmd.color('acc_red', acc_name + '_isosurf')
                    cmd.color('acc_red', acc_name + '_CV_isosurf')

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
        msg = 'Busy...'
        self.statusBar().showMessage(msg, 3000)
        param_filename = '{}/parameters.json'.format(self.settings['root_path'])
        export.save_labels(param_filename, self.labels)
        self.av[self.labelName] = cloud.Volume(self.struct, self.labelName, self.labels)
        if self.av[self.labelName].acv is None:
            msg = 'ACV could not be calculated!'
            self.statusBar().showMessage(msg, 3000)
            print(msg)
        else:
            av_name = self.labelName.replace('\'', 'p')
            av_filename = '{}/{}.pdb'.format(self.settings['root_path'], av_name)
            self.av[self.labelName].save_acv(av_filename, format='pdb')

            self.addLabelToList(self.av[self.labelName])
            self.define_DA()
            msg = 'ACV successfully calculated!'
            self.statusBar().showMessage(msg, 3000)
            if self._pymol_running:
                cmd.delete('{}*'.format(av_name))
                cmd.load(av_filename)
                self.doubleSpinBox_contourValue.setEnabled(True)
                self.spinBox_bfactor.setEnabled(True)
                self.spinBox_gaussRes.setEnabled(True)
                self.doubleSpinBox_gridBuffer.setEnabled(True)
                self.pymol_update_isosurface()

    def calculateFRET(self):
        self.update_labelDict()
        self.traj[(self.donorName, self.acceptorName)] = cloud.FRET_Trajectory(self.av[self.donorName], self.av[self.acceptorName], self.distanceName, self.labels)
        self.addDistanceToList(self.traj[(self.donorName, self.acceptorName)])
        self.define_DA()
        self.traj[(self.donorName, self.acceptorName)].save_fret('{}/{}-{}.json'.format(self.settings['root_path'], self.donorName, self.acceptorName))

    def deleteFRET(self):
        DA = '{} -> {}'.format(self.donorName, self.acceptorName)
        r = 0
        rowCount = self.tableWidget_FRET.rowCount()
        while r < rowCount:
            if DA in self.tableWidget_FRET.item(r, 0).text():
                self.tableWidget_FRET.removeRow(r)
                rowCount -= 1
            else:
                r += 1

    def pymol_update_isosurface(self):
        """
        Update the isosurface when spin box value is changed.
        The spin box is only active when PYMOL is running.
        """
        if self.av:
            av_name = self.labelName.replace('\'', 'p')
            contour_level = self.doubleSpinBox_contourValue.value()
            bfactor = self.spinBox_bfactor.value()
            gaussRes = self.spinBox_gaussRes.value()
            gridBuffer = self.doubleSpinBox_gridBuffer.value()
            grid_spacing = self.labels['Position'][self.labelName]['grid_spacing']
            isosurf.smooth_map_from_xyz(av_name, av_name, contour_level, grid_spacing, bfactor, gaussRes, gridBuffer)
            if any(self.av[self.labelName].acv.tag_1d > 1):
                self.doubleSpinBox_contourValue_CV.setEnabled(True)
                contour_level_CV = self.doubleSpinBox_contourValue_CV.value()
                sele_CV = '{} and resn CV'.format(av_name)
                isosurf.smooth_map_from_xyz(av_name + '_CV', sele_CV, contour_level_CV, grid_spacing, bfactor, gaussRes, gridBuffer)
                cmd.set('transparency', 0.4, av_name + '_isosurf')
            else:
                self.doubleSpinBox_contourValue_CV.setEnabled(False)
                cmd.set('transparency', 0, av_name + '_isosurf')

    def setRootDirectory(self):
        rootDir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Set root directory')
        if rootDir:
            self.lineEdit_rootDirectory.setText(rootDir)
            self.settingsWindow.lineEdit_rootDirectory.setText(self.settings['root_path'])
            self.settings['root_path'] = rootDir
            with open('{}/.fretraj_settings.conf'.format(package_directory), 'w') as f:
                json.dump(self.settings, f)

    def openDocumentation(self):
        if not self.readTheDocsURL:
            if not self.settings['local_docs']:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Information)
                msg.setWindowTitle("Location of docs not configured")
                msg.setStandardButtons(QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel)
                msg.setText('Press <OK> and specify the path of the local docs.')
                returnValue = msg.exec_()
                if returnValue == QtWidgets.QMessageBox.Ok:
                    docs_path = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select docs directory')
                    self.settings['local_docs'] = re.escape(docs_path)
                    with open('{}/.fretraj_settings.conf'.format(package_directory), 'w') as f:
                        json.dump(self.settings, f)

        if not self.settings['browser']:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.setWindowTitle("Web browser not configured")
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel)
            msg.setText('For the documentation to be displayed, the path to a web browser needs to be configured. Press <OK> and search for the browser executable (Firefox, Chrome or Edge).')
            returnValue = msg.exec_()
            if returnValue == QtWidgets.QMessageBox.Ok:
                browser_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Search for default browser')
                self.settings['browser'] = re.escape(browser_path)
                with open('{}/.fretraj_settings.conf'.format(package_directory), 'w') as f:
                    json.dump(self.settings, f)

        if self.settings['browser']:
            try:
                browser = webbrowser.get('{} %s'.format(self.settings['browser']))
                browser.open('file://{}/index.html'.format(self.settings['local_docs']))
            except webbrowser.Error:
                self.settings['browser'] = None
                print('Browser not found!')
                with open('{}/.fretraj_settings.conf'.format(package_directory), 'w') as f:
                    json.dump(self.settings, f)
            except FileNotFoundError:
                self.settings['local_docs'] = None
                print('Local docs not found!')
                with open('{}/.fretraj_settings.conf'.format(package_directory), 'w') as f:
                    json.dump(self.settings, f)

    def openAbout(self):
        """
        Open About Window
        """
        msg = QtWidgets.QMessageBox()
        pixmap = utils.QtGui.QPixmap(self.uiIconPath)
        msg.setIconPixmap(pixmap.scaledToWidth(64))
        msg.setWindowTitle("About FRETraj")
        msg.setText('{} {}\n\n{}\n\n(C) {}'.format(about['__title__'], str(about['__version__']), about['__description__'], about['__copyright__']))
        msg.exec_()

    def openExample(self, name, fileformat):
        """
        Load an example file and calculate an ACV
        """
        fileNamePath_pdb = '{}/{}{}'.format(self.exampleDataPath, name, fileformat)
        fileNamePath_param = '{}/{}_labels{}'.format(self.exampleDataPath, name, '.json')
        pdb_load = self.loadPDB(fileNamePath_pdb)
        if pdb_load:
            param_load = self.loadParameterFile(fileNamePath_param)
            if param_load:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Information)
                msg.setWindowTitle("Calculate ACV?")
                msg.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
                msg.setText('Would you like to calculate an ACV?')
                returnValue = msg.exec_()
                if returnValue == QtWidgets.QMessageBox.Yes:
                    n_pos = 0
                    for pos in self.labels['Position'].keys():
                        if pos != self.labelName_default:
                            index = self.comboBox_labelName.findText(pos, QtCore.Qt.MatchExactly)
                            self.comboBox_labelName.setCurrentIndex(index)
                            self.computeACV()
                            n_pos += 1
                    if n_pos > 1:
                        self.comboBox_donorName.setStyleSheet('QComboBox {background-color: rgb(229,134,145);}')
                        self.comboBox_acceptorName.setStyleSheet('QComboBox {background-color: rgb(229,134,145);}')
                        msg = QtWidgets.QMessageBox()
                        msg.setIcon(QtWidgets.QMessageBox.Information)
                        msg.setWindowTitle("Calculate FRET?")
                        msg.setText('To calculate FRET select the donor and acceptor position from the highlighted dropdown menu and press on <Calculate FRET>')
                        returnValue = msg.exec_()
                        self.comboBox_donorName.setStyleSheet('')
                        self.comboBox_acceptorName.setStyleSheet('')

    def openSettings(self):
        self.settingsWindow.lineEdit_rootDirectory.setText(self.settings['root_path'])
        self.settingsWindow.lineEdit_browser.setText(self.settings['browser'])
        self.settingsWindow.lineEdit_localdocs.setText(self.settings['local_docs'])
        isOK = self.settingsWindow.exec_()
        if isOK:
            self.settings['root_path'] = self.settingsWindow.lineEdit_rootDirectory.text()
            self.settings['browser'] = self.settingsWindow.lineEdit_browser.text()
            self.settings['local_docs'] = self.settingsWindow.lineEdit_localdocs.text()
            self.lineEdit_rootDirectory.setText(self.settings['root_path'])
            with open('{}/.fretraj_settings.conf'.format(package_directory), 'w') as f:
                json.dump(self.settings, f)

    def openPDBFile(self):
        self.textWindow.textBrowser_pdbFile.setText(self.pdbText)
        self.textWindow.setWindowTitle("FRETraj - {}".format(self.fileName_pdb))
        isOK = self.textWindow.exec_()

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
