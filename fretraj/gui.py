import sys
import os
from pymol import cmd
from pymol.Qt import QtWidgets, utils, QtCore
import mdtraj as md

from fretraj import cloud
from fretraj import isosurf


class App(QtWidgets.QWidget):

    def __init__(self, _pymol_running=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        uifile = os.path.join(os.path.dirname(__file__), 'fretraj.ui')
        uiIcon = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'hplc2spect_icon.png')
        utils.loadUi(uifile, self)
        self.setWindowIcon(utils.QtGui.QIcon(uiIcon))

        self._pymol_running = _pymol_running

        self.labels = {'Position': {}, 'Distance': {}}
        self.simulationType()
        self.push_computeACV.setEnabled(False)

        self.push_computeACV.clicked.connect(self.computeACV)
        self.push_loadPDB.clicked.connect(self.loadPDB)
        self.comboBox_simType.activated.connect(self.simulationType)
        self.push_calculateFRET.clicked.connect(self.calculateFRET)
        self.push_deleteLabel.clicked.connect(self.deleteLabelFromList)

    def update_label_dict(self):

        pos = self.lineEdit_labelName.text()
        self.labelName = pos
        self.labels['Position'][pos] = {}
        self.labels['Distance'][pos] = {}
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
        self.labels['Position'][pos]['dye_radius3'] = self.doubleSpinBox_dyeRadius3.value()
        self.labels['Position'][pos]['mol_selection'] = self.lineEdit_molSelection.text()
        self.labels['Distance'][pos]['R0'] = self.doubleSpinBox_R0.value()

        cloud.check_labels(self.labels)

    def loadPDB(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        self.fileNamePath, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load PDB', '', "PDB Files (*.pdb);;All Files (*)", options=options)
        if self.fileNamePath:
            self.in_filePDB = self.fileNamePath.split("/")[-1]

            try:
                self.struct = md.load_pdb(self.fileNamePath)
                if self._pymol_running:
                    cmd.load(self.fileNamePath)

            except IndexError:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setWindowTitle("PDB File Error")
                msg.setText("The specified file \"{}\" can be not be loaded".format(self.in_filePDB))
                msg.exec_()
            else:
                self.lineEdit_pdbFile.setText(self.in_filePDB)
                self.push_computeACV.setEnabled(True)

    def simulationType(self):
        if self.comboBox_simType.currentText() == 'AV1':
            self.doubleSpinBox_dyeRadius2.setEnabled(False)
            self.doubleSpinBox_dyeRadius3.setEnabled(False)
        else:
            self.doubleSpinBox_dyeRadius2.setEnabled(True)
            self.doubleSpinBox_dyeRadius3.setEnabled(True)

    def computeACV(self):
        """
        Create a cloud.Volume object at the specified label position. If PyMOL is running the ACV is rendered in the current window
        """
        self.update_label_dict()

        av = cloud.Volume(self.struct, 0, self.labelName, self.labels)
        if av.acv is None:
            print('-> ACV could not be calculated')
        else:
            av_filename = '{}.xyz'.format(self.labelName)
            av.save_acv(av_filename, format='xyz')
            self.addLabelToList(av)
            if self._pymol_running:
                cmd.load(av_filename)
                isoval = self.spinBox_atomID.value()
                isosurf.smooth_map_from_xyz(self.labelName, self.labelName, isoval)

    def addLabelToList(self, av):
        self.deleteLabelFromList()
        rowCount = self.tableWidget_MeanPos.rowCount()
        self.tableWidget_MeanPos.insertRow(rowCount)
        self.tableWidget_MeanPos.setItem(rowCount, 0, QtWidgets.QTableWidgetItem(self.labelName))
        self.tableWidget_MeanPos.setItem(rowCount, 1, QtWidgets.QTableWidgetItem(', '.join('{:.2f}'.format(p) for p in av.acv.mp)))
        #print(self.tableWidget_MeanPos.findItems(self.labelName, QtCore.Qt.MatchExactly))
        #print(self.tableWidget_MeanPos.findItems('label_2', QtCore.Qt.MatchExactly))

    def deleteLabelFromList(self):
        self.labelName = self.lineEdit_labelName.text()
        rowCount = self.tableWidget_MeanPos.rowCount()
        for r in range(rowCount):
            if self.labelName == self.tableWidget_MeanPos.item(r, 0).text():
                self.tableWidget_MeanPos.removeRow(r)
                break

    def calculateFRET(self):
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setWindowTitle("FRET calculation")
        msg.setText("FRET is cool!")
        msg.exec_()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = App()
    window.show()
    app.exec_()
