# -*- coding: utf-8 -*-
from PyQt5 import QtWidgets, QtGui, QtCore
from project_tab_ui import Ui_MyProjectTab
from comp_entry_ui import Ui_CompEntry
from hapi import ISO_ID
import os
import json


class CompEntry(QtWidgets.QWidget):
    _counter = 0
    # closed = QtCore.pyqtSignal()

    def __init__(self):
        super(CompEntry, self).__init__()
        type(self)._counter += 1
        self.ui = Ui_CompEntry()
        self.ui.setupUi(self)
        self.molec = self._setComboBoxes()
        self.ui.molec_val.addItems(self.molec)

        # Set Validators
        dval = QtGui.QDoubleValidator()
        self.ui.amount_val.setValidator(dval)

        # Connections
        self.ui.toolButton.clicked.connect(self.onRemove)
        self.ui.molec_val.activated.connect(self.setIso_val)

    def __del__(self):
        type(self)._counter -= 1

    def getValues(self):
        # Convert amount to
        return self.ui.molec_val.currentText(),\
               self.ui.iso_val.currentText(),\
               float(self.ui.amount_val.displayText()),\
               self.ui.amount_unit.currentText()

    def onRemove(self):
        # self.closed.emit()
        self.close()

    def setIso_val(self):
        molec = self.ui.molec_val.currentText()
        ans = ["Natural"]
        for k, v in ISO_ID.items():
            if v[-1] == molec:
                ans.append(v[2])
        self.ui.iso_val.addItems(ans)

    def _setComboBoxes(self):
        ans = []
        for k, v in ISO_ID.items():
            if v[-1] not in ans:
                ans.append(v[-1])
        return ans

    def setCompVal(self, data):
        _uconv = {'ppm': 1e6, 'ppb': 1e9, 'ppt': 1e12, 'mb': 1.01325e3}
        self.ui.molec_val.setCurrentText(data['molec'])
        self.setIso_val()
        self.ui.iso_val.setCurrentText(data['isotop'])
        self.ui.amount_val.setText("{:.1f}".format(
            data['vmr']*_uconv[data['amountU']]))
        self.ui.amount_unit.setCurrentText(data['amountU'])


class MyProjectTab(QtWidgets.QWidget):
    filesUpdate = QtCore.pyqtSignal(list)
    fitFiles = QtCore.pyqtSignal(dict, dict)
    # fitFiles = QtCore.pyqtSignal()

    def __init__(self):
        super(MyProjectTab, self).__init__()
        self.ui = Ui_MyProjectTab()
        self.ui.setupUi(self)
        self.projectName = None
        self.projectDir = None
        self.projectDB = None
        self.files = None
        self.cwd = os.getcwd()
        self.components = []
        self.components.append(CompEntry())
        self.ui.Component_Layout.addWidget(self.components[-1])
        # Set Validators
        dval = QtGui.QDoubleValidator()
        self.ui.temp_val.setValidator(dval)
        self.ui.press_val.setValidator(dval)
        self.ui.pathL_val.setValidator(dval)
        self.ui.roiLow_val.setValidator(dval)
        self.ui.roiHigh_val.setValidator(dval)
        self.ui.fov_val.setValidator(dval)
        self.ui.res_val.setValidator(dval)

        # Connections
        self.ui.addComponentButton.clicked.connect(self.on_AddComponent)
        self.ui.fitSpectraButton.clicked.connect(self.on_FitSpectra)
        self.ui.pdir_Button.clicked.connect(self.getDir)
        self.ui.loadfiles_Button.clicked.connect(self.getFiles)
        self.ui.project_name.editingFinished.connect(self.onProjectEnter)
        self.ui.saveProj_Button.clicked.connect(self.saveProject)
        self.ui.loadProj_Button.clicked.connect(self.loadProject)
        # CompEntry.closed.connect(self.on_RemoveComponent)

    def _collectVal(self):
        # Check if there are values
        _unit_conv = {'ppm': 1e-6, 'ppb': 1e-9, 'ppt': 1e-12,
                      'mb': 0.986923266716e-3, 'm': 1e2, 'cm': 1, 'km': 1e5,
                      'K': 1, 'C': 1, 'atm': 1, 'bar': 0.986923266716}
        projectKey = ['pname', 'pdir', 'p_db', 'pcomments', 'pfiles']
        projectDict = dict(zip(projectKey, self.getProjectVal()))
        # Build sample dictonary
        sampleKey = ['temp', 'tempU', 'press', 'pressU', 'pathL', 'pathLU',
                     'ROI', 'res', 'apod', 'fov', 'ftype', 'bgfile']
        sampleDict = dict(zip(sampleKey, self.getSampleVal()))

        # Change the units to HITRAN units
        for v, u in zip(['temp', 'press', 'pathL'],
                        ['tempU', 'pressU', 'pathLU']):
            sampleDict[v] = _unit_conv[sampleDict[u]]*sampleDict[v]
        # Build components dictonaries
        compkey = ['molec', 'isotop', 'vmr', 'amountU']
        components = []
        for comp in self.components:
            if not comp.isHidden():
                components.append(dict(zip(compkey, comp.getValues())))
        # Calculate the vmr of each component
        for comp in components:
            tv = _unit_conv[comp['amountU']]*comp['vmr']
            if comp['amountU'] in ['mb']:
                tv = tv/sampleDict['press']
            comp['vmr'] = tv

        sampleDict['comp'] = components
        return projectDict, sampleDict

    def saveProject(self):
        projDict, sampleDict = self._collectVal()
        data = {'projectDict': projDict, 'sampleDict': sampleDict}
        data_json = json.dumps(data)
        ftmp = projDict['pname']+".json"
        fname = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', ftmp,
                                                      'Project;(*.json)', self.cwd)
        with open(fname[0], "w") as ofile:
            ofile.write(data_json)

    def loadProject(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', '',
                                                      'Project;(*.json)', self.cwd)
        with open(fname[0], 'r') as infile:
            data = json.load(infile)
        projDict = data['projectDict']
        sampleDict = data['sampleDict']
        self.setProjectVal(projDict)
        self.setSampleVal(sampleDict)

    def setSampleVal(self, data):
        _uconv = {'ppm': 1e6, 'ppb': 1e9, 'ppt': 1e12,
                  'mb': 1.01325e3, 'm': 1e-2, 'cm': 1, 'km': 1e-5,
                  'K': 1, 'C': 1, 'atm': 1, 'bar': 1.01325}
        self.ui.temp_val.setText(str(data['temp']))
        self.ui.temp_unit.setCurrentText(data['tempU'])
        self.ui.press_val.setText(str(data['press']*_uconv[data['pressU']]))
        self.ui.press_unit.setCurrentText(data['pressU'])
        self.ui.pathL_val.setText(str(data['pathL']*_uconv[data['pathLU']]))
        self.ui.pathL_unit.setCurrentText(data['pathLU'])
        self.ui.roiLow_val.setText(str(data['ROI'][0]))
        self.ui.roiHigh_val.setText(str(data['ROI'][1]))
        self.ui.res_val.setText(str(data['res']))
        self.ui.apod_val.setCurrentText(data['apod'])
        self.ui.fov_val.setText(str(data['fov']))
        self.ui.ftype_val.setCurrentText(data['ftype'])
        self.ui.bgfile_val.setText(data['bgfile'])

        # Build the components
        # Add components and then change them
        components = data['comp']
        for n in range(len(components)-1):
            self.components.append(CompEntry())
            self.ui.Component_Layout.addWidget(self.components[-1])
        for comp, data in zip(self.components, components):
            comp.setCompVal(data)

    def setProjectVal(self, data):
        self.projectName = data['pname']
        self.ui.project_name.setText(self.projectName)
        self.projectDir = data['pdir']
        self.ui.project_dir.setText(os.path.split(self.projectDir)[-1])
        self.projectDB = data['p_db']
        self.ui.project_DB_dir.setText(
            os.path.join(os.path.split(self.projectDir)[-1],
                         os.path.split(self.projectDB)[-1]))
        self.ui.comments_val.setPlainText(data['pcomments'])
        self.files = data['pfiles']
        if self.files:
            self.filesUpdate.emit(self.files)

    def on_AddComponent(self):
        self.components.append(CompEntry())
        self.ui.Component_Layout.addWidget(self.components[-1])

    def on_FitSpectra(self):
        self.fitFiles.emit(*self._collectVal())

    def getProjectVal(self):
        return self.projectName,\
               self.projectDir,\
               self.projectDB,\
               self.ui.comments_val.toPlainText(),\
               self.files

    def getSampleVal(self):
        return float(self.ui.temp_val.text()),\
               self.ui.temp_unit.currentText(),\
               float(self.ui.press_val.text()),\
               self.ui.press_unit.currentText(),\
               float(self.ui.pathL_val.text()),\
               self.ui.pathL_unit.currentText(),\
               (float(self.ui.roiLow_val.text()),
                float(self.ui.roiHigh_val.text())),\
               float(self.ui.res_val.text()),\
               self.ui.apod_val.currentText(),\
               float(self.ui.fov_val.text()),\
               self.ui.ftype_val.currentText(),\
               self.ui.bgfile_val.text()

    def getDir(self):
        diag = QtWidgets.QFileDialog(self)
        diag.setFileMode(QtWidgets.QFileDialog.Directory)
        diag.setOption(QtWidgets.QFileDialog.ShowDirsOnly)
        if diag.exec_():
            dirName = diag.selectedFiles()
            self.ui.project_dir.setText(dirName[0])

    def getFiles(self):
        diag = QtWidgets.QFileDialog(self)
        diag.setDirectory(self.cwd)
        diag.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        if diag.exec_():
            if not self.files:
                self.files = []
            for file in diag.selectedFiles():
                if file not in self.files:
                    self.files.append(file)
            self.filesUpdate.emit(self.files)

    def onProjectEnter(self):
        self.projectName = "{}".format(self.ui.project_name.text())
        self.projectDir = os.path.join(self.cwd, self.projectName.upper())
        self.projectDB = os.path.join(self.cwd, self.projectName+'.db')
        self.ui.project_dir.setText(os.path.split(self.projectDir)[-1])
        self.ui.project_DB_dir.setText(
            os.path.join(os.path.split(self.projectDir)[-1],
                         os.path.split(self.projectDB)[-1]))
