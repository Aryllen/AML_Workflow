import os
import glob
import sys
import functools
import jsonpickle
from collections import OrderedDict
from Orange.widgets import widget, gui, settings
import Orange.data
from Orange.data.io import FileFormat
from DockerClient import DockerClient
from BwBase import OWBwBWidget, ConnectionDict, BwbGuiElements, getIconName, getJsonName
from PyQt5 import QtWidgets, QtGui

class OWVar_Cleaning(OWBwBWidget):
    name = "Var_Cleaning"
    description = "Verification of clean variant data"
    priority = 1
    icon = getIconName(__file__,"dtoxs-analysis2.svg")
    want_main_area = False
    docker_image_name = "aml-clean-gather"
    docker_image_tag = "latest"
    inputs = [("varDir",str,"handleInputsvarDir")]
    outputs = [("varDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    varDir=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Var_Cleaning")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsvarDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("varDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"varDir"):
            outputValue=getattr(self,"varDir")
        self.send("varDir", outputValue)
