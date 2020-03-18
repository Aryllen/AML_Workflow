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

class OWDS_Rm_gnumeric(OWBwBWidget):
    name = "DS_Rm_gnumeric"
    description = "Open source spreadsheet"
    priority = 2
    icon = getIconName(__file__,"gnumeric.png")
    want_main_area = False
    docker_image_name = "biodepot/gnumeric"
    docker_image_tag = "1.12.18__debian-jessie-slim__080418"
    inputs = [("dsDir",str,"handleInputsdsDir")]
    outputs = [("dsDir",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    dsDir=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"DS_Rm_gnumeric")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsdsDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("dsDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"dsDir"):
            outputValue=getattr(self,"dsDir")
        self.send("dsDir", outputValue)
