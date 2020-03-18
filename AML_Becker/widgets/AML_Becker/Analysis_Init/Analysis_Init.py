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

class OWAnalysis_Init(OWBwBWidget):
    name = "Analysis_Init"
    description = "Initial drug sensitivity and variant analysis."
    priority = 1
    icon = getIconName(__file__,"bioc-r.png")
    want_main_area = False
    docker_image_name = "aml-clean-gather"
    docker_image_tag = "latest"
    inputs = [("varDir",str,"handleInputsvarDir"),("dsDir",str,"handleInputsdsDir"),("outputDir",str,"handleInputsoutputDir"),("alpha",str,"handleInputsalpha")]
    outputs = [("outputDir",str),("activeMax",str),("alpha",str)]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    varDir=pset(None)
    dsDir=pset(None)
    outputDir=pset(None)
    frequency=pset(None)
    siftDel=pset(False)
    siftTol=pset(False)
    polyProb=pset(False)
    polyPoss=pset(False)
    polyBenign=pset(False)
    activeMin=pset(None)
    activeMax=pset("0.000001")
    alpha=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Analysis_Init")) as f:
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
    def handleInputsdsDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("dsDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsoutputDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleInputsalpha(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("alpha", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
    def handleOutputs(self):
        outputValue=None
        if hasattr(self,"outputDir"):
            outputValue=getattr(self,"outputDir")
        self.send("outputDir", outputValue)
        outputValue="0.000001"
        if hasattr(self,"activeMax"):
            outputValue=getattr(self,"activeMax")
        self.send("activeMax", outputValue)
        outputValue="0.05"
        if hasattr(self,"alpha"):
            outputValue=getattr(self,"alpha")
        self.send("alpha", outputValue)
