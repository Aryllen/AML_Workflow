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

class OWAnalysis_Jupyter(OWBwBWidget):
    name = "Analysis_Jupyter"
    description = "Jupyter notebook installation with packages for AML analysis"
    priority = 103
    icon = getIconName(__file__,"jupyter_image.png")
    want_main_area = False
    docker_image_name = "aml-analysis-jupyter"
    docker_image_tag = "latest"
    inputs = [("outputDir",str,"handleInputsoutputDir")]
    pset=functools.partial(settings.Setting,schema_only=True)
    runMode=pset(0)
    exportGraphics=pset(False)
    runTriggers=pset([])
    triggerReady=pset({})
    inputConnectionsStore=pset({})
    optionsChecked=pset({})
    startingNotebook=pset(None)
    type=pset("notebook")
    outputNotebook=pset(None)
    debug=pset(False)
    generateConfig=pset(False)
    autoyes=pset(True)
    allowRoot=pset(True)
    loglevel=pset("30")
    port=pset(8888)
    config=pset(None)
    transport=pset(None)
    keyfile=pset(None)
    certfile=pset(None)
    clientca=pset(None)
    outputDir=pset(None)
    browser=pset(None)
    def __init__(self):
        super().__init__(self.docker_image_name, self.docker_image_tag)
        with open(getJsonName(__file__,"Analysis_Jupyter")) as f:
            self.data=jsonpickle.decode(f.read())
            f.close()
        self.initVolumes()
        self.inputConnections = ConnectionDict(self.inputConnectionsStore)
        self.drawGUI()
    def handleInputsoutputDir(self, value, *args):
        if args and len(args) > 0: 
            self.handleInputs("outputDir", value, args[0][0], test=args[0][3])
        else:
            self.handleInputs("inputFile", value, None, False)
