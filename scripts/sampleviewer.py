#!/usr/bin/env python
import sys
from PyQt5.QtWidgets import QApplication
from sampleviewer import SampleViewer

app = QApplication(sys.argv)
filename = None
if len(sys.argv) > 1:
    filename = sys.argv[1]
window = SampleViewer(filename)
sys.exit(app.exec_())
