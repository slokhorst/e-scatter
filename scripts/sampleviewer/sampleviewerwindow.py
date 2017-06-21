from PyQt5.QtWidgets import (QSizePolicy, QWidget, QMainWindow, QAction, QFileDialog, QCheckBox, QPushButton, QHBoxLayout, QVBoxLayout, QLabel)
from .sampleviewerwidget import SampleViewerWidget


class SampleViewer(QMainWindow):
    def __init__(self, initial_geometry_fn=None):
        super().__init__()

        self.initUI()

        self.show()

        if initial_geometry_fn is not None:
            self.loadFile(initial_geometry_fn)

    def loadFile(self, filename):
        self.viewer.loadFile(filename)
        self.dimension_lbl.setText('size: {}'.format(self.viewer.max_dim))
        self.triangles_lbl.setText('{} triangles\n(rendering {})'.format(int(self.viewer.num_triangles/2), self.viewer.num_triangles))

    def initUI(self):
        self.setWindowTitle('Sample viewer')
        self.resize(1280, 720)

        self.statusBar()

        # menubar
        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)

        openAction = QAction('&Open', self)
        openAction.setShortcut('Ctrl+O')
        openAction.setStatusTip('Load sample')
        openAction.triggered.connect(self.buttonHandler)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openAction)
        fileMenu.addAction(exitAction)

        # viewer
        self.viewer = SampleViewerWidget()
        self.viewer.setSizePolicy(
            QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred))

        # sidebar
        self.camreset_btn = QPushButton('reset camera', self)
        self.camreset_btn.clicked.connect(self.buttonHandler)
        self.wireframe_cb = QCheckBox('wireframe', self)
        self.wireframe_cb.clicked.connect(self.buttonHandler)
        self.lighting_cb = QCheckBox('lighting', self)
        self.lighting_cb.setCheckState(2)
        self.lighting_cb.clicked.connect(self.buttonHandler)
        self.dimension_lbl = QLabel('', self)
        self.triangles_lbl = QLabel('', self)

        self.sidebar = QWidget(self)
        self.sidebar.setSizePolicy(
            QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Preferred))
        sidebar_layout = QVBoxLayout(self.sidebar)
        sidebar_layout.addWidget(self.camreset_btn)
        sidebar_layout.addWidget(self.wireframe_cb)
        sidebar_layout.addWidget(self.lighting_cb)
        sidebar_layout.addStretch(1)
        sidebar_layout.addWidget(self.dimension_lbl)
        sidebar_layout.addWidget(self.triangles_lbl)

        # window layout
        window = QWidget(self)
        layout = QHBoxLayout(window)
        layout.addWidget(self.viewer)
        layout.addWidget(self.sidebar)

        self.setCentralWidget(window)

    def buttonHandler(self):
        sender = self.sender()
        if sender.text() == '&Open':
            filename, _ = QFileDialog.getOpenFileName(
                self, 'Open file', '', 'Geometry files (*.tri)')
            if filename is None or filename == '':
                return
            self.loadFile(filename)

        elif sender == self.camreset_btn:
            self.viewer.resetCamera()

        elif sender == self.wireframe_cb:
            self.viewer.setWireframe(self.sender().isChecked())

        elif sender == self.lighting_cb:
            self.viewer.setLighting(self.sender().isChecked())

        self.viewer.update()
