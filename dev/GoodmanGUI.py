#!/usr/bin/env python2
import sys
from PyQt4 import QtGui, QtCore


class MainGuiApp(object):

    def __init__(self):
        self.app = QtGui.QApplication(sys.argv)
        self.tabs = QtGui.QTabWidget()
        self.tabs.resize(600, 200)
        self.ccd_tab = QtGui.QWidget()
        self.spec_tab = QtGui.QWidget()

        # data path initialization
        self.raw_path = None
        self.red_path = None
        self.source_path = None
        self.destiny_path = None

        # CCD REDUCTION TAB
        self.ccd_grid = QtGui.QGridLayout()
        self.ccd_grid.setSpacing(1)
        # RAW PATH
        # text
        self.raw_path_label = QtGui.QLabel()
        self.raw_path_label.setText("Raw Data Path ")
        # textBox
        self.raw_path_textbox = QtGui.QLineEdit()
        self.raw_path_textbox.resize(200, 10)
        # pushButton
        self.raw_browse_button = QtGui.QPushButton("Browse")
        # add widgets to browse horizontal layout
        self.ccd_grid.addWidget(self.raw_path_label, 1, 0)
        self.ccd_grid.addWidget(self.raw_path_textbox, 1, 1)
        self.ccd_grid.addWidget(self.raw_browse_button, 1, 2)

        # RED PATH
        # text
        self.red_path_label = QtGui.QLabel()
        self.red_path_label.setText("Reduced Data Path ")
        # textBox
        self.red_path_textbox = QtGui.QLineEdit()
        # self.red_path_textbox.resize(200, 10)
        # pushButton
        self.red_browse_button = QtGui.QPushButton("Browse")
        # add widgets to browse horizontal layout
        self.ccd_grid.addWidget(self.red_path_label, 2, 0)
        self.ccd_grid.addWidget(self.red_path_textbox, 2, 1)
        self.ccd_grid.addWidget(self.red_browse_button, 2, 2)

        # Add checkboxes and go button

        # Final buttons
        self.ccd_reset_button = QtGui.QPushButton("Reset")
        self.ccd_go_button = QtGui.QPushButton("Go")

        # self.ccd_grid.addWidget(self.ccd_reset_button, 4, 2)
        # self.ccd_grid.addWidget(self.ccd_go_button, 5, 2)

        small_grid = QtGui.QGridLayout()
        small_grid.setSpacing(10)

        small_grid.addWidget(self.ccd_reset_button, 1, 0)
        small_grid.addWidget(self.ccd_go_button, 1, 1)

        vlayout = QtGui.QVBoxLayout()
        vlayout.addItem(self.ccd_grid)
        vlayout.addItem(small_grid)



        # self.ccd_tab.setLayout(self.ccd_grid)
        self.ccd_tab.setLayout(vlayout)
        # END CCD REDUCTION TAB

        # SPECTROSCOPIC REDUCTION TAB
        self.spec_grid = QtGui.QGridLayout()
        self.spec_grid.setSpacing(1)
        # SOURCE PATH
        # text
        self.source_path_label = QtGui.QLabel()
        self.source_path_label.setText("Source Data Path ")
        # textBox
        self.source_path_textbox = QtGui.QLineEdit()
        self.source_path_textbox.resize(200, 10)
        # pushButton
        self.source_browse_button = QtGui.QPushButton("Browse")
        # add widgets to browse horizontal layout
        self.spec_grid.addWidget(self.source_path_label, 1, 0)
        self.spec_grid.addWidget(self.source_path_textbox, 1, 1)
        self.spec_grid.addWidget(self.source_browse_button, 1, 2)

        # DESTINATION PATH
        # text
        self.destiny_path_label = QtGui.QLabel()
        self.destiny_path_label.setText("Destination Data Path ")
        # textBox
        self.destiny_path_textbox = QtGui.QLineEdit()
        # self.destiny_path_textbox.resize(200, 10)
        # pushButton
        self.destiny_browse_button = QtGui.QPushButton("Browse")
        # self.destiny_browse_button.connect(self.spec_tab, self.on_browse_click)
        # add widgets to browse horizontal layout
        self.spec_grid.addWidget(self.destiny_path_label, 2, 0)
        self.spec_grid.addWidget(self.destiny_path_textbox, 2, 1)
        self.spec_grid.addWidget(self.destiny_browse_button, 2, 2)

        # Add checkboxes and go button

        # Final buttons
        self.ccd_reset_button = QtGui.QPushButton("Reset")
        self.ccd_go_button = QtGui.QPushButton("Go")
        self.spec_grid.addWidget(self.ccd_reset_button, 4, 2)
        self.spec_grid.addWidget(self.ccd_go_button, 5, 2)

        self.spec_tab.setLayout(self.spec_grid)

        # connect buttons to events
        # ccd mode
        self.ccd_tab.connect(self.raw_browse_button, QtCore.SIGNAL('clicked()'), self.on_crb_click)
        self.ccd_tab.connect(self.red_browse_button, QtCore.SIGNAL('clicked()'), self.on_crrb_click)

        # spectroscopic mode
        self.spec_tab.connect(self.source_browse_button, QtCore.SIGNAL('clicked()'), self.on_ssb_click)
        self.spec_tab.connect(self.destiny_browse_button, QtCore.SIGNAL('clicked()'), self.on_sdb_click)


        # Add tabs
        self.tabs.addTab(self.ccd_tab, "CCD Reduction")
        self.tabs.addTab(self.spec_tab, "Spectroscopic Reduction")

        # set title and show
        self.tabs.setWindowTitle("Goodman Data Reduction")
        self.tabs.show()

    def __call__(self, *args, **kwargs):
        sys.exit(self.app.exec_())

    # Event handlers
    def on_crb_click(self):
        """Event handler for raw data directory selection in ccd mode"""
        self.raw_path = QtGui.QFileDialog.getExistingDirectory()
        self.raw_path_textbox.setText(self.raw_path)
        if self.red_path_textbox.text() == '' and self.raw_path != '':
            self.red_path = self.raw_path + '/RED'
            self.red_path_textbox.setText(self.red_path)

    def on_crrb_click(self):
        """Event handler for reduction data directory selection in ccd mode"""
        self.red_path = QtGui.QFileDialog.getExistingDirectory()
        self.red_path_textbox.setText(self.red_path)

    def on_ssb_click(self):
        """Event handler for source directory selection in spectroscopic mode"""
        self.source_path = QtGui.QFileDialog.getExistingDirectory()
        self.source_path_textbox.setText(self.source_path)
        if self.destiny_path_textbox.text() == '' and self.source_path != '':
            self.destiny_path = self.source_path
            self.destiny_path_textbox.setText(self.destiny_path)

    def on_sdb_click(self):
        """Event handler for destination directory selection in spectroscopic mode"""
        self.destiny_path = QtGui.QFileDialog.getExistingDirectory()
        self.destiny_path_textbox.setText(self.destiny_path)




if __name__ == '__main__':
    GUI_APP = MainGuiApp()
    GUI_APP()
