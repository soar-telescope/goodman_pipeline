import pygtk
import gtk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
pygtk.require('2.0')


class Base:

    def __init__(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_position(gtk.WIN_POS_CENTER)
        self.window.set_size_request(1200, 800)
        self.window.set_title("Interactive Wavelength Calibration")
        self.window.set_tooltip_text("Goodman's Interactive Wavelength Calibration")
        # self.window.fullscreen()

        # self.menu = gtk.Menu()

        # self.

        self.button1 = gtk.Button('EXIT')
        self.button1.connect("clicked", self.destroy)
        self.button1.set_tooltip_text("This button will close this window")

        self.button2 = gtk.Button('Hide')
        self.button2.connect("clicked", self.myhide)

        self.button3 = gtk.Button('Show')
        self.button3.connect("clicked", self.myshow)

        self.button4 = gtk.Button("Rename Label")
        self.button4.connect("clicked", self.relabel)

        self.button5 = gtk.Button("Clear Text")
        self.button5.connect("clicked", self.clear_text)

        self.label1 = gtk.Label("New Label")

        self.textbox = gtk.Entry()
        self.textbox.connect("changed", self.textchange)

        # fixed = gtk.Fixed()
        # fixed.put(self.button1, 20, 30)
        # fixed.put(self.button2, 100, 30)
        # fixed.put(self.button3, 180, 30)

        self.vbox1 = gtk.VBox(False, 0)
        # self.box1.pack_start(self.button1)
        # self.box1.pack_start(self.button2)
        # self.box1.pack_start(self.button3)
        # self.box1.pack_start(self.label1)
        # self.box1.pack_start(self.button4)
        # self.box1.pack_start(self.textbox)
        # self.box1.pack_start(self.button5)

        self.menu = gtk.Menu()
        menu_item = gtk.MenuItem('Load')
        menu_item2 = gtk.MenuItem('Help')

        self.menu.append(menu_item)
        self.menu.append(menu_item2)

        root_menu = gtk.MenuItem('Root Menu')
        root_menu.show()

        root_menu.set_submenu(self.menu)

        self.menu_bar = gtk.MenuBar()
        self.menu_bar.append(root_menu)
        self.vbox1.pack_start(self.menu_bar, False, False, 2)
        # self.window.add(self.menu_bar)
        # self.menu_bar.show()

        file_item = gtk.MenuItem("File")
        # file_item.show()

        self.fig1 = plt.figure(1)
        x_axis = np.range(1,25)
        self.fig1.plot(x_axis, np.sin(x_axis))



        self.window.add(self.vbox1)
        # self.window.add(fixed)
        self.window.show_all()
        self.window.connect("destroy", self.destroy)

    def destroy(self, widget, data=None):
        gtk.main_quit()

    def myhide(self, widget):
        self.button1.hide()

    def myshow(self, widget):
        self.button1.show()

    def relabel(self, widget):
        self.label1.set_text("New Text")

    def textchange(self, widget):
        self.window.set_title(self.textbox.get_text())
        self.label1.set_text(self.textbox.get_text())

    def clear_text(self, widget):
        self.textbox.set_text("")

    def main(self):
        gtk.main()

if __name__ == '__main__':
    base = Base()
    base.main()