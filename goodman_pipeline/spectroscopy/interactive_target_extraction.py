import copy
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PySimpleGUI as sg


from astropy.visualization import ZScaleInterval
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.widgets import Button

from ..core import extract_fractional_pixel

#
# class InteractiveExtractionGUI(object):
#     def __init__(self):
#         self.window = None
#         self.fig_agg = None
#
#         self.extraction_width = 2
#         self.background_spacing = 3
#
#         self.fig = None
#         self.ax1 = None
#         self.ax2 = None
#         self.ax3 = None
#         self.ax4 = None
#         self.image_bb = None
#         self.spatial_profile_bb = None
#
#         self.ccd = None
#
#         self.trace = None
#         self.profile = None
#
#         self.scale = ZScaleInterval()
#         plt.style.use('dark_background')
#
#         self.layout = [[sg.Button('Settings'), sg.Button('Save & Exit'), sg.Button('Help')],
#                        [sg.Canvas(key='fig_canvas')]]
#         self.window = sg.Window(title="Interactive Target Extraction",
#                                 layout=self.layout,
#                                 finalize=True,
#                                 resizable=True, location=(100, 100))
#
#     def __call__(self, ccd, trace, *args, **kwargs):
#         self.ccd = ccd
#         self.fig = None
#         self.fig, ((self.ax1, self.ax2), (self.ax3, self.ax4)) = plt.subplots(
#             nrows=2,
#             ncols=2,
#             sharey='row',
#             sharex='col',
#             gridspec_kw={'width_ratios': [3, 1]},
#             figsize=(16, 9))
#
#         spatial_profile = np.median(ccd.data, axis=1)
#         self.ax1.set_title(f"2D Reduced Image")
#         self.ax1.imshow(self.ccd.data, clim=self.scale.get_limits(self.ccd.data), cmap='gray', aspect='auto')
#         self.ax2.set_title(f"Spatial Profile")
#         self.ax2.plot(spatial_profile, range(len(spatial_profile)))
#
#         for self.trace, self.profile, trace_info in trace:
#             self.ax1.plot(self.trace(range(ccd.data.shape[1])), color='r')
#             try:
#                 self.ax2.axhline(self.profile.mean.value, color='r')
#             except AttributeError as error:
#                 self.ax2.axhline(self.profile.x_0.value, color='r')
#
#             # self._perform_extraction(offsets=[0])
#         # this to a function??
#         self.fig_agg = draw
#
#         while True:
#             event, values = self.window.read(timeout=200)
#             if event == sg.WIN_CLOSED or event == 'Exit':
#                 break
#             if event == 'Help':
#                 print("Help")
#
#     def _draw_figure(self, canvas, figure):
#         figure_canvas_agg = FigureCanvasQTAgg(figure, canvas)
#         figure_canvas_agg.draw()
#         print(dir(figure_canvas_agg))






class InteractiveExtraction(object):

    def __init__(self):
        # settings
        self.extraction_width = 2
        self.background_spacing = 3
        self.use_trace_in_preview = False
        self.use_single_background = False
        self.extraction_width_preview = 10 # pixels
        self.profile_section = 'full'
        self.profile_section_size = 100

        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ax4 = None
        self.image_bb = None
        self.spatial_profile_bb = None
        self.preview_mode = False
        self.preview_line = None
        self.preview_line_image = None
        self.preview_line_profile = None
        self.color_index = 0

        # Buttons
        self.btn_preview = None

        self.ccd = None
        self.spatial_profile = None

        self.all_traces = None
        self.trace = None
        self.profile = None
        self.trace_info = None

        self.scale = ZScaleInterval()
        plt.style.use('dark_background')

        self.targets = []
        self.extracted = []

    class Targets:
        def __init__(self, trace=None, profile=None, bkg_info=None, target_type='initial'):
            self.color_index = 1
            self.trace = trace
            self.profile = profile
            self.bkg_info = bkg_info
            self.target_type = target_type
            self.extracted = None
            self.background = None

            # lines
            self.trace_line = None
            self.profile_center_line = None
            self.extracted_line = None
            self.background_line = None
            self.background_spans = []

        @property
        def profile_center(self):
            try:
                return self.profile.mean.value
            except AttributeError as error:
                return self.profile.x_0.value

        def clear_lines(self):
            for l in [self.trace_line,
                      self.profile_center_line,
                      self.extracted_line,
                      self.background_line]:
                if l is not None:
                    l.remove()
            for b in self.background_spans:
                b.remove()
            self.trace_line = None
            self.profile_center_line = None
            self.extracted_line = None
            self.background_line = None
            self.background_spans = []


    def __call__(self, ccd, lamps, traces, *args, **kwargs):
        self.ccd = ccd
        self.fig = None

        self.fig, ((self.ax1, self.ax2), (self.ax3, self.ax4)) = plt.subplots(
            nrows=2,
            ncols=2,
            sharey='row',
            sharex='col',
            gridspec_kw={'width_ratios': [3, 1]},
            figsize=(16, 9))
        self.ax4.remove()
        self.fig.canvas.set_window_title(f"Interactive Target Extraction - {self.ccd.header['OBJECT']} {self.ccd.header['WAVMODE']}")

        manager = plt.get_current_fig_manager()
        if plt.get_backend() == u'GTK3Agg':
            manager.window.maximize()
        elif plt.get_backend() == u'Qt5Agg':
            manager.window.showMaximized()
        self.spatial_profile = np.median(self.ccd.data, axis=1)

        for self.trace, self.profile, self.bkg_info in traces:
            target = self.Targets(trace=self.trace, profile=self.profile, bkg_info=self.bkg_info, target_type='initial')
            target.color_index = len(self.targets) + 1
            self.targets.append(target)

        self._perform_extraction()

        # Help Button
        #  [left, bottom, width, height]
        ax_help = plt.axes([0.76, 0.43, 0.07, 0.05])
        btn_help = Button(ax_help, "Help", color='dimgray', hovercolor='k')
        btn_help.on_clicked(self.__display_help)

        ax_settings = plt.axes([0.84, 0.43, 0.07, 0.05])
        btn_settings = Button(ax_settings, "Settings", color='dimgray', hovercolor='k')
        btn_settings.on_clicked(self.__settings)

        ax_save = plt.axes([0.92, 0.43, 0.07, 0.05])
        btn_save = Button(ax_save, "Save", color='dimgray', hovercolor='k')
        btn_save.on_clicked(self.__display_help)

        # second row

        # ax_clear = plt.axes([0.76, 0.36, 0.07, 0.05])
        # btn_clear = Button(ax_clear, "Clear", color='dimgray', hovercolor='k')
        # btn_clear.on_clicked(self.__display_help)
        #
        ax_preview = plt.axes([0.84, 0.36, 0.07, 0.05])
        self.btn_preview = Button(ax_preview, "Preview Mode", color='dimgray', hovercolor='k')
        self.btn_preview.on_clicked(self.__preview_mode)

        ax_clear = plt.axes([0.92, 0.36, 0.07, 0.05])
        btn_clear = Button(ax_clear, "Clear", color='dimgray', hovercolor='r')
        btn_clear.on_clicked(self.__clear)

        # Third row

        ax_identify = plt.axes([0.76, 0.29, 0.07, 0.05])
        btn_identify = Button(ax_identify, "Identify", color='dimgray', hovercolor='k')
        btn_identify.on_clicked(self.__display_help)

        ax_trace = plt.axes([0.84, 0.29, 0.07, 0.05])
        btn_trace = Button(ax_trace, "Trace", color='dimgray', hovercolor='k')
        btn_trace.on_clicked(self.__preview_mode)

        ax_extract = plt.axes([0.92, 0.29, 0.07, 0.05])
        btn_extract = Button(ax_extract, "Extract", color='dimgray', hovercolor='k')
        btn_extract.on_clicked(self.__extract)

        plt.subplots_adjust(left=0.05,
                            right=0.99,
                            top=0.96,
                            bottom=0.04,
                            hspace=0.17,
                            wspace=0.11)



        self.image_bb = self.ax1.get_position()
        self.spatial_profile_bb = self.ax2.get_position()

        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_press_event', self.key_pressed)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_over)
        self.fig.tight_layout()
        plt.show()
        return True

    def __draw_plots(self):


        self.ax1.set_title(f"2D Image")
        self.ax2.set_title(f"Spatial Profile")

        self.ax1.imshow(self.ccd.data, clim=self.scale.get_limits(self.ccd.data), cmap='gray', aspect='auto')
        self.ax2.plot(self.spatial_profile, range(len(self.spatial_profile)), color=f"C0")

        self.ax3.set_title(f"Extraction Preview")
        for i in range(len(self.targets)):
            self.targets[i].clear_lines()
            self.targets[i].trace_line, = self.ax1.plot(self.targets[i].trace(range(self.ccd.data.shape[1])), color=f"C{self.targets[i].color_index}")
            self.targets[i].profile_center_line = self.ax2.axhline(self.targets[i].profile_center, color=f"C{self.targets[i].color_index}")

            self.targets[i].extracted_line, = self.ax3.plot(self.targets[i].extracted.data, label="Extracted")
            self.targets[i].background_line, = self.ax3.plot(self.targets[i].background, label="Background")
            for bkg in self.targets[i].bkg_info:
                if bkg is not None:
                    bkg_values = bkg.split(' ')[0]
                    lim_1, lim_2 = [float(x) for x in bkg_values.split(':')]
                    bk1 = self.ax1.axhspan(lim_1, lim_2, label=f"Background zone: {bkg_values}", facecolor=f"C{self.targets[i].color_index}",
                                     alpha=0.5)
                    bk2 = self.ax2.axhspan(lim_1, lim_2, label=f"Background zone: {bkg_values}", facecolor=f"C{self.targets[i].color_index}",
                                     alpha=0.5)
                    self.targets[i].background_spans.extend([bk1, bk2])
        self.ax1.legend(loc=2)
        self.ax3.legend(loc=2)
        self.fig.canvas.draw()


    @staticmethod
    def __display_help(event):
        file_path = os.path.join(os.path.dirname(__file__), 'interactive_extraction_help.txt')
        with open(file_path, 'r') as help:
            help_text = help.read()
            sg.theme('DarkTeal6')
            sg.set_options(font=("sans-serif", 13))
            sg.popup_scrolled(help_text, title=f"Interactive Extraction Help", size=(80, None), location=(100, 100))

    def __extract(self, event):
        self._perform_extraction()

    def __preview_mode(self, event):
        self.preview_mode = not self.preview_mode
        if self.preview_mode:
            self.btn_preview.color = 'red'
            for t in self.targets:
                t.extracted_line.remove()
                t.background_line.remove()
        else:
            self.btn_preview.color = 'dimgray'
            for t in self.targets:
                self.ax3.add_line(t.extracted_line)
                self.ax3.add_line(t.background_line)
                self.fig.canvas.draw()
                self.ax3.relim()
                self.ax3.autoscale()
        if not self.preview_mode and self.preview_line is not None:
            self.preview_line.remove()
            self.preview_line_image.remove()
            self.preview_line_profile.remove()
            # self.ax3.relim()

    def __settings(self, event):
        sg.theme('DarkTeal6')
        sg.set_options(font=("sans-serif", 13))
        layout = [[[sg.Text("Extraction Settings:",  font=("Helvetica", 15), justification='left')],
                   [sg.Text("Extraction Width: ", size=(20, 1)), sg.Spin(values=list(range(1, 10, 1)), key='-EXT-WIDTH-', initial_value=self.extraction_width)],
                   [sg.Text("Backround Spacing: ", size=(20, 1)), sg.Spin(values=list(range(1, 10, 1)),key='-BKG-SPACE-', initial_value=self.background_spacing)],
                   [sg.Checkbox("Use Single Background", default=self.use_single_background)],
                   [sg.Text('_' * 30)],
                   [sg.Text("Preview Mode:", font=("Helvetica", 15), justification='left')],
                   [sg.Checkbox("Use Trace", default=self.use_trace_in_preview)],
                   [sg.Text("Extraction Width Pixels"), sg.Spin(values=list(range(1, 50, 1)), key='-EXT-WIDTH-PREVIEW-', initial_value=self.extraction_width_preview)],
                   [sg.Text('_' * 30)],
                   [sg.Text("Profile Calculation:", font=("Helvetica", 15), justification='left')],
                   [sg.Text("Section to use: ", size=(20, 1)), sg.Spin(values=['full', 'center', 'start'], size=(7,1), key='-PROFILE-SECTION-', initial_value=self.profile_section)],
                   [sg.Text("Size (if not full): ", size=(20, 1)),
                    sg.Spin(values=list(range(1,100,1)), size=(7, 1), key='-PROFILE-SECTION-SIZE-',
                            initial_value=self.profile_section_size)],
                   [sg.Button("Ok")]]]

        window = sg.Window("Settings", layout, location=(100, 100))
        while True:
            event, values = window.read()
            self.extraction_width = values['-EXT-WIDTH-']
            self.background_spacing = values['-BKG-SPACE-']
            self.extraction_width_preview = values['-EXT-WIDTH-PREVIEW-']
            if event == sg.WIN_CLOSED or event == 'Ok':
                break

        window.close()
        self._perform_extraction()

    def __clear(self, event):
        for t in self.targets:
            if t.target_type != 'initial':
                t.clear_lines()
                self.targets.remove(t)

        self.__draw_plots()


    def on_click(self, event):
        print(event)
        if event.xdata is not None and event.ydata is not None:
            figure_x, figure_y = \
                self.fig.transFigure.inverted().transform((event.x, event.y))
            if self.image_bb.contains(figure_x, figure_y) and event.button == 2:
                self.shift_trace(event)
                self.color_index += 1

    def on_mouse_over(self, event):
        if self.preview_mode:
            ax1_x, ax1_y = self.fig.transFigure.inverted().transform((event.x, event.y))

            if self.image_bb.contains(ax1_x, ax1_y):
                for line in [self.preview_line, self.preview_line_image, self.preview_line_profile]:
                    if line is not None:
                        try:
                            line.remove()
                            self.ax3.relim()
                        except ValueError as error:
                            pass
                if event.ydata is not None:
                    y0 = np.max([0, int(event.ydata - 0.5 * self.extraction_width_preview)])
                    y1 = np.min([self.ccd.data.shape[0],  int(event.ydata + 0.5 * self.extraction_width_preview)])
                    sample = np.sum(self.ccd.data[y0:y1, :], axis=0)
                    self.preview_line, = self.ax3.plot(sample, color='r')
                    self.preview_line_image = self.ax1.axhline(event.ydata, color='r')
                    self.preview_line_profile = self.ax2.axhline(event.ydata, color='r')
                    self.fig.canvas.draw()



    def key_pressed(self, event):
        print(event.key)
        if event.key == 'ctrl+q':
            sys.exit()
        elif event.key == 'ctrl+p':
            print("Preview Mode")
        elif event.key == 'h':
            print("Help")
            self.__display_help()
        elif event.key == 's':
            self.__settings()

    def shift_trace(self, event):
        offset = event.ydata - self.trace(event.xdata)

        shifted_trace = copy.deepcopy(self.trace)
        shifted_profile = copy.deepcopy(self.profile)
        try:
            shifted_profile.mean.value += offset
            profile_center = shifted_profile.mean.value
        except AttributeError as error:
            shifted_profile.x_0.value += offset
            profile_center = shifted_profile.x_0.value
        try:
            if offset != 0:
                shifted_trace.c0.value += offset
        except AttributeError as error:
            print(error)

        target = self.Targets(trace=shifted_trace, profile=shifted_profile, target_type='additional')
        target.color_index = len(self.targets) + 1
        self.targets.append(target)

        self._perform_extraction()


    def _perform_extraction(self):

        for t in self.targets:
            extracted, background, bkg_info = extract_fractional_pixel(
                ccd=self.ccd,
                target_trace=t.trace,
                target_fwhm=t.profile.fwhm,
                extraction_width=self.extraction_width,
                background_spacing=self.background_spacing)

            t.extracted = extracted
            t.background = background
            t.bkg_info = bkg_info

        self.__draw_plots()







if __name__ == '__main__':
    ie = InteractiveExtraction()
    ie()