import pyinotify
from ccdproc import CCDData
import os

path = '/data/simon/Downloads/testing_watchdog'

file_events = pyinotify.IN_OPEN | pyinotify.IN_CLOSE_WRITE | pyinotify.IN_CLOSE_NOWRITE

wm = pyinotify.WatchManager()


class EventHandler(pyinotify.ProcessEvent):
    def process_IN_DELETE(self, event):
        pass

    def process_IN_CREATE(self, event):
        print("File {:s} was created".format(event.pathname))

    def process_IN_CLOSE_WRITE(self, event):
        # print("This is what I want: {:s}.".format(event.pathname))
        self._print_object(full_path=event.pathname)

    def _print_object(self, full_path):
        ccd = CCDData.read(full_path, unit='adu')
        print(ccd.header['OBJECT'], ccd.header['OBSTYPE'], os.path.basename(full_path))

event_handler = EventHandler()

notifier = pyinotify.Notifier(wm, event_handler)

wm.add_watch(path, file_events)

notifier.loop()