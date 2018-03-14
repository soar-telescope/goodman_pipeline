import glob
import os
from watchdog.observers import Observer
from watchdog.events import (LoggingEventHandler,
                             FileSystemEventHandler,
                             FileModifiedEvent, FileCreatedEvent)
import time
import logging
from ccdproc import CCDData

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')


class EventHandler(FileSystemEventHandler):

    def __init__(self):
        self.log = logging.getLogger(__name__)

    def dispatch(self, event):
        if isinstance(event, FileCreatedEvent):
            time.sleep(1)
            self.print_object(event)
        else:
            if event.is_directory:
                self.log.debug(event.__class__.__name__)
            else:
                self.log.info("File {:s}: {:s}".format(event.event_type,
                                                       os.path.basename(event.src_path)))
            # print(dir(event))
            # self.log.info()
            # self.log.info(event.__class__.__name__)

    def print_object(self, event):
        # print(event)
        # print(event.src_path)
        if '.fits' in event.src_path:
            # print(time.time())
            try:
                self.log.info("Opening File: {:s}".format(os.path.basename(event.src_path)))
                ccd = CCDData.read(event.src_path, unit='adu')
                self.log.info("{:s} {:s} {:d}".format(ccd.header['OBJECT'],
                                                      ccd.header['OBSTYPE'],
                                                      ccd.data.max()))

            except FileNotFoundError:
                self.log.error("Problem finding File: {:s}".format(event.src_path))
        else:
            self.log.warning('Not a fits file')
        # ccd = CCDData()


class MainApp(object):

    def __init__(self, folder):
        self.path = folder
        # self.event_handler = LoggingEventHandler()
        self.handler = EventHandler()
        self.observer = Observer()
        self.observer.schedule(self.handler, self.path)

    def __call__(self, *args, **kwargs):
        # print(glob.glob(os.path.join(self.path, '*')))
        self.observer.start()
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            self.observer.stop()
        self.observer.join()


if __name__ == '__main__':
    path = '/data/simon/Downloads/testing_watchdog'
    live = MainApp(path)
    live()
