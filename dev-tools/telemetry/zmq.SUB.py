#!/usr/bin/python

#ZMQ client

import zmq
import sys
import time

server = "tcp://localhost:5556"

context = zmq.Context()
socket = context.socket(zmq.SUB)
print("Initializing subscription to %s"%server)
socket.connect(server)
socket.setsockopt(zmq.SUBSCRIBE,'INSERT')



while True:
    try:
        message = socket.recv()
        print message
        time.sleep(0.5)
    except KeyboardInterrupt:
        print("")
        print("The Ending subscription to %s"%server)
        sys.exit(0)
