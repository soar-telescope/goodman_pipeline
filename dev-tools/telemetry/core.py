import zmq
import time
import sys
import os

# import MySQLdb
import json

import random

import pkg_resources


def load_config(name, config_path=None):
    """Load configuration dictionary from json file

    Args:
        name (str): name of configuration stored in the json file
        config_file (str): json configuration file name

    Returns:
        A dictionary of the configuration defined by **name**. If it doesn't
          exist it will not return anything.

    """
    if config_path is None:
        config_path = pkg_resources.resource_filename('pipeline',
                                                      'data/params/config.json')

    if os.path.isfile(config_path):
        with open(config_path) as json_config:
            config = json.load(json_config)
            try:
                return config[name]
            except KeyError:
                print('Configuration for {:s} does not exist.'.format(name))
    else:
        print("Configuration file {:s} does not exist.".format(config_path))
        raise FileNotFoundError


class ZmqPublisher(object):
    """ZMQ Publisher implementation"""

    def __init__(self, ip=None, port=None):
        """Initializes the publisher instance

        A publisher instance requires an ip range and a port where to
        broadcast messages

        Args:
            ip (str): IP address range, can be a single ip also. For the full
              range use '*'
            port (str): Port number usually something like '5555'
        """

        if ip is None and port is None:
            self.config = load_config('publisher')
            ip = self.config['ip']
            port = self.config['port']

        self.server = "tcp://{:s}:{:s}".format(ip, port)
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.PUB)
        self.socket.bind(self.server)

    def broadcast(self, message):
        """Broadcast the message

        Args:
            message (str): The message to be sent.
        """
        self.socket.send(message)

class ZmqSubscriber(object):
    """ZMQ Subscriber implementation"""

    def __init__(self, host=None, port=None):
        """Initializes the subscriber instance

        A subscriber needs an IP address and a port. I also needs a filter
        defined by a string. For this implementation it is fixed to use
        messages containing the string INSERT because it is expected to capture
        MySQL queries.

        Args:
            host (str): Publisher server IP.
            port (str): Listening port.
        """
        if host is None and port is None:
            self.config = load_config('subscriber')
            # TODO (simon): do something if self.config is None
            host = self.config['server_ip']
            port = self.config['listening_port']

        self.server = "tcp://{:s}:{:s}".format(host, port)
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.SUB)
        self.socket.connect(self.server)
        # limit subscription to INSERT queries
        # self.socket.setsockopt(zmq.SUBSCRIBE, 'INSERT')

    def listen_and_print(self):
        """Listen to the publisher and prints any incomming message."""

        while True:
            try:
                message = self.socket.recv()
                print(message)
                time.sleep(0.5)
            except KeyboardInterrupt:
                sys.exit("End of subscription to {:s}".format(self.server))

    def listen_and_save(self, db_config=None):
        """Listen to the publisher and forwards the query to the MariaDB server.

        This method creates an infinite loop and once a message arrives it
        forwards the query to the MariaDB server.

        Args:
            db_config (dict): Access information for the database server. The
              fields are: 'host', 'user', 'password' and 'database'.

        """
        if db_config is None:
            db_config = load_config(name='mysql')
        database = DatabaseHandler(host=db_config['host'],
                                   user=db_config['user'],
                                   password=db_config['password'],
                                   database=db_config['database'])
        while True:
            try:
                mysql_query = self.socket.recv()
                database.execute(sql_query=mysql_query)

                time.sleep(0.5)
            except KeyboardInterrupt:
                sys.exit("End of subscription to {:s}".format(self.server))


class DatabaseHandler(object):
    """Database handler implementation"""
    # create table telemetry
    # (ID int NOT NULL AUTO_INCREMENT,
    # INPUT_DATE TIMESTAMP,
    # MODEL_NAME varchar(30) NOT NULL,
    # MODEL_ORDER int,
    # C0 float,
    # C1 float,
    # C2 float,
    # C3 float,
    # C4 float,
    # C5 float,
    # C6 float,
    # C7 float,
    # C8 float,
    # PRIMARY KEY (ID));

    def __init__(self, host=None, user=None, password=None, database=None):
        """Instantiate a MySQLdb connection

        Args:
            host (str): Server IP Address.
            user (str): MariaDB server user name.
            password (str): Password for the user
            database (str): Name of the database
        """
        self.connection = MySQLdb.connect(host=host,
                                          user=user,
                                          passwd=password,
                                          db=database)
        self.cursor = self.connection.cursor()

    def execute(self, sql_query):
        """Executes a query

        Args:
            sql_query (str): Full SQL query.
        """
        try:
            self.cursor.execute(sql_query)
            self.connection.commit()
        except:
            self.connection.rollback()


if __name__ == '__main__':
    json_config = 'config.json'
    with open('config.json') as json_config:
        config = json.load(json_config)

        database = DatabaseHandler(host=config['mysql']['host'],
                                   user=config['mysql']['user'],
                                   password=config['mysql']['password'],
                                   database=config['mysql']['database'])

        # query = ("INSERT INTO telemetry ( MODEL_NAME, MODEL_ORDER, C0, C1, C2) "
        #          "VALUES ('{:s}', '{:d}', '{:.3f}', '{:.3f}', '{:.3f}');").format('Chebyshev1D',
        #                                                                           3,
        #                                                       random.randint(0, 10) + random.random(),
        #                                                       random.randint(0, 10) + random.random(),
        #                                                       random.randint(0, 10) + random.random())

        # print(query)
        # database.execute(query=query)


    publisher = ZmqPublisher('*', '5556')

    try:
        while 1:
            query = ("INSERT INTO telemetry "
                     "(MODEL_NAME, MODEL_ORDER, C0, C1, C2) "
                     "VALUES ('{:s}', '{:d}', '{:.3f}', '{:.3f}', '{:.3f}');"
                     ).format('Chebyshev1D',
                              3,
                              random.randint(0, 10) + random.random(),
                              random.randint(0, 10) + random.random(),
                              random.randint(0, 10) + random.random())

            publisher.broadcast(query)
            time.sleep(3)
    except KeyboardInterrupt:
        sys.exit('KeyboardInterrupt Exit')