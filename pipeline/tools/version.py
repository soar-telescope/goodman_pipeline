
import requests
from bs4 import BeautifulSoup

api = 1
feature = 0
bug = 0

month = 0
year = 0

__str__ = "v{:d}.{:d}.{:d} - {:d}-{:d}".format(api, feature, bug, month, year)


def check_last(url="https://github.com/soar-telescope/goodman/releases/latest"):

    html = requests.get(url)

    return api, feature, bug