"""
v1.0.0
    - First public release.
"""

import logging
import requests
import os
import re

from distutils.version import LooseVersion

logger = logging.getLogger(__name__)

try:
    ACCESS_TOKEN = os.environ['GITHUB_ACCESS_TOKEN']
    LATEST_URL = \
        'https://api.github.com/repos/soar-telescope/goodman/releases/latest' \
        '?access_token={:s}'.format(ACCESS_TOKEN)
except KeyError:
    LATEST_URL = \
        'https://api.github.com/repos/soar-telescope/goodman/releases/latest'


def get_last(url=LATEST_URL):
    """
    Returns the version of the last release on GitHub.

    Parameters
    ----------
        url (str, optinal) : the URL that is used to retrieve the information

    Returns
    -------
        version (LooseVersion) : the last version of the pipeline.
    """

    response = requests.get(url)

    if response.status_code != 200:
        raise ConnectionRefusedError('Number of tests reached maximum for now.')

    tag_name = response.json()['tag_name'].replace('v', '')
    _version = LooseVersion(tag_name)

    return _version.vstring


def am_i_updated(version):

    version = LooseVersion(version.replace('v', ''))
    last_version = get_last()

    return last_version <= version
