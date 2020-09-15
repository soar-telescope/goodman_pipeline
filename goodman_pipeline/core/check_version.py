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

API_URL = 'https://api.github.com/repos/soar-telescope/goodman/releases/latest'


def get_last(github_api_token='GITHUB_ACCESS_TOKEN'):
    """
    Returns the version of the last release on GitHub.

    Parameters
    ----------
        github_api_token (str, optional) : Name of the environment variable
        holding the github access token for the API

    Returns
    -------
        version (LooseVersion) : the last version of the pipeline.
    """
    try:
        access_token = os.environ[github_api_token]
        headers = {'Authorization': 'token {}'.format(access_token)}
    except KeyError:
        headers = {}

    response = requests.get(API_URL, headers=headers)

    if response.status_code != 200:  # pragma: no cover
        raise ConnectionRefusedError('Number of tests reached maximum for now.')

    tag_name = response.json()['tag_name'].replace('v', '')
    _version = LooseVersion(tag_name)

    return _version.vstring


def am_i_updated(version):

    version = LooseVersion(version.replace('v', ''))
    last_version = get_last()

    return last_version <= version
