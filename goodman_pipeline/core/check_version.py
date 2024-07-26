"""
v1.0.0
    - First public release.
"""

import logging
import requests
import os

from packaging.version import Version

logger = logging.getLogger(__name__)

API_URL = 'https://api.github.com/repos/soar-telescope/goodman/releases/latest'


def get_last(github_api_token: str = 'GITHUB_ACCESS_TOKEN') -> Version:
    """
    Returns the version of the last release on GitHub.

    Parameters
    ----------
        github_api_token (str, optional) : Name of the environment variable
          holding the github access token for the API

    Returns
    -------
        version (object) : A :class:`pkg_resources.extern.packaging.version.Version` the last version of the pipeline.
    """
    try:
        access_token = os.environ[github_api_token]
        headers = {'Authorization': 'token {}'.format(access_token)}
    except KeyError:
        headers = {}

    response = requests.get(API_URL, headers=headers, timeout=3)

    if response.status_code != 200:  # pragma: no cover
        raise ConnectionRefusedError('Number of tests reached maximum for now.')

    tag_name = response.json()['tag_name'].replace('v', '')
    _version = Version(tag_name)

    return _version


def am_i_updated(version: str) -> bool:

    version = Version(version)
    last_version = get_last()

    return last_version <= version
