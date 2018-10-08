# -*- coding: utf-8 -*-
"""
SOAR Telescope - Goodman Pipeline.

Goodman High Throughput Spectrograph Data Reduction Pipeline.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

import os
import sys

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

from sphinx.setup_command import BuildDoc
# To use a consistent encoding
from codecs import open


here = os.path.abspath(os.path.dirname(__file__))


def create_version_py(packagename, version, source_dir='.'):
    package_dir = os.path.join(source_dir, packagename)
    version_py = os.path.join(package_dir, 'version.py')

    version_str = "# This is an automatic generated file please do not edit\n" \
                  "__version__ = '{:s}'".format(version)

    with open(version_py, 'w') as f:
        f.write(version_str)

# Get configuration information from setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser
conf = ConfigParser()

# conf.read([os.path.join(os.path.dirname(__file__), '..', 'setup.cfg')])
conf.read([os.path.join(os.path.dirname(__file__), 'setup.cfg')])
metadata = dict(conf.items('metadata'))

PACKAGENAME = metadata['package_name']

VERSION = metadata['version']

LICENSE = metadata['license']

DESCRIPTION = metadata['description']

LONG_DESCRIPTION = metadata['long_description']

AUTHOR = metadata['author']

AUTHOR_EMAIL = metadata['author_email']

# freezes version information in version.py
create_version_py(PACKAGENAME, VERSION)


cmdclassd = {'build_sphinx': BuildDoc,
             'build_docs': BuildDoc}

setup(
    name=metadata['package_name'],

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=VERSION,

    description=DESCRIPTION,

    long_description=LONG_DESCRIPTION,

    # The project's main homepage.
    url='https://github.com/soar-telescope/goodman_pipeline',

    # Author details
    author=u'Simon Torres R., '
           u'Bruno Quint, '
           u'Cesar Briceño, '
           u'David Sanmartin, '
           ,
    cmdclass=cmdclassd,

    author_email='storres@ctio.noao.edu, bquint@ctio.noao.edu, '
                 'cbriceno@ctio.noao.edu',

    # Choose your license
    license=LICENSE,

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',

        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',

        'Natural Language :: English',

        'Operating System :: POSIX :: Linux',
        'Operating System :: POSIX :: Other',
        'Operating System :: MacOS :: MacOS X',

        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',

    ],

    # What does your project relate to?
    keywords='soar pipelines astronomy images spectroscopy',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().

    packages=['goodman_pipeline',
              'goodman_pipeline.core',
              'goodman_pipeline.images',
              'goodman_pipeline.spectroscopy',
              'goodman_pipeline.wcs', ],

    package_dir={'goodman_pipeline': 'goodman_pipeline'},

    package_data={'goodman_pipeline': ['data/params/dcr.par',
                               'data/params/*.json',
                               'data/ref_comp/*fits',
                               'data/dcr-source/README.md',
                               'data/dcr-source/dcr/*',
                               'data/test_data/master_flat/*',
                               'data/test_data/wcs_data/*']},

    scripts=['goodman_pipeline/scripts/redccd',
             'goodman_pipeline/scripts/redspec', ],

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    #install_requires=[''],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # extras_require={
    #    'docs': ['astropy_sphinx_theme'],
    #    'test': ['coverage'],
    # },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'sample': ['package_data.dat'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('dcr-source', ['dcr-source'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
)


