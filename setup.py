from distutils.core import setup

setup(
    name='goodman',
    version='1.0b2',
    packages=['ccd', 'spectroscopy'],
    package_dir={'ccd': 'ccd',
                 'spectroscopy': 'spectroscopy'},
    package_data={'ccd': ['files/dcr.par'],
                  'spectroscopy': ['ref_comp/*fits']},
    scripts=['bin/redccd', 'bin/redspec'],
    url='https://github.com/soar-telescope/goodman',
    license='BSD 3-Clause',
    author='Simon Torres R.',
    author_email='storres@ctio.noao.edu',
    description='Pipelines for CCD and Spectroscopic Reduction of Goodman Data'
)
