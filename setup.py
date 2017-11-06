from distutils.core import setup

setup(
    name='goodman',
    version='1.0b2',
    packages=['goodman', 'goodman.core', 'goodman.images', 'goodman.spectroscopy', 'goodman.telemetry'],
    package_dir={'goodman': 'goodman'},
    package_data={'goodman': ['data/params/dcr.par', 'data/ref_comp/*fits']},
    scripts=['bin/redccd', 'bin/redspec'],
    url='https://github.com/soar-telescope/goodman',
    license='BSD 3-Clause',
    author='Simon Torres R.',
    author_email='storres@ctio.noao.edu',
    description='Pipelines for CCD and Spectroscopic Reduction of Goodman Data'
)
