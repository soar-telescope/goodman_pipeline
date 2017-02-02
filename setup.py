from distutils.core import setup

setup(
    name='goodman',
    version='0.1',
    packages=['goodman_ccd', 'goodman_spec'],
    scripts=['bin/redccd', 'bin/redspec'],
    url='https://github.com/simontorres/goodman',
    license='BSD 3-Clause',
    author='Simon Torres R.',
    author_email='storres@ctio.noao.edu',
    description='Pipelines for CCD and Spectroscopic Reduction of Goodman Data'
)
