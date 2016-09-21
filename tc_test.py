import tc_read_IRAF_spec
from specutils.io import read_fits
filename = '/data/simon/data/soar/work/goodman/test/extraction-tests/cuhear600nonlinearli.fits'
#tc_read_IRAF_spec.read_IRAF_spec(filename)

data = read_fits.read_fits_spectrum1d(filename)