from linelist import ReferenceData
from astropy.io import fits

filename = '/user/simon/data/soar/work/goodman/test/extraction-tests/eefc_0047.SO2016A-019_0320.fits'
#tc_read_IRAF_spec.read_IRAF_spec(filename)

header = fits.getheader(filename)

refdata = ReferenceData(args=[])

name = header['OBJECT']
print name
refdata.get_ref_spectrum_from_linelist(3087, 6224, name)
# with argon fails, too many lines????