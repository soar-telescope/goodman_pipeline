"""Test for Process

Returns:
    True value

"""
from astropy.io import fits
from redspec import ScienceObject
from process import Process
from sphinx import cmdline

single = True

if single:
    sourcepath = '/data/simon/data/soar/work/goodman/test/'
    obj = ScienceObject('R339704',
                        'fc_0045.SO2016A-019_0320.fits',
                        '2016-03-20T23:54:15.96',
                        65.69922083333333,
                        -76.57505888888889)
    # obj.add_lamp('fc_0046.SO2016A-019_0320.fits', 'HgAr', 65.69922083333333, -76.57505888888889)
    obj.add_lamp('fc_0047.SO2016A-019_0320.fits', 'CuHeAr', 65.69922083333333, -76.57505888888889)
else:
    sourcepath = '/user/simon/data/SOAR/work/multi-objects/'
    obj = ScienceObject('CVSO-34',
                        'bc_0099.cvso34_400M2_GG455.fits',
                        '2015-01-04T05:33:57.22',
                        81.41224583333333,
                        1.7139411111111111)
    obj.add_lamp('bc_0110.near_400M2_GG455.fits', 'NeAr', 81.45332916666668, 1.8253911111111112)
    obj.add_lamp('bc_0111.near_400M2_GG455.fits', 'NeAr', 81.45332916666668, 1.8253911111111112)

obj.print_all()


p = Process(sourcepath,obj)
#targets = p.identify_spectra()
#for target in targets:
#    print(target)