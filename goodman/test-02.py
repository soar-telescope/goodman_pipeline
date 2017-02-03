from astropy.io import fits
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import glob

# image = '20161125_multiobject_n2_s2-fwhm.fits'
image_list = glob.glob('f*fits')
for image in image_list:
    # image = '20161125_multiobject_n2_s1-fwhm.fits'
    data = fits.getdata(image)

    sample = np.median(data[:, 1900:2100], axis=1)

    # foo = np.where(np.abs(data > data.min() + 0.75 * data.max()), data, 0)
    foo = np.where(np.abs(sample > sample.min() + 0.75 * sample.max()), sample, 0)

    maxima = signal.argrelmax(foo, axis=0, order=10)[0]
    voigt_init = None
    fit_voigt = fitting.LevMarLSQFitter()
    x_axis = range(len(sample))
    print dir(maxima)

    if maxima != []:
        print maxima, len(list(maxima))
        for val in maxima:
            if voigt_init is None:
                voigt_init = models.Voigt1D(amplitude_L=sample[val], x_0=val)
            else:
                voigt_init += models.Voigt1D(amplitude_L=sample[val], x_0=val)
            plt.axvline(val)
        voigt = fit_voigt(voigt_init, x_axis, sample)
        for sub_fit in voigt:
            plt.plot(x_axis, sub_fit(x_axis), linestyle='--', color='g')
        # print voigt[0]
        # plt.plot(x_axis, voigt(x_axis), color='r')


    else:
        print('No targets detected')



    print sample.shape

    # import astropy.io.fits as pyfits
    # import numpy as np
    # spec = pyfits.getdata(filename)
    # foo = np.where(np.abs(spec > spec.min() + 0.75 * spec.max()), spec, 0)
    # maxima = signal.argrelmax(foo, axis=0, order=10)[0]

    plt.plot(sample)
    plt.title(image)
    plt.show()

