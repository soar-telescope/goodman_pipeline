import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
from matplotlib import gridspec
from numpy.polynomial.chebyshev import chebfit, chebval
from numpy.polynomial.legendre import legfit, legval
from scipy.interpolate import interp1d
import multiprocessing
from multiprocessing import Pool


def fitting(x, y, function='Legendre', order=41):
    # prepare for masking arrays - 'conventional' arrays won't do it
    # data = np.ma.array(data)

    # Assuring that the output yfit has always same number of points of the input
    x = np.arange(0, x.size, 1)

    yfit = []
    if function.lower() == 'legendre':
        yfit = chebval(x, (chebfit(x, y, deg=order)))
    elif function.lower() == 'chebychev':
        yfit = legval(x, (legfit(x, y, deg=order)))
    elif function.lower() == 'spline3':
        nsum = 5

        y_resampled = np.asarray([np.median(y[i:i + nsum]) for i in range(0, len(y) - len(y) % nsum, nsum)])
        x_resampled = np.linspace(0, y.size, y_resampled.size)

        # Masking invalid elements (e.g. with nan)
        y_resampled = ma.masked_invalid(y_resampled)

        # Fitting
        f = interp1d(x_resampled, y_resampled, kind=order, bounds_error=False, fill_value=0.0)

        # Function with the original x array
        yfit = np.asarray(f(x))

        # TODO eliminar esse plots quando a rotina estiver funcionando
        # plotting Spline to see results
        plt.plot(x, y, 'k-', label='Original')
        plt.plot(x_resampled, y_resampled, 'bo', label='Resampled')
        plt.plot(x, yfit, label='Fit')
        plt.legend(loc='best')
        plt.show()

    # Calculating chi2
    npar = order + 1
    sigma2 = 1. / (y.size - 1.) * ((y - y.mean()) ** 2).sum()
    chi2_dof = ((y - yfit) ** 2 / sigma2).sum() / (y.size - npar)

    return yfit, chi2_dof


def sigma_clipping(data, low_rej=2.5, high_rej=2.5):
    # lower = - low_rej * data.std()
    # upper = + high_rej * data.std()
    # Just the good data (it masks vector outside the interval -lower and upper)
    # good_data = np.ma.masked_outside(data, lower, upper)
    # Just the clipped data
    # bad_data = np.ma.masked_inside(data, lower, upper)

    good_data = sigma_clip(data, sigma_lower=low_rej, sigma_upper=high_rej, iters=1, cenfunc=ma.mean)
    bad_data = good_data.mask * data
    bad_data = ma.masked_where(bad_data == 0.0, bad_data)

    # Creating masks
    good_mask = good_data / good_data
    bad_mask = bad_data / bad_data

    return good_mask, bad_mask


def cor_flat(flatname, function='Chebychev', order=55, clipping='sigma_clipping', niter=1,
             low_rej=5.0, high_rej=2.5):
    # Input and trimming data
    ccddata, hdr = fits.getdata(flatname, header=True, ignore_missing_end=True)
    axis0, axis1 = np.size(ccddata, axis=0), np.size(ccddata, axis=1)

    ccdfit = np.empty([axis0, axis1])
    # ccdfit = np.zeros(axis1)

    # Collapsing 2D image into a single vector (a big aperture of all lines)
    xpix = np.arange(0, np.size(ccddata, axis=1), 1)
    for j in np.arange(0, np.size(ccddata, axis=0)):
    #for j in np.arange(977, 1300, 1):

        spec = ccddata[j, :]

        # Fitting function to the spectrum of the flat before any clipping
        # fit, chi2_dof = fitting(xpix, spec, function=function, order=order)
        chi2_vet = np.empty(order)
        fit_vet = np.empty([order, len(spec)])
        for i in np.arange(0, order, 1):
            fit, chi2_dof = fitting(xpix, spec, function=function, order=i)
            fit_vet[i, :] = fit
            chi2_vet[i] = chi2_dof

        # Considering min of chi2
        index = int(np.mean(np.where(chi2_vet == min(chi2_vet))))
        order = index + 1
        fit = fit_vet[index, :]

        Plot = False
        if Plot is True:
            xorder = np.arange(1, len(chi2_vet)+1, 1)
            plt.plot(xorder, chi2_vet, label='Res')
            plt.legend(loc='best')
            plt.show()
        # Residual
        residual = spec - fit

        # Buffering quantities before clipping
        firstspec = spec
        firstresidual = residual

        bad_mask = ma.asanyarray(np.empty(len(residual)))

        if clipping.lower() == 'sigma_clipping' and niter > 0:

            niter += 1
            for iters in np.arange(1, niter, 1):
                good_mask, bad_mask = sigma_clipping(residual, low_rej=low_rej, high_rej=high_rej)

                # Filtering spec for good data
                spec = spec * good_mask

                # Fitting after  clipping iteration
                fit, chi2_dof = fitting(xpix, spec, function=function, order=order)

                # print fit - firstfit
                residual = spec - fit
        print 'Line: %s' % (j + 1)

        # newrow = np.asarray(fit)
        # ccdfit = np.vstack([ccdfit, newrow])

        ccdfit[j, :] = fit

        # Plotting things
        plot = False
        if plot is True:
            plt.figure()
            plt.clf()

            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1], sharex=ax0)

            firstresidual = 100 * (firstresidual / firstspec)
            residual = 100 * (residual / spec)

            if clipping.lower() == 'sigma_clipping' and niter > 0:
                ax0.plot(xpix, firstspec * bad_mask, 'rx', label='Clipped', lw=1)
                ax1.plot(xpix, firstresidual * bad_mask, 'rx', alpha=0.5, lw=2)

            ax0.plot(xpix, firstspec, 'k-', label='Flat Spectrum', lw=2)
            ax0.plot(xpix, fit, 'g-', label=str(function.title()) + ' - Order ' + str(order))
            ax1.plot(xpix, residual, '-', color='grey', alpha=0.5, lw=2)

            # Labelling things
            ax0.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=12)
            ax0.set_title('Chi2 = ' + str(chi2_dof))
            ax0.set_ylabel("ADU's", fontsize=12)
            ax1.set_xlabel('Column [pixel]', fontsize=12)
            ax1.set_ylabel('residual [%]', fontsize=11)

            ax0.set_xlim(xpix[0] - 150, xpix[-1] + 150)
            ax1.set_ylim(residual.min() - 4 * residual.std(), residual.max() + 4 * residual.std())

            plt.show()

    # ccdfit = ccdfit[1:,:]

    dir = '/home/davidsanm/PyCharmProjects/GoodmanDataReduction/2016-03-20/RED/TST/'
    fits.writeto(dir + 'c_master_flat_600.fits', ccdfit, hdr, clobber=True)


if __name__ == '__main__':

    # flat = '/home/davidsanm/PyCharmProjects/GoodmanDataReduction/2016-03-20/RED/TST/master_flat_600.fits'
    # flat_nogrt = '/home/davidsanm/PyCharmProjects/GoodmanDataReduction/2016-03-20/RED/master_flat_nogrt.fits'
    # cor_flat(flat, function='Legendre', order=41, clipping='sigma_clipping', niter=1, low_rej=2.5, high_rej=2.5)

    flat = '/home/davidsanm/PyCharmProjects/GoodmanDataReduction/2016-03-20/RED/TST/master_flat_600.fits'
    function = 'Legendre'
    order = 35
    clipping = 'sigma_clipping'
    niter = 1
    low_rej = 2.5
    high_rej = 1.5
    input_params = flat

    try:
        pool = Pool(processes=multiprocessing.cpu_count())
        result = pool.map_async(cor_flat(flat))
    finally:
        pool.close()
        pool.join()

'''
1) With pool.map_async
real	1m15.996s
user	1m23.126s
sys	1m9.093s

2) With pool.map

+ 5m

3) Without multiprocessing
real	4m37.276s
user	8m1.631s
sys	27m15.178s

'''
