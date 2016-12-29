import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d


# Fitting Sline3


def fit_spline3(y, x, order=3, nsum=3):


    y_resampled = [np.median(y[i:i + nsum]) for i in range(0, len(y) - len(y) % nsum, nsum)]
    x_resampled = np.linspace(0, len(y), len(y_resampled))

    # Fitting
    f = interp1d(x_resampled, y_resampled, kind=order, bounds_error=True)

    # Return function to be constructed with any other x array
    return f


# Local Minima and Maxima
def local_minmax(data, nmin=2, nmax=2):
    # Identifying indices of local minima-maxima points
    id_min = (np.gradient(np.sign(np.gradient(data))) > 0).nonzero()[0]  # index of local min
    id_max = (np.gradient(np.sign(np.gradient(data))) < 0).nonzero()[0]  # index of local max

    # Taking values at min/max points
    list_min, list_max = data[id_min], data[id_max]

    # Sorting minima-maxima values (bigger --> lower)
    list_min, id_min = (list(p) for p in zip(*sorted(zip(list_min, id_min), reverse=False)))
    list_max, id_max = (list(p) for p in zip(*sorted(zip(list_max, id_max), reverse=True)))

    # Taking the desired number of local minima-maxima points
    list_min, list_max, id_min, id_max = list_min[0:nmin], list_max[0:nmax], id_min[0:nmin], id_max[0:nmax]

    return list_min, list_max, id_min, id_max


def trim_slitedge(flat, plot=True):
    # Getting input data
    ccddata = fits.getdata(flat, ignore_missing_end=True)

    # Collapse flat in the dispersion direction
    flat_collapsed = fits.getdata(flat, ignore_missing_end=True).sum(axis=1) / ccddata.shape[1]
    lines = np.arange(0, flat_collapsed.size, 1)

    # Excluding first pixels in the spatial direction
    cut = 3
    c_flat = flat_collapsed[cut:-cut]
    c_lines = np.arange(0, c_flat.size, 1)

    # Fittin cubic spline. It's working very well with order=5, nsum=2
    func_splin3 = fit_spline3(c_flat, c_lines, order=5, nsum=2)
    smooth_flat = func_splin3(c_lines)

    # Compute 1st and flat smoothed
    dy = np.gradient(smooth_flat)
    dy2 = np.gradient(dy)

    # Regions to compute local minina-maxima
    # Region one: it represent first 40 percent of all data
    # Region two: ... last 40%
    pixa, pixb = int(len(c_flat) * 0.4), int(len(c_flat) * 0.6)
    dy2_one, dy2_two = dy2[0:pixa], dy2[pixb:]

    # Reg. 1: Compute local min/max of the 2nd derivative
    list_min_1, list_max_1, id_min_1, id_max_1 = local_minmax(dy2_one, nmin=1, nmax=1)
    list_min_2, list_max_2, id_min_2, id_max_2 = local_minmax(dy2_two, nmin=1, nmax=1)

    # Indice have to be reshifted to the original indices of the function dy2
    id_min_2 = np.array(id_min_2) + pixb

    # Slit edges are the local maxima/minima 1/2 [accounting the cutted pixels]
    slit_1, slit_2 = int(np.array(id_min_1) + cut), int(np.array(id_min_2) + cut)

    print slit_1, slit_2

    if plot is True:
        import matplotlib.pyplot as plt
        c_lines += cut
        plt.plot(lines, flat_collapsed, 'k-', label='Flat Collapsed')
        plt.plot(lines[slit_1:slit_2], flat_collapsed[slit_1:slit_2], 'r-', label = 'Cutted Flat')
        plt.plot(c_lines, dy, 'g-', label="Dy/dx")
        plt.plot(c_lines, dy2, 'y-', label="Dy2/dx")
        plt.plot(slit_1, list_min_1, 'bo', label='Slit Edge 1 ')
        plt.plot(slit_2, list_min_2, 'ro', label='Slit Edge 2')
        plt.xlim(lines.min() - 50, lines.max() + 50)
        plt.legend(loc='best')
        plt.show()

    return slit_1, slit_2

flat = '/home/davidsanm/PyCharmProjects/GoodmanDataReduction/2016-03-20/RED/master_flat_600.fits'
trim_slitedge(flat, plot = True)