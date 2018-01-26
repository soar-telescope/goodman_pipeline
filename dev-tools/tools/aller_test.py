from ccdproc import (CCDData, ImageFileCollection)
from scipy import (signal, interpolate)
import os
import glob
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io.fits.header import Header

import sys
sys.path.append('/user/simon/development/soar/goodman')

# import goodman
from pipeline.core import (read_fits, write_fits, identify_targets, trace ,extraction)
from astropy.modeling import (models, fitting)
from pipeline.wcs.wcs import WCS

plt.rcParams["figure.figsize"] = [16, 9]


def normalize_data(data):
    return data / np.mean(data)


def aller_test(image_list):
    wcs = WCS()
    assert isinstance(wcs, WCS)
    all_data = []
    all_dates = []
    for file_name in image_list:
        ccd = CCDData.read(file_name, unit=u.adu)
        print("Master Flat used: {:s}".format(ccd.header['GSP_FLAT']))
        print("{:s} : {:s}".format(ccd.header['DATE-OBS'], file_name))
        # print(ccd.header["GSP_*"])
        wav, intens = wcs.read(ccd=ccd)
        all_data.append([wav, normalize_data(intens), ccd.header])
        all_dates.append(ccd.header['DATE'])
        plt.title(file_name)
        plt.plot(wav, normalize_data(intens))
        plt.xlabel("Wavelength (Angstrom)")
        plt.ylabel("Intensity (ADU)")
        plt.show()

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    plt.title("Normalized Spectrum")
    ax1.axvline(4822.3787, label='4822.3787A', color='k', linestyle='--',
                alpha=0.7)
    ax2.axvline(4955.91, label='4955.91A', color='k', linestyle='--', alpha=0.7)
    for i in range(len(all_data)):
        ax1.plot(all_data[i][0], all_data[i][1], label=all_dates[i])
        ax1.set_xlim(4817, 4826)
        # ax1.axvline(4822.3787, label='4822.3787A', color='k', linestyle='--')
        ax2.plot(all_data[i][0], all_data[i][1], label=all_dates[i])
        ax2.set_xlim(4951, 4958)

        ax1.set_xlabel("Wavelength (Angstrom)")
        ax1.set_ylabel("Intensity (ADU)")
        ax1.legend(loc='best')

    plt.legend(loc='best')
    plt.show()
    return all_data





def extract_from_file(image_list, clim=None):
    all_traces = []
    dispersion = 0
    for image in image_list:
        ccd = read_fits(image)
        spatial, dispersion = ccd.data.shape
        print(spatial, dispersion)

        target = identify_targets(ccd=ccd)
        print(target[0])
        trace_model = models.Polynomial1D(degree=2)
        trace_fitter = fitting.LevMarLSQFitter()

        traces = trace(ccd=ccd,
                       model=target[0],
                       trace_model=trace_model,
                       fitter=trace_fitter,
                       sampling_step=5)
        all_traces.append([traces, image])
        print(traces)
        extracted = extraction(ccd=ccd,
                               trace=traces,
                               spatial_profile=target[0],
                               extraction='fractional')
        # print(traces)


        print(np.median(ccd.data), np.mean(ccd.data))
        plt.title(image)
        plt.plot(traces(range(dispersion)), color='r')
        if clim is None:
            clim = (0.3 * np.median(ccd.data), 0.3 * np.mean(ccd.data))
        plt.imshow(ccd.data, clim=clim)
        plt.show()

        plt.plot(extracted.data)
        plt.show()
    for single_trace, image in all_traces:
        plt.plot(single_trace(range(dispersion)), label=image)
    plt.legend(loc='best')
    plt.title("Target Traces")
    plt.xlabel("Dispersion Axis")
    plt.ylabel("Spatial Axis")
    plt.show()


def plot_and_split(all_data):
    fig2, axes = plt.subplots(8, 2, figsize=(16, 30))
    limits = np.linspace(0, len(all_data[0][1]) - 1, 17)
    all_limits = [[int(limits[i]), int(limits[i + 1])] for i in
                  range(0, len(limits) - 1, 1)]
    print(all_limits)

    for data in all_data:
        z = 0
        for i in range(8):
            for e in range(2):
                # print(i, e)
                # print(i * e, i * e +1)
                # print(all_limits[z])
                axes[i, e].plot(data[0], data[1],
                                label="{:s}".format(data[2]['DATE']))
                axes[i, e].set_xlim((all_data[0][0][all_limits[z][0]],
                                     all_data[0][0][all_limits[z][1]]))
                axes[i, e].legend(loc='best')

                z += 1
    # for ax in sub_ax:
    #         ax.plot(all_datas[0])

    plt.show()


if __name__ == '__main__':
    path = '/user/simon/jupyter'
    images = glob.glob(os.path.join(path, "data/c*fits"))
    extract_from_file(image_list=images)

    files = glob.glob(os.path.join(path, 'data/g*HD*fits'))
    all_data = aller_test(image_list=files)

    lamp_files = glob.glob(os.path.join(path, "data/g*Cu*fits"))
    all_lamps = aller_test(lamp_files)

    plot_and_split(all_data=all_lamps)






