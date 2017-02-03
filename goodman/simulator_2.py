from astropy.io import fits
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
import time


class DataSimulator(object):

    def __init__(self, data_dimensions=1, x_size=None, y_size=None):
        self.dimensions = data_dimensions
        if self.dimensions == 1:
            self.x_length = x_size
            self.y_length = y_size
        elif self.dimensions == 2 and (x_size is not None and y_size is not None):
            self.x_length = x_size
            self.y_length = y_size
        else:
            raise NotImplementedError

    def __call__(self, *args, **kwargs):
        if self.dimensions == 1:
            if self.x_length is None and self.y_length is not None:
                print('object detection')
                zero_array = np.zeros(self.y_length)
            elif self.y_length is None:
                pass


def intensity(axis):
    return np.exp(axis/100.) *(3 * axis ** 2 + 5 * axis + 20000) + axis ** 3 + 1000 * np.sin(axis)
    # return 20000  * abs(np.sin(axis) +np.exp(axis) *np.sin(axis + 0.75 * np.pi)) +


def make_2d_spectra(n_targets=1, fwhm=8., snr=None, noise_level=0., plots=False):
    x = 4056
    y = 1550

    header_copy = fits.getheader('/data/simon/data/soar/work/20161114_eng/reduced_data/fzh.0298_CVSO166_400m2_gg455.fits')
    header_copy['OBJECT'] = 'SNR-%s'%str(snr)
    header_copy['HISTORY'] = 'Simulated spectrum Signal To Noise Ratio %s' % (snr)

    image = []
    if n_targets > 1:
        target_location = [y / 2., y / 2.]
    else:
        target_location = [y / 2.]

    sub_x = np.linspace(0, 10, x)
    y_axis = range(y)

    spectrum = [[] for i in range(int(n_targets))]
    # signal_peak = snr * noise_amplitude
    signal_peak = 7000.
    noise_amplitude = signal_peak / snr

    print('Signal Peak %.2f' % signal_peak)
    print('Noise Amplitude %.2f' % noise_amplitude)


    for i in range(x):
        if noise_level == -1:
            data = np.ones(y)
        else:
            data = np.random.normal(noise_level, noise_amplitude, y)
            # print(np.mean(data), np.std(data))
        for tar in range(int(n_targets)):
            # amplitude = intens * intensity(sub_x[i])

            # gauss = models.Gaussian1D(amplitude=amplitude, mean=target_location[tar], stddev=fwhm) 8.24687326842
            voigt = models.Voigt1D(amplitude_L=signal_peak, x_0=target_location[tar], fwhm_L=0.942561669206, fwhm_G=fwhm)
            # gauss2 = models.Gaussian1D(amplitude=amplitude, mean=target_location[1], stddev=fwhm)


            # spectrum[tar].append(amplitude)
            sd = voigt(y_axis)
            spectrum[tar].append(np.sum(sd[int(target_location[tar] - 2.5 * fwhm):int(target_location[tar] + 2.5 * fwhm)]))
            # gauss.amplitude.value = amplitude
            data += voigt(y_axis)
            # plt.plot(data)
            # plt.plot(voigt(y_axis) + noise_level, color='r')
            # plt.show()

        # signal_peak = xxxx
        # snr = 5
        # noise_level = 400.
        # noise_amplitude = signal_peak / snr
        # data =  np.random.normal(noise_level, noise_amplitude, data.shape)

        # data2 = gauss2(y_axis)
        # plt.plot(data)
        # plt.plot(data2)
        # plt.show()

        if i == int(x / 2.) and plots:
            # plt.title('FWHM Separation: %s' % separation_fwhm)
            plt.title('SNR: %s' % (snr))
            # plt.axvline(int(target_location[tar] - 2.5 * fwhm), color='r')
            # plt.axvline(int(target_location[tar] + 2.5 * fwhm), color='r')

            plt.plot(data)
            plt.plot(voigt(y_axis), color='r')
            plt.show()

        image.append(data)

        # plt.plot(y_axis, data)
        # plt.show()
        # rotated = zip(*original[::-1])
    rotated_image = zip(*image[::-1])
    hdu_name_file = '20161221_single-object_snr_%s.fits'% (snr)
    print(hdu_name_file)
    new_hdu = fits.PrimaryHDU(rotated_image, header=header_copy)
    new_hdu.writeto(hdu_name_file, clobber=True)
    for part_index in range(len(spectrum)):
        # f = open('obj.save', 'wb')
        # cPickle.dump(my_obj, f, protocol=cPickle.HIGHEST_PROTOCOL)
        # f.close()
        f = open(hdu_name_file.replace('.fits','_%s.pkl' % str(part_index + 1)), 'wb')
        pickle.dump(spectrum[part_index][::-1], f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()
        plt.plot(spectrum[part_index][::-1])
        plt.title('Simulated spectrum')
        plt.xlabel('Pixel Axis')
        plt.ylabel('Peak Intensity')
        plt.savefig('img/' + hdu_name_file.replace('.fits','.png'), dpi=300)
        if plots:
            plt.show()
        else:
            plt.clf()

    if True:
        median_rotated_image = np.median(rotated_image, axis=1)
        fig1, (ax1, ax2) = plt.subplots(1,2, gridspec_kw={'width_ratios':[4, 1]}, sharey=True, figsize=(25,10))
        # fig1.tight_layout()gridspec_kw = {'width_ratios':[3, 1]}
        # ax1 = plt.subplot(1, 3, gridspec_kw={'width_ratios':[3, 1]})
        ax1.set_title('Signal To Noise %s - n_targets %s' % (snr, n_targets))
        ax1.set_ylabel('Spatial Axis')
        ax1.set_xlabel('Dispersion Axis')
        ax1.imshow(rotated_image, clim=(5, 150), cmap='cubehelix', origin='lower', interpolation='nearest')

        # ax2 = plt.subplot(133)
        ax2.set_title('Median across Dispersion Axis')
        ax2.plot(median_rotated_image, range(len(median_rotated_image)))
        ax2.set_xlabel('Intensity (counts)')
        fig1.tight_layout()
        fig1.show()
        # fig1.clear()
        fig1.savefig('img/' + hdu_name_file.replace('.fits', '_im.png'), dpi=300)
        #plt.show()

        if plots:
            fig1.show()
        else:
            fig1.clf()
        plt.close(fig1)



if __name__ == '__main__':
    snr = [i for i in range(5, 20)]
    #for sep in range(1,15):
    for s in snr:
        print('SNR %s' % s)
        make_2d_spectra(n_targets=1, fwhm=8.24687326842, snr=s, noise_level=10., plots=False)
