from astropy.io import fits
from astropy.modeling import models, fitting
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle


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


def make_2d_spectra(separation_fwhm=0, n_targets=1, fwhm=8., intens=0., noise_level=1., plots=False):
    x = 4056
    y = 1550

    header_copy = fits.getheader('/data/simon/data/soar/work/20161114_eng/reduced_data/fzh.0298_CVSO166_400m2_gg455.fits')
    header_copy['OBJECT'] = 'Test-%s'%str(separation_fwhm)
    header_copy['HISTORY'] = 'Simulated spectrum N-sources %s separation_fwhm %s FWHM %s' % (n_targets, separation_fwhm, fwhm)

    targets = n_targets
    target_separation_fwhm = float(separation_fwhm)
    image = []
    if targets > 1:
        target_location = [y / 2. - target_separation_fwhm / float(targets) * fwhm,
                           y / 2. + target_separation_fwhm / float(targets) * fwhm]
        print separation_fwhm, int(y / 2.), target_location, (target_separation_fwhm / targets) * fwhm
    else:
        target_location = [y / 2.]

    sub_x = np.linspace(0, 10, x)
    y_axis = range(y)

    spectrum = [[] for i in range(int(targets))]

    for i in range(x):
        if noise_level == 0:
            data = np.ones(y)
        else:
            data = np.random.normal(10, noise_level, y)
        for tar in range(int(targets)):
            amplitude = intens * intensity(sub_x[i])
            # theo_snr = amplitude / noise_level
            # print theo_snr
            # gauss = models.Gaussian1D(amplitude=amplitude, mean=target_location[tar], stddev=fwhm) 8.24687326842
            voigt = models.Voigt1D(amplitude_L=amplitude, x_0=target_location[tar], fwhm_L=0.942561669206, fwhm_G=fwhm)
            # gauss2 = models.Gaussian1D(amplitude=amplitude, mean=target_location[1], stddev=fwhm)

            # spectrum[tar].append(amplitude)
            sd = voigt(y_axis)
            spectrum[tar].append(np.sum(sd[int(target_location[tar] - 2.5 * fwhm):int(target_location[tar] + 2.5 * fwhm)]))
            # gauss.amplitude.value = amplitude
            data += voigt(y_axis)

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
            plt.title('Intensity: %s' % (intens))
            plt.axvline(int(target_location[tar] - 2.5 * fwhm), color='r')
            plt.axvline(int(target_location[tar] + 2.5 * fwhm), color='r')
            plt.plot(data)
            plt.show()

        image.append(data)

        # plt.plot(y_axis, data)
        # plt.show()
        # rotated = zip(*original[::-1])
    rotated_image = zip(*image[::-1])
    hdu_name_file = '20161128_single-object_n%s_s%s-fwhm_%1.3f_int.fits'% (str(int(n_targets)), str(int(separation_fwhm)), intens)
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

    if plots:
        plt.title('Target Separation %s - N targets %s' % (str(separation_fwhm), targets))
        plt.imshow(rotated_image)
        plt.show()



if __name__ == '__main__':
    inten = np.linspace(0.1, 1.3, 20)
    #for sep in range(1,15):
    for inte in inten:
        make_2d_spectra(n_targets=1, fwhm=8.24687326842, intens=inte, noise_level=3.5)
