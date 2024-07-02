import numpy as np
from astropy.io import fits
from psgspectrum import PsgSpectrum as Spectrum


class PsgSky:

    trans_spectrum, flux_spectrum = None, None

    def __init__(self):
        return

    def resample(self, new_waves):
        """ Resample the transmission and flux spectra onto a new wavelength scale. This is done by interpolation.
        The new scale should/must be smaller in extent than the current scale and increase monotonically.
        New spectra objects are returned (to allow pre-/-post comparison).
        """
        waves = PsgSky.trans_spectrum.waves
        t_vals = PsgSky.trans_spectrum.vals
        f_vals = PsgSky.flux_spectrum.vals
        nt_vals, nf_vals = np.zeros(new_waves.shape), np.zeros(new_waves.shape)
        i = 0
        fmt = "{:10s}{:10s}{:10s}{:10s}{:10s}{:10s}"
        print(fmt.format('New Wave', 'New T', 'W1', 'W2', 'T1', 'T2'))
        for j, new_wave in enumerate(new_waves):
            while waves[i] < new_wave:
                i += 1
            nt_vals[j] = np.interp(new_wave, waves[i-1:i+1], t_vals[i-1:i+1])
            nf_vals[j] = np.interp(new_wave, waves[i-1:i+1], f_vals[i-1:i+1])
            i += 1
        new_trans_spectrum = Spectrum(self.trans_spectrum)
        new_trans_spectrum.waves, new_trans_spectrum.vals = new_waves, nt_vals
        new_flux_spectrum = Spectrum(self.flux_spectrum)
        new_flux_spectrum.waves, new_flux_spectrum.vals = new_waves, nf_vals
        PsgSky.trans_spectrum = new_trans_spectrum
        PsgSky.flux_spectrum = new_flux_spectrum
        return

    @staticmethod
    def load_spectra():
        path = './psg_in/elt_sky.fits'
        hdu_list = fits.open(path, mode='readonly')
        data_table = hdu_list[1].data
        wave = data_table['lam'] / 1000.
        trans, trans_errs = data_table['trans'], None
        units, label, colour = '-', 'Transmission', 'black'
        PsgSky.trans_spectrum = Spectrum(wave, trans, trans_errs, units, label, colour)
        flux, flux_errs = data_table['flux'], None
        flux_units, label, colour = 'ph/s/m2/um/arcsec2', 'Emission', 'orange'
        PsgSky.flux_spectrum = Spectrum(wave, flux, flux_errs, flux_units, label, colour)
        print("Loaded sky transmission and emission spectrum with units {:s}".format(flux_units))
        return

    @staticmethod
    def rad_trans(model_spectrum, sky_trans_spec, sky_flux_spec):
        sky_wav, sky_flux, sky_flux_err, sky_units, _, _ = sky_flux_spec.get()
        _, sky_trans, _, _, _, _ = sky_trans_spec.get()
        exo_wav, exo_flux, exo_flux_err, exo_units, _, _ = model_spectrum.get()
        tgt_flux = sky_trans * exo_flux + sky_flux
        tgt_flux_errs = np.zeros(tgt_flux.shape)
        bgd_flux = sky_flux
        bgd_flux_errs = np.zeros(bgd_flux.shape)
        tgt_spectrum = Spectrum(exo_wav, tgt_flux, tgt_flux_errs, exo_units, 'LMS tgt flux', 'green')
        bgd_spectrum = Spectrum(exo_wav, bgd_flux, bgd_flux_errs, exo_units, 'LMS bgd flux', 'grey')
        return tgt_spectrum, bgd_spectrum
