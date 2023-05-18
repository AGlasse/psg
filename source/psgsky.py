import numpy as np
from astropy.io import fits
from psgspectrum import PsgSpectrum as Spectrum


class PsgSky:

    trans_spectrum, flux_spectrum = None, None

    def __init__(self):
        return

    def load_spectra(self):
        path = '../data/elt_sky.fits'
        hdu_list = fits.open(path, mode='readonly')
        data_table = hdu_list[1].data
        wave = data_table['lam'] / 1000.
        trans = data_table['trans']
        trans_units = '-'
        PsgSky.trans_spectrum = Spectrum(wave, trans, None, '-', 'Transmission', 'black')
        flux = data_table['flux']
        flux_units = 'ph/s/m2/um/arcsec2'
        PsgSky.flux_spectrum = Spectrum(wave, flux, None, flux_units, 'Flux', 'blue')
        print("Loaded sky transmission and emission spectrum with units {:s}".format(flux_units))
        return

    def get_spectrum(self, ord):
        if ord == 'flux':
            return PsgSky.flux_spectrum
        else:
            return PsgSky.trans_spectrum

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

