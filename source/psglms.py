import math
import numpy as np
from numpy import random
from psgspectrum import PsgSpectrum as Spectrum


class PsgLms:

    flux_spectrum = None, None, '-'
    det_spectrum = None, None, '-'
    d_elt = 39.
    read_interval = 0.84        # Read interval (seconds)
    # Read noise as function of n_reads (2x values quoted for best detector in Finger et al. 2008)
    read_noise = {2: 16., 4: 12., 8: 9., 16: 7., 32: 5.6, 64: 5.,
                  128: 5., 256: 5.6, 512: 6.4, 1024: 8., 2048: 10.6}
    nreads_v_read_noise = np.array([[2, 16.], [4, 12.], [8, 9.],
                                    [16, 7.], [32, 5.6], [64, 5.],
                                    [128, 5.], [256, 5.6], [512, 6.4],
                                    [1024, 8.], [2048, 10.6]])

    def __init__(self):
        return

    @staticmethod
    def configure(**kwargs):
        PsgLms.dit = kwargs.get('dit', 100.)
        PsgLms.pce = kwargs.get('pce', 0.1)
        PsgLms.rnoise = kwargs.get('rnoise', 10.)
        PsgLms.n_pix_aper = kwargs.get('n_pix_aper', 6)
        return

    def detect(self, flux_spectrum, **kwargs):
        """ Perform a simulated detection on a spectrum.  If apply_noise = True the spectrum will have a normal
        noise distribution derived from the error array added to it.
        """
        apply_noise = kwargs.get('apply_noise', False)
        area_elt = math.pi * math.pow(PsgLms.d_elt / 2., 2.)
        solid_angle_pix = 0.018 * 0.018
        waves = flux_spectrum.waves
        flux = flux_spectrum.vals
        bw = np.zeros(waves.shape)
        bw[:-1] = np.array(waves[1:] - waves[:-1])
        bw[-1] = bw[-2]
        n_pix = PsgLms.n_pix_aper
        sa_pix = n_pix * solid_angle_pix
        rad_pix = flux * area_elt * bw * sa_pix             # ph/sec/pixel
        sig_pix = rad_pix * PsgLms.pce * PsgLms.dit       # el/pixel
        noise_pix = np.sqrt(sig_pix) + PsgLms.rnoise
        signal = sig_pix * n_pix
        error = noise_pix * math.sqrt(n_pix)
        rng = np.random.default_rng()
        shape = signal.shape
        zero, sigma = np.zeros(shape), error
        noise = rng.normal(zero, sigma, shape)
        if apply_noise:
            signal = signal + noise

        det_units = 'el'
        det_spectrum = Spectrum(waves, signal, error, det_units, 'Signal', 'red')
        return det_spectrum
