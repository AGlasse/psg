import sys
import math
import numpy as np
from psg_spectrum import PSG_spectrum as Spectrum


class PSG_globals:

    data_dir, output_dir = '../data/', '../output/'


    def __init__(self):
        return

    @staticmethod
    def convert_flux_to_radiance(flux_spectrum):
        """ Convert a spectrum from sky flux units (ph/s/m2/micron/arcsec2 to PSG units watt/m2/sr/micron)
        """
        wav, flux, flux_units, label = flux_spectrum
        arcsec_steradian = 4.2545166E10
        h_planck = 6.6260693E-34
        c_light = 2.99792458E14
        joule_ph = h_planck * c_light / wav
        radiance = flux * joule_ph * arcsec_steradian
        return wav, radiance, 'W/sr/m2/um', label

    @staticmethod
    def convert_radiance_to_flux(radiance_spectrum):
        """ Convert a spectrum from sky flux units (ph/s/m2/um/arcsec2 to PSG units watt/m2/sr/micron)
        """
        wav, radiance, radiance_error, radiance_units, label, colour = radiance_spectrum.get()
        arcsec_steradian = 4.2545166E10
        h_planck = 6.6260693E-34
        c_light = 2.99792458E14
        joule_ph = h_planck * c_light / wav
        k = 1. / (joule_ph * arcsec_steradian)
        flux = k * radiance
        flux_error = k * radiance_error
        label = 'Flux'
        flux_spectrum = Spectrum(wav, flux, flux_error, 'ph/s/m2/um/arcsec2', label, colour)
        return flux_spectrum

    @staticmethod
    def down_sample(spec, template_spec):
        """ Create a new spectrum 'down_spec' holding those values in spec which are closest (< 1 pm) in wavelength to
        the wavelength points in 'template_spec'.
        """
        waves, vals, val_errs, units, label, colour = spec.get()
        twaves, _, _, _, _, _ = template_spec.get()
        dwaves, dvals, dval_errs = np.zeros(twaves.shape), np.zeros(twaves.shape), np.zeros(twaves.shape)
        ite, i = 0, 0
        for twav in twaves:
            for wav in waves[i:]:
                is_equal = math.isclose(wav, twav, rel_tol=1.e-7, abs_tol=1.e-7)
                if is_equal:
                    dwaves[ite], dvals[ite], dval_errs[ite] = wav, vals[i], val_errs[i]
                    break
                i += 1
            ite += 1
        down_spec = Spectrum(dwaves, dvals, dval_errs, units, label, colour)
        return down_spec

    @staticmethod
    def trim_wavelengths(sp_in, wave_range):
        """ Copy the passed spectrum and trim to include points within the passed wavelength range. """
        sp_out = Spectrum(sp_in)
        waves, vals, val_errs = sp_out.waves, sp_out.vals, sp_out.val_errs
        a = np.nonzero(waves >= wave_range[0])
        idx1 = np.nonzero(waves >= wave_range[0])[0][0]
        idx2 = np.nonzero(waves <= wave_range[1])[0][-1]
        waves_out = waves[idx1: idx2 + 1]
        vals_out = vals[idx1: idx2 + 1]
        val_errs_out = val_errs[idx1: idx2 + 1]

        sp_out.waves, sp_out.vals, sp_out.val_errs = waves_out, vals_out, val_errs_out
        wave_range_out = [waves_out[0], waves_out[-1]]
        return sp_out, wave_range_out
