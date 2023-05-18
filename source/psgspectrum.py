import math
import numpy as np


class PsgSpectrum:

    def __init__(self, *args, **kwargs):
        n_args = len(args)
        if n_args > 0 and isinstance(args[0], PsgSpectrum):
            sp = args[0]
            self.waves, self.vals, self.val_errs, self.units, self.label, self.colour = \
                sp.waves, sp.vals, sp.val_errs, sp.units, sp.label, sp.colour
            return
        # Create default parameters
        self.waves, self.vals, self.val_errs = np.zeros(2), np.array([1., 1.]), np.array([.1, .1])
        self.units, self.label, self.colour =  '-', 'default', 'black'
        if n_args > 3:
            self.waves, self.vals, self.units = args[0], args[1], args[3]
            self.val_errs = np.zeros(self.vals.shape) if args[2] is None else args[2]
        if n_args > 5:
            self.label, self.colour = args[4], args[5]
        return

    def get(self):
        return self.waves, self.vals, self.val_errs, self.units, self.label, self.colour

    def set(self, waves, vals, val_errs, units, label, colour):
        self.waves, self.vals, self.val_errs, self.units = waves, vals, val_err, units
        self.label, self.colour = label, colour
        return

    def get_wave_range(self):
        wmin = np.amin(self.waves)
        wmax = np.amax(self.waves)
        return [wmin, wmax]

    def minus(self, spec2):
        """ Subtract a spectrum from this one, including noise propagation.  Return as a new spectrum
        """
        diff_spec = PsgSpectrum(self)
        diff_spec.vals = diff_spec.vals - spec2.vals
        diff_vars = np.square(diff_spec.val_errs)
        spec2_vars = np.square(spec2.val_errs)
        diff_spec.val_errs = np.sqrt(diff_vars + spec2_vars)
        return diff_spec

