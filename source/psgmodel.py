import subprocess
import pickle
import os
import shutil
import numpy as np
from psgglobals import PsgGlobals as Globals
from psgspectrum import PsgSpectrum as Spectrum


class PsgModel:

    rad_spectrum, flux_spectrum = None, None
    run_config_file = 'config.txt'

    def __init__(self):
        return

    @staticmethod
    def make_run_config_file(base_config_file, **kwargs):
        """ PSG generates a model with delta-wave = wave1 / srp.  However, the spectrum starts at wave1 + delta-wave
        and only the first 1 or 2 (?) decimal places of wave1/2 are read, so we .
        so we """

        wave1 = kwargs.get('wave1', '2.7')
        wave2 = kwargs.get('wave2', '2.9')
        spec_res = kwargs.get('spec_res', '0.0000135')
        spec_res_units = kwargs.get('spec_res_units', 'um')
        rad_units = kwargs.get('rad_units', 'Wsrm2um')
        v_obs = kwargs.get('v_obs', '-27.000')

        config_symbol = {'wave1': ('<GENERATOR-RANGE1>', wave1),
                         'wave2': ('<GENERATOR-RANGE2>', wave2),
                         'spec_res': ('<GENERATOR-RESOLUTION>', spec_res),
                         'spec_res_units': ('<GENERATOR-RESOLUTIONUNIT>', spec_res_units),
                         'rad_units': ('<GENERATOR-RADUNITS>', rad_units),
                         'v_obs': ('<OBJECT-OBS-VELOCITY>', v_obs)
                         }

        data_dir = Globals.data_dir
        config_file = PsgModel.run_config_file
        config_path = data_dir + PsgModel.run_config_file
        os.chdir(data_dir)
        shutil.copy(base_config_file, config_file)
        new_lines = []
        with open(config_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                new_line = line
                for key in config_symbol:
                    symbol, val = config_symbol[key]
                    if symbol in line:
                        new_line = symbol + val + '\n'
                new_lines.append(new_line)
        file.close()

        os.remove(config_path)
        with open(config_path, 'w') as file:
            for new_line in new_lines:
                file.write(new_line)
        file.close()
        os.chdir('../source')
        return

    @staticmethod
    def run_model():
        config_path = Globals.data_dir + 'config.txt'
        os.chdir('../data')
        curl_command = 'curl --data-urlencode file@config.txt https://psg.gsfc.nasa.gov/api.php'
        process = subprocess.Popen(curl_command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        data = output.decode('utf-8')
        pickle_path = '../output/pickle.pkl'
        file = open(pickle_path, 'wb')
        pickle.dump(data, file)
        file.close()
        os.chdir('../source')
        return

    @staticmethod
    def load_spectrum(pickle_path, **kwargs):
        lms = kwargs.get('lms', None)
        rad_spectrum = PsgModel._load_pickle(pickle_path)
        flux_spectrum = Globals.convert_radiance_to_flux(rad_spectrum)
        PsgModel.rad_spectrum, PsgModel.flux_spectrum = rad_spectrum, flux_spectrum
        print("Loaded model spectra with units of {:s}".format(rad_spectrum.units))
        return rad_spectrum

    @staticmethod
    def _load_pickle(pickle_path):
        file = open(pickle_path, 'rb')
        data = pickle.load(file)
        file.close()

        # Parse data file
        wav_list, rad_list, rad_err_list = [], [], []
        records = data.split('\n')
        print("Model flux units = {:s}".format(records[11]))
        for record in records[14:]:
            if len(record) < 2:
                continue
#            print(record)
            tokens = record.split(' ')
            floats = []
            for token in tokens:
                if len(token) > 2:
                    floats.append(float(token))
            wav, rad = floats[0], floats[1]
            wav_list.append(wav)
            rad_list.append(rad)
            rad_err_list.append(0.)
        rad_units = 'W/sr/m2/um'
        print("Loaded model radiance spectrum with units {:s}".format(rad_units))
        spectrum = Spectrum(np.array(wav_list), np.array(rad_list), np.array(rad_err_list), rad_units,
                            'Radiance', 'orange')
        return spectrum
