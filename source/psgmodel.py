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
    def make_run_config_file(base_config_file_name, **kwargs):
        """ PSG generates a model with delta-wave = wave1 / srp.  However, the spectrum starts at wave1 + delta-wave
        and only the first 1 or 2 (?) decimal places of wave1/2 are read, so we .
        so we
        """
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
        # config_file_path = PsgModel.run_config_file
        # config_path = data_dir + PsgModel.run_config_file
        base_config_path = data_dir + base_config_file_name
        run_config_path = data_dir + PsgModel.run_config_file
        shutil.copy(base_config_path, run_config_path)
        new_lines = []
        with open(run_config_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                new_line = line
                for key in config_symbol:
                    symbol, val = config_symbol[key]
                    if symbol in line:
                        new_line = symbol + val + '\n'
                new_lines.append(new_line)
        file.close()

        os.remove(run_config_path)
        with open(run_config_path, 'w') as file:
            for new_line in new_lines:
                file.write(new_line)
        file.close()
        return

    @staticmethod
    def run_model():
        """ Run PSG and save spectrum to pickle file, labelled using v_obs.
        """
        os.chdir('./psg_in')
        curl_command = 'curl --data-urlencode file@config.txt https://psg.gsfc.nasa.gov/api.php'
        process = subprocess.Popen(curl_command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        os.chdir('..')
        text_block = output.decode('utf-8')
        wav_list, rad_list, rad_err_list = [], [], []
        lines = text_block.split('\n')
        for line in lines:
            tokens = line.split(' ')
            if '#' in tokens[0]:
                continue
            if len(tokens) < 2:
                break
            # print(tokens)
            wav, stellar, exo = float(tokens[0]), float(tokens[4]), float(tokens[6])
            wav_list.append(wav)
            rad_list.append(exo)
            rad_err_list.append(0.)

        waves, vals, val_errs = np.array(wav_list), np.array(rad_list), np.array(rad_err_list)
        units = 'Wsrm2um'
        label = 'Exoplanet'
        colour = 'magenta'
        rad_spec = Spectrum(waves, vals, val_errs, units, label, colour)
        return rad_spec

    @staticmethod
    def store_spectrum(spectrum, v_obs):
        pickle_path = PsgModel._make_pickle_path(v_obs)
        file = open(pickle_path, 'wb')
        pickle.dump(spectrum, file)
        file.close()
        return

    @staticmethod
    def load_spectrum(v_obs, **kwargs):
        pickle_path = PsgModel._make_pickle_path(v_obs)
        file = open(pickle_path, 'rb')
        spectrum = pickle.load(file)
        file.close()
        print("Loaded model spectrum with units of {:s}".format(spectrum.units))
        return spectrum

    @staticmethod
    def _make_pickle_path(v_obs):
        v_mps = int(1000. * v_obs)
        pickle_file_name = "pickle_v_{:06d}_mps".format(v_mps)       # Encode filename with v_obs in metre/sec
        pickle_path = './psg_out/' + pickle_file_name
        return pickle_path
