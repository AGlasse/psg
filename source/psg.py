from psgsky import PsgSky as Sky
from psgplot import PsgPlot as Plot
from psgmodel import PsgModel as Model
from psglms import PsgLms as Lms
from psgglobals import PsgGlobals as Globals
import numpy as np

if __name__ == "__main__":
    print('psg starting')
    lms = Lms()
    model = Model()
    globals = Globals()
    plot = Plot()

    # Define PSG model parameters
    reference_wavelength = 2.7
    srp = 200000                        # SRP at reference wavelength, with dw = wave_range[0] / srp (nom 200,000)
    spec_res = reference_wavelength / srp
    model_wave_range = [3.0, 3.2]       # Extract limited range in a single run (to avoid truncation by PSG)
    v_obs = 0.                          # System velocity km/sec

    plot_wave_range = [3.00, 3.05]
    plot.set_wave_range(plot_wave_range)

    # Load sky transmission and emission spectra
    sky = Sky()
    sky.load_spectra()
    w_sky = sky.trans_spectrum.waves
    n_w = w_sky.shape[0]
    w_idx = 2.7 + 0.0000135 * np.arange(0, n_w, 1, dtype='int')
    dw = 1E7 * (w_sky - w_idx)

    run_psg = True
    if run_psg:
        # Adjust the model wavelength range to align with the sky spectrum.
        sky_waves = sky.trans_spectrum.waves
        sky_hi_indices = np.where(sky_waves > model_wave_range[1])
        sky_lo_indices = np.where(sky_waves < model_wave_range[0])
        lo_idx, hi_idx = sky_lo_indices[0][-1], sky_hi_indices[0][0]
        mod_wav_lo, mod_wav_hi = sky_waves[lo_idx], sky_waves[hi_idx]

        base_config_file = 'ross458c_config.txt'
        kwargs = {'wave1': "{:.9f}".format(mod_wav_lo),     # PSG first point is at wave1 + dw..
                  'wave2': "{:.9f}".format(mod_wav_hi),
                  'spec_res': "{:.7f}".format(spec_res),
                  'v_obs': "{:.2f}".format(v_obs)
                  }
        model.make_run_config_file(base_config_file, **kwargs)
        model.run_model()

    pickle_path = '../output/pickle.pkl'
    rad_spectrum = model.load_spectrum(pickle_path, lms=lms)
    plot.spectra([rad_spectrum])
    wave_range = model.flux_spectrum.get_wave_range()

    mod_flux_spec, wave_range = Globals.trim_wavelengths(model.flux_spectrum, wave_range)

    # Down sample the model spectrum onto the sky transmission wavelength grid
    sky_trans_spec, wave_range = Globals.trim_wavelengths(sky.trans_spectrum, wave_range)
    sky_flux_spec, wave_range = Globals.trim_wavelengths(sky.flux_spectrum, wave_range)
    plot.spectra([sky_trans_spec], title='Sky transmission')
    plot.spectra([sky_flux_spec], title='Sky (plus telescope?) emission')
    ds_mod_flux_spec = Globals.down_sample(mod_flux_spec, sky_trans_spec)
    plot.spectra([mod_flux_spec], ymin=0., title='PSG model emission')
    lms.configure(dit=100., pce=0.1, rnoise=10.)
    lms_tgt_spectrum, lms_bgd_spectrum = sky.rad_trans(ds_mod_flux_spec, sky_trans_spec, sky_flux_spec)
#    plot.spectra([lms_tgt_spectrum, lms_bgd_spectrum])
    ideal_tgt_spectrum = lms.detect(lms_tgt_spectrum, apply_noise=False)
    ideal_bgd_spectrum = lms.detect(lms_bgd_spectrum, apply_noise=False)
    ideal_spectrum = ideal_tgt_spectrum.minus(ideal_bgd_spectrum)
    det_tgt_spectrum = lms.detect(lms_tgt_spectrum, apply_noise=True)
    plot.spectra([det_tgt_spectrum],
                 title='Target + background signal', plot_errors=[True])
    det_bgd_spectrum = lms.detect(lms_bgd_spectrum, apply_noise=True)
    det_spectrum = det_tgt_spectrum.minus(det_bgd_spectrum)
    det_spectrum.colour = 'red'
    ideal_spectrum.colour = 'blue'
    plot.spectra([det_spectrum, ideal_spectrum],
                 title='Background subtracted target', plot_errors=[False, True])

    print('Done')
