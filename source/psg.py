from psg_sky import PSG_sky as Sky
from psg_plot import PSG_Plot as Plot
from psg_model import PSG_model as Model
from psg_lms import PSG_lms as Lms
from psgglobals import PsgGlobals as Globals

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

    run_psg = True
    if run_psg:
        base_config_file = 'ross458c_config.txt'
        kwargs = {'wave1': "{:.9f}".format(model_wave_range[0]),     # PSG first point is at wave1 + dw..
                  'wave2': "{:.9f}".format(model_wave_range[1]),
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
