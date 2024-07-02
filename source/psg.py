from psgsky import PsgSky as Sky
from psgplot import PsgPlot as Plot
from psgmodel import PsgModel as Model
from psglms import PsgLms as Lms
from psgglobals import PsgGlobals as Globals
from psgcorrelate import PsgCorrelate as Correlate

if __name__ == "__main__":
    print('psg starting')
    Globals()
    lms = Lms()
    model = Model()
    plot = Plot()
    correlate = Correlate()

    # Define PSG model parameters
    reference_wavelength = 2.7
    srp = 200000                        # SRP at reference wavelength, with dw = wave_range[0] / srp (nom 200,000)
    spec_res = reference_wavelength / srp

    plot_wave_range = [3.15, 3.20]
    plot.set_wave_range(plot_wave_range)

    model_wave_range = [3.10, 3.25]          # Extract limited range in a single run (to avoid truncation by PSG)
    v_obs = 0.                              # System velocity km/sec
    v_obs_list = [0., 0.005, 5.000]         # Velocity shift in km/sec

    run_psg = True
    if run_psg:
        for v_obs in v_obs_list:
            base_config_file_name = 'hd189733b_config.txt'    # 'ross458c_config.txt'
            kwargs = {'wave1': "{:.9f}".format(model_wave_range[0] - spec_res),     # PSG first point is at wave1 + dw..
                      'wave2': "{:.9f}".format(model_wave_range[1] + spec_res),
                      'spec_res': "{:.7f}".format(spec_res),
                      'v_obs': "{:.2f}".format(v_obs)
                      }
            model.make_run_config_file(base_config_file_name, **kwargs)
            rad_spec = model.run_model()
            flux_spec = Globals.convert_radiance_to_flux(rad_spec)
            model.store_spectrum(flux_spec, v_obs)

    # Plot radiance change wrt first (0 km/sec) spectrum.  Illustrates that shifts <~ 5 m/s are not modelled by PSG.
    ref_spec = None
    sky = None
    obs_list = []
    for v_obs in v_obs_list:

        # Load sky transmission and emission spectra and resample them onto the exo-planet model wavelength scale
        if sky is None:
            exo_spec = model.load_spectrum(v_obs)
            title = "Noiseless target spectrum $\Delta$V = {:8.3f}".format(v_obs)
            plot.spectra([exo_spec], ymin=0., title=title)

            sky = Sky()
            sky.load_spectra()
            w_sky = sky.trans_spectrum.waves
            model_waves = exo_spec.waves
            sky.resample(model_waves)
            plot_sky = False
            if plot_sky:
                plot.spectra([sky.trans_spectrum], title='Sky transmission')
                plot.spectra([sky.flux_spectrum], title='Sky (plus telescope) emission')

        if ref_spec is None:
            ref_spec = exo_spec
        do_dv_check = False
        if do_dv_check:
            fig, plt, ax_list = plot.set_plot_area()
            ax = ax_list[0, 0]
            title = "Radiance difference for $\Delta$V = {:6.3f} km/s".format(v_obs)
            ax.set_title(title)
            ax.plot(ref_spec.vals, exo_spec.vals - ref_spec.vals, marker='.', ls='none', color='blue')
            plt.show()

        lms.configure(dit=100., pce=0.1, rnoise=10.)
        lms_tgt_spectrum, lms_bgd_spectrum = sky.rad_trans(exo_spec, sky.trans_spectrum, sky.flux_spectrum)
        lms_tgt_spectrum.colour, lms_bgd_spectrum.colour = 'brown', 'orange'
        # plot.spectra([lms_tgt_spectrum, lms_bgd_spectrum],
        #              title='Target + background flux at ELT', plot_errors=[True, False])

        ideal_obs_spectrum = lms.detect(lms_tgt_spectrum, apply_noise=False)
        ideal_bgd_spectrum = lms.detect(lms_bgd_spectrum, apply_noise=False)
        ideal_tgt_spectrum = ideal_obs_spectrum.minus(ideal_bgd_spectrum)
        det_tgt_spectrum = lms.detect(lms_tgt_spectrum, apply_noise=True)
        det_tgt_spectrum.colour = 'magenta'
        det_tgt_spectrum.label = 'target + bgd'

        det_bgd_spectrum = lms.detect(lms_bgd_spectrum, apply_noise=True)
        det_bgd_spectrum.colour = 'orange'
        det_bgd_spectrum.label = 'background'
        # plot.spectra([det_tgt_spectrum, det_bgd_spectrum],
        #              title='Target + background signal', plot_errors=[True, False])
        det_spectrum = det_tgt_spectrum.minus(det_bgd_spectrum)
        det_spectrum.colour = 'red'
        det_spectrum.label = det_spectrum.units
        ideal_tgt_spectrum.colour = 'blue'
        ideal_tgt_spectrum.label = 'noiseless observation'
        # plot.spectra([det_spectrum, ideal_tgt_spectrum],
        #              title="Bgd subtracted, V = {:6.3f} km/s".format(v_obs), plot_errors=[False, True])
        obs_list.append(det_spectrum)

    ref_spec = obs_list[0]
    ref_spec.colour = 'green'
    corr_spec_list = []
    for i, spec in enumerate(obs_list):
        colours = ['black', 'orange']
        v_obs = v_obs_list[i]
        title = "Observation v reference, $\Delta$V = {:6.3f} km/s".format(v_obs)
        plot.spectra([spec, ref_spec],
                     plot_points=True, title=title)
        corr_spec, delta_v, corr_max = Correlate.cross_correlate(spec, ref_spec)
        print(delta_v, corr_max)

        corr_spec.colour = 'royalblue'
        title = "Correlation, $\Delta$V = {:6.3f} km/s".format(v_obs)
        plot.spectra([corr_spec], wlim=[-14900, 14900],
                     plot_points=True, title=title)
        corr_spec_list.append(corr_spec)

    colours = ['black', 'orange', 'royalblue']
    plot.spectra(obs_list, colours=colours,
                 plot_points=True, title='All observations')

    print('Done')
