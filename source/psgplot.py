import matplotlib.pyplot as plt
import numpy as np
import sys


class PsgPlot:

    wave_range = [0., 1.]

    def __init__(self):
        return

    @staticmethod
    def set_plot_area(**kwargs):
        figsize = kwargs.get('figsize', [12, 9])
        xlim = kwargs.get('xlim', None)            # Common limits for all plots
        ylim = kwargs.get('ylim', None)            # Common limits for all plots
        xlabel = kwargs.get('xlabel', '')          # Common axis labels
        ylabel = kwargs.get('ylabel', '')
        ncols = kwargs.get('ncols', 1)             # Number of plot columns
        nrows = kwargs.get('nrows', 1)
        remplots = kwargs.get('remplots', None)
        aspect = kwargs.get('aspect', 'auto')      # 'equal' for aspect = 1.0
        fontsize = kwargs.get('fontsize', 16)
        plt.rcParams.update({'font.size': fontsize})

        sharex = xlim is not None
        sharey = ylim is not None
        fig, ax_list = plt.subplots(nrows, ncols, figsize=figsize,
                                    sharex=sharex, sharey=sharey,
                                    squeeze=False)
        fig.patch.set_facecolor('white')

        for i in range(0, nrows):
            for j in range(0, ncols):
                ax = ax_list[i, j]
                ax.set_aspect(aspect)       # Set equal axes
                if xlim is not None:
                    ax.set_xlim(xlim)
                if ylim is not None:
                    ax.set_ylim(ylim)
                if i == nrows-1 and j == 0:
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel(ylabel)
        if remplots is not None:
            rps = np.atleast_2d(remplots)
            for i in range(0, len(rps)):
                ax_list[rps[i, 0], rps[i, 1]].remove()
        return fig, plt, ax_list

    @staticmethod
    def set_wave_range(wave_range):
        PsgPlot.wave_range = wave_range
        return

    @staticmethod
    def spectra(spectra, **kwargs):

        title = kwargs.get('title', '')
        png_path = kwargs.get('png_path', None)
        colours = kwargs.get('colours', None)
        plot_errors = kwargs.get('plot_errors', [False]*len(spectra))
        plot_points = kwargs.get('plot_points', False)
        ls = kwargs.get('ls', 'none')
        xlabel = 'Wavelength [$\mu$m]'
        label0 = "{:s} [{:s}]".format(spectra[0].label, spectra[0].units)
        ylabel = kwargs.get('ylabel', label0)

        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        ax.set_title(title)
        waves = spectra[0].waves
        wlim = kwargs.get('wlim', PsgPlot.wave_range)
        ax.set_xlim(wlim)
        # Get (Y) plot limits within selected wavelength limits
        ylim = [sys.float_info.max, sys.float_info.min]
        indices = np.argwhere(np.logical_and(waves > wlim[0], waves < wlim[1]))
        idxmin, idxmax = indices[0, 0], indices[-1, 0]
        for spectrum in spectra:
            vals = spectrum.vals[idxmin:idxmax]
            ymin, ymax = np.amin(vals), np.amax(vals)
            ylim[0] = ymin if ymin < ylim[0] else ylim[0]
            ylim[1] = ymax if ymax > ylim[1] else ylim[1]

        ylim[0] = kwargs.get('ymin', ylim[0])
        dylim = ylim[1] - ylim[0]
        ylim = ylim if dylim > .01 * ylim[0] else [ylim[0] - dylim, ylim[1] + dylim]
        ylim = kwargs.get('ylim', ylim)
        ax.set_ylim(ylim)

        for i, spectrum in enumerate(spectra):
            waves, vals, val_errs, units, label, spec_colour = spectrum.get()
            color = spec_colour if colours is None else colours[i]

            ax.set_xlabel(xlabel, fontsize=16.0)
            ax.set_ylabel(ylabel, fontsize=16.0)

            fillstyle = 'full' if plot_points else 'none'
            marker = '|' if plot_points else 'none'
            key_words = {'marker': marker, 'ms': 10., 'mew': 2.0, 'fillstyle': fillstyle,
                         'label': label, 'color': color}
            ax.plot(waves, vals, **key_words)
            if plot_errors[i]:
                if plot_points:
                    ax.ploterrors(waves, vals, yerr=val_errs)
                else:
                    key_words = {'marker': marker, 'ms': 1., 'mew': 2.0, 'fillstyle': fillstyle,
                                 'label': label, 'color': color}
                    ax.fill_between(waves, vals - val_errs, vals + val_errs, alpha=0.5)

        plt.legend()
        if png_path is not None:
            plt.savefig(png_path, bbox_inches='tight')
            plt.close(fig)
        plt.show()
        return
