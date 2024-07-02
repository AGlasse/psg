import math
import numpy as np
import scipy.stats as stats
from scipy.stats import spearmanr
from scipy.signal import correlate
from scipy.fft import fft, ifft
# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LinearRegression
# from sklearn.metrics import mean_squared_error
from scipy.stats import kendalltau
# from sklearn.feature_selection import mutual_info_regression
from scipy.spatial.distance import pdist
from scipy.spatial import procrustes
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
# import dcor
import pingouin as pg
import pandas as pd
from scipy.signal import csd
from sklearn.feature_selection import mutual_info_regression

from psgspectrum import PsgSpectrum as Spectrum


class PsgCorrelate:

    def __init__(self):
        PsgCorrelate.fwhm_sigma = 2 * math.sqrt(2.0 * math.log(2.0))
        return

    @staticmethod
    def calculate_velocity_difference(spectrum1, spectrum2, reference_shift, reference_velocity):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.cross_correlate(spectrum1, spectrum2)
        # distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift

        return abs(velocity_difference)

    @staticmethod
    def pearson_correlation(spec1, spec2, velocity_list):
        vals1 = spec1.vals
        vals2 = spec2.vals
        correlation_coefficient, _ = stats.pearsonr(vals1, vals2)

        max_index = np.argmax(correlation_coefficient)
        velocity_difference = velocity_list[max_index % len(velocity_list)]

        return correlation_coefficient, velocity_difference

    @staticmethod
    def spearman_correlation(spec1, spec2):
        vals1 = spec1.vals
        vals2 = spec2.vals
        correlation_coefficient, _ = spearmanr(vals1, vals2)
        return correlation_coefficient

    @staticmethod
    def cross_correlation(spec1, spec2):
        vals1 = spec1.vals
        vals2 = spec2.vals

        cross_correlation = correlate(vals1, vals2, mode='same') # from scipy
        cross_corr2 = np.correlate(vals1, vals2, mode='full') # from numpy

        velocities = np.arange(-len(vals2) + 1, len(vals1))
        max_index = np.argmax(cross_corr2)
        vel_difference = velocities[max_index]

        correlation_coefficient = cross_corr2[max_index] / (np.linalg.norm(vals1) * np.linalg.norm(vals2))

        peak_index = np.argmax(cross_correlation)
        delta_v = peak_index - len(vals1) // 2

        return correlation_coefficient, vel_difference

    @staticmethod
    def fit_peak_curve(x, y):

        def gaussian(x, amplitude, mean, sigma):
            return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

        # Initial guess for the Gaussian parameters
        initial_params = [np.max(y), np.mean(x), 1.0]


        # Perform curve fitting
        fitted_params, _ = curve_fit(gaussian, x, y, p0=initial_params, maxfev=50000)

        # Extract the center of the fitted curve
        center = fitted_params[1]

        return center

    @staticmethod
    def calculate_x_distance(spectrum1, spectrum2):
        peaks1 = PsgCorrelate.identify_peaks(spectrum1.vals)  # Identify peaks in spectrum 1
        peaks2 = PsgCorrelate.identify_peaks(spectrum2.vals)  # Identify peaks in spectrum 2

        x_distance = []

        for peak1 in peaks1:
            for peak2 in peaks2:
                center1 = PsgCorrelate.fit_peak_curve(spectrum1.waves, spectrum1.vals[peak1])
                center2 = PsgCorrelate.fit_peak_curve(spectrum2.waves, spectrum2.vals[peak2])
                distance = center1 - center2
                x_distance.append(distance)

        return x_distance

    @staticmethod
    def calculate_velocity_difference(spectrum1, spectrum2, reference_shift=0.01, reference_velocity=990):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.cross_correlate(spectrum1, spectrum2)
        #distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift

        return abs(velocity_difference)

    @staticmethod
    def calculate_velocity_difference2(spectrum1, spectrum2, reference_shift=0.01, reference_velocity=990):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.calculate_phase_difference_phase_correlation(spectrum1, spectrum2)
        # distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift
        return abs(velocity_difference)

    @staticmethod
    def calculate_velocity_difference3(spectrum1, spectrum2, reference_shift=0.01, reference_velocity=990):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.calculate_phase_difference_fourier(spectrum1, spectrum2)
        # distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift

        return abs(velocity_difference)

    @staticmethod
    def normalize_spectrum(spectrum):
        # Normalize the spectrum to have zero mean and unit variance
        normalized_spectrum = (spectrum.vals - np.mean(spectrum.vals)) / np.std(spectrum.vals)
        return normalized_spectrum # which will also be vals

    @staticmethod
    def cross_correlate(spectrum1, spectrum2):
        """ """
        # Perform cross-correlation
        cross_corr = np.correlate(spectrum1.vals, spectrum2.vals, mode='full') # with all spectra
        # Find the lag with the highest correlation
        lag_with_highest_corr = np.argmax(cross_corr) - (len(spectrum1.vals) - 1)
        n_pts, = cross_corr.shape
        c = 2.997E8
        dv = spectrum1.dv_c * c
        dv_hw = 0.5 * n_pts * dv
        delta_vs = np.arange(0, n_pts) * dv - dv_hw
        corr_spectrum = Spectrum(delta_vs, cross_corr, 0.*cross_corr, spectrum2.label, 'arb. units', 'green')
        return corr_spectrum

    @staticmethod
    def find_corr_max(corr_spectrum):
        delta_vs, cross_corr, _, _, _, _ = corr_spectrum.get()
        max_corr = np.amax(cross_corr)
        return None     # dv, dv_err

    @staticmethod
    def calculate_phase_difference_phase_correlation(spectrum1, spectrum2):
        nperseg = len(spectrum1.vals)  # Number of psg_in points per segment for cross-spectral density
        fs = 1.0 / (spectrum1.waves[1] - spectrum1.waves[0])

        # Compute the cross-spectral density (CSD)
        _, csd_spectrum = csd(spectrum1.vals, spectrum2.vals, fs=fs, nperseg=nperseg, scaling='spectrum')

        # Find the phase difference in radians using the phase of the complex-valued CSD
        phase_difference_radians = np.angle(csd_spectrum)

        # Convert phase difference to microns
        wavelength_difference = 1 / 0.0000135
        phase_difference_microns = phase_difference_radians * wavelength_difference / (2 * np.pi)

        # Calculate the average phase difference
        average_phase_difference_microns = np.mean(phase_difference_microns)

        return abs(average_phase_difference_microns)

    @staticmethod
    def calculate_phase_difference_fourier(spectrum1, spectrum2):
        # Compute the Fourier transforms
        fft_spectrum1 = fft(spectrum1.vals)
        fft_spectrum2 = fft(spectrum2.vals)

        # Find the phase difference in radians using the phase of the complex-valued spectra
        phase_difference_radians = np.angle(fft_spectrum2) - np.angle(fft_spectrum1)

        # Convert phase difference to microns
        wavelength_difference = 1 / 0.0000135 # should this be the diff between two points on x axis
        phase_difference_microns = phase_difference_radians * wavelength_difference / (2 * np.pi)

        return abs(phase_difference_microns)
import numpy as np
import scipy.stats as stats
from scipy.stats import spearmanr
from scipy.signal import correlate
from scipy.fft import fft, ifft
# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LinearRegression
# from sklearn.metrics import mean_squared_error
from scipy.stats import kendalltau
# from sklearn.feature_selection import mutual_info_regression
from scipy.spatial.distance import pdist
from scipy.spatial import procrustes
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
# import dcor
import pingouin as pg
import pandas as pd
from scipy.signal import csd
from sklearn.feature_selection import mutual_info_regression

from psgspectrum import PsgSpectrum as Spectrum


class PsgCorrelate:

    def __init__(self):
        return

    @staticmethod
    def calculate_velocity_difference(spectrum1, spectrum2, reference_shift, reference_velocity):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.cross_correlate(spectrum1, spectrum2)
        #distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift

        return abs(velocity_difference)

    @staticmethod
    def pearson_correlation(spec1, spec2, velocity_list):
        vals1 = spec1.vals
        vals2 = spec2.vals
        correlation_coefficient, _ = stats.pearsonr(vals1, vals2)

        max_index = np.argmax(correlation_coefficient)
        velocity_difference = velocity_list[max_index % len(velocity_list)]

        return correlation_coefficient, velocity_difference

    @staticmethod
    def spearman_correlation(spec1, spec2):
        vals1 = spec1.vals
        vals2 = spec2.vals
        correlation_coefficient, _ = spearmanr(vals1, vals2)
        return correlation_coefficient

    @staticmethod
    def cross_correlation(spec1, spec2):
        vals1 = spec1.vals
        vals2 = spec2.vals

        cross_correlation = correlate(vals1, vals2, mode='same') # from scipy
        cross_corr2 = np.correlate(vals1, vals2, mode='full') # from numpy

        velocities = np.arange(-len(vals2) + 1, len(vals1))
        max_index = np.argmax(cross_corr2)
        vel_difference = velocities[max_index]

        correlation_coefficient = cross_corr2[max_index] / (np.linalg.norm(vals1) * np.linalg.norm(vals2))

        peak_index = np.argmax(cross_correlation)
        delta_v = peak_index - len(vals1) // 2

        return correlation_coefficient, vel_difference

    @staticmethod
    def fourier_transforms_correlation(spectrum1, spectrum2):
        # Step 1: Compute the Fourier transforms
        fft_spectrum1 = fft(spectrum1.vals)
        fft_spectrum2 = fft(spectrum2.vals)

        # Step 2: Compute the complex conjugate of one of the spectra
        conjugate_spectrum2 = np.conj(fft_spectrum2)

        # Step 3: Multiply the two complex-valued spectra
        cross_correlation = fft_spectrum1 * conjugate_spectrum2

        # Step 4: Perform inverse Fourier transform on the multiplied spectrum
        cross_correlation_result = ifft(cross_correlation)

        # Step 5: Find the wavelength corresponding to the peak of the cross-correlation function
        peak_wavelength_index = np.argmax(np.abs(cross_correlation_result))
        peak_wavelength = spectrum1.waves[peak_wavelength_index]

        # Step 6: Compute the phase difference in microns
        wavelength_difference = peak_wavelength - spectrum2.waves[0]  # Assuming both spectra have the same x-values
        phase_difference_microns = wavelength_difference

        print(phase_difference_microns)

    @staticmethod
    def ml_linear_regression_corr(spec1, spec2):
        X = spec1.vals
        Y = spec2.vals

        #Perform linear regression
        #coefficients = np.linalg.inv(X.T @ X) @ X.T @ y

        # Combine the spectra into a feature matrix
        Z = np.vstack((X, Y)).T

        # Calculate the correlation coefficient between the spectra
        correlation_coefficient = np.corrcoef(Z[:, 0], Z[:, 1])[0, 1]

        return correlation_coefficient

    @staticmethod
    def kendalltau_correlation(spec1, spec2):
        tau, _ = kendalltau(spec1.vals, spec2.vals)

        # Tau = the correlation coefficient;

        return tau

    @staticmethod
    def mutual_information_corr(spec1, spec2):
        vals1 = spec1.vals.reshape(-1,1)
        vals2 = spec2.vals

        mi = mutual_info_regression(vals1, vals2)

        return mi

    @staticmethod
    def mantel_test_corr(spec1, spec2):
        vals1 = spec1.vals
        vals2 = spec2.vals
        dissimilarity1 = pdist(vals1.reshape(-1,1), metric='euclidean')
        dissimilarity2 = pdist(vals2.reshape(-1,1), metric='euclidean')

        correlation, p_value = spearmanr(dissimilarity1, dissimilarity2)
        # mantel test correlation
        return correlation, p_value

    @staticmethod
    def procrustes_analysis_corr(spec1, spec2):
        vals1 = np.unique(spec1.vals, axis=0).reshape(-1,2)
        vals2 = np.unique(spec2.vals, axis=0).reshape(-1,2)

        mtx1, mtx2, disparity = procrustes(vals1, vals2)
        # returns procrustes disparity
        return disparity

    @staticmethod
    def distance_correlation(spec1, spec2):
        vals1 = spec1.vals
        vals2 = spec2.vals

        # correlation = dcor.distance_correlation(vals1, vals2)

        return None     # correlation

    @staticmethod
    def partial_correlation(spec1, spec2):
        vals1 = spec1.vals
        vals2 = spec2.vals

        df = pd.DataFrame({'Spectrum1':vals1, 'Spectrum2':vals2})

        partial_corr = pg.partial_corr(data=df, x='Spectrum1', y='Spectrum2')

        partial_corr_coeff = partial_corr['r'][0]
        p_value = partial_corr['p-val'][0]

        return partial_corr_coeff, p_value

    @staticmethod
    def identify_peaks(spectrum):
        y = spectrum  # Spectrum psg_in (y-values)

        # Find peaks using the find_peaks function
        peaks, _ = find_peaks(y)

        return peaks

    @staticmethod
    def fit_peak_curve(x, y):

        def gaussian(x, amplitude, mean, sigma):
            return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)

        # Initial guess for the Gaussian parameters
        initial_params = [np.max(y), np.mean(x), 1.0]


        # Perform curve fitting
        fitted_params, _ = curve_fit(gaussian, x, y, p0=initial_params, maxfev=50000)

        # Extract the center of the fitted curve
        center = fitted_params[1]

        return center

    @staticmethod
    def calculate_x_distance(spectrum1, spectrum2):
        peaks1 = PsgCorrelate.identify_peaks(spectrum1.vals)  # Identify peaks in spectrum 1
        peaks2 = PsgCorrelate.identify_peaks(spectrum2.vals)  # Identify peaks in spectrum 2

        x_distance = []

        for peak1 in peaks1:
            for peak2 in peaks2:
                center1 = PsgCorrelate.fit_peak_curve(spectrum1.waves, spectrum1.vals[peak1])
                center2 = PsgCorrelate.fit_peak_curve(spectrum2.waves, spectrum2.vals[peak2])
                distance = center1 - center2
                x_distance.append(distance)

        return x_distance

    @staticmethod
    def calculate_velocity_difference(spectrum1, spectrum2, reference_shift=0.01, reference_velocity=990):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.cross_correlate(spectrum1, spectrum2)
        #distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift

        return abs(velocity_difference)

    @staticmethod
    def calculate_velocity_difference2(spectrum1, spectrum2, reference_shift=0.01, reference_velocity=990):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.calculate_phase_difference_phase_correlation(spectrum1, spectrum2)
        # distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift
        return abs(velocity_difference)

    @staticmethod
    def calculate_velocity_difference3(spectrum1, spectrum2, reference_shift=0.01, reference_velocity=990):
        # Calculate the x-axis shift (distance) between the two spectra
        distance = PsgCorrelate.calculate_phase_difference_fourier(spectrum1, spectrum2)
        # distance_km = distance * 1e-11 # conversion between um to km
        # Calculate the velocity difference based on the x-axis shift and reference values
        velocity_difference = distance * reference_velocity / reference_shift

        return abs(velocity_difference)

    @staticmethod
    def normalize_spectrum(spectrum):
        # Normalize the spectrum to have zero mean and unit variance
        normalized_spectrum = (spectrum.vals - np.mean(spectrum.vals)) / np.std(spectrum.vals)
        return normalized_spectrum # which will also be vals

    @staticmethod
    def cross_correlate(spectrum1, spectrum2):
        """ """
        # Perform cross-correlation
        cross_corr = np.correlate(spectrum1.vals, spectrum2.vals, mode='full') # with all spectra
        # Find the lag with the highest correlation
        idx_max = np.argmax(cross_corr)
        lag_with_highest_corr = idx_max - (len(spectrum1.vals) - 1)
        n_pts, = cross_corr.shape
        c = 2.997E8
        # dv = spectrum1.dv_c * c
        # dv_hw = 0.5 * n_pts * dv
        # delta_vs = np.arange(0, n_pts) * dv - dv_hw
        label = "Lag = {:d}".format(lag_with_highest_corr)
        n_waves, = spectrum1.vals.shape
        xcorr = np.arange(-n_waves, n_waves-1)
        corr_spectrum = Spectrum(xcorr, cross_corr, 0.*cross_corr, label, 'arb. units', 'green')
        vals = corr_spectrum.vals
        corr_max = np.amax(vals)
        idx_max = np.argmax(vals)
        # Find the sub-sample position of maximum correlation by fitting a Gaussian to the central peak
        peak = np.array(vals[idx_max-19:idx_max+20])
        x = np.arange(-19, +20, 1)
        xpeak = PsgCorrelate.fit_peak_curve(x, peak)
        return corr_spectrum, xpeak, corr_max

    @staticmethod
    def find_corr_max(corr_spectrum):
        delta_vs, cross_corr, _, _, _, _ = corr_spectrum.get()
        max_corr = np.amax(cross_corr)
        return None     # dv, dv_err

    @staticmethod
    def calculate_phase_difference_phase_correlation(spectrum1, spectrum2):
        nperseg = len(spectrum1.vals)  # Number of psg_in points per segment for cross-spectral density
        fs = 1.0 / (spectrum1.waves[1] - spectrum1.waves[0])

        # Compute the cross-spectral density (CSD)
        _, csd_spectrum = csd(spectrum1.vals, spectrum2.vals, fs=fs, nperseg=nperseg, scaling='spectrum')

        # Find the phase difference in radians using the phase of the complex-valued CSD
        phase_difference_radians = np.angle(csd_spectrum)

        # Convert phase difference to microns
        wavelength_difference = 1 / 0.0000135
        phase_difference_microns = phase_difference_radians * wavelength_difference / (2 * np.pi)

        # Calculate the average phase difference
        average_phase_difference_microns = np.mean(phase_difference_microns)

        return abs(average_phase_difference_microns)

    @staticmethod
    def calculate_phase_difference_fourier(spectrum1, spectrum2):
        # Compute the Fourier transforms
        fft_spectrum1 = fft(spectrum1.vals)
        fft_spectrum2 = fft(spectrum2.vals)

        # Find the phase difference in radians using the phase of the complex-valued spectra
        phase_difference_radians = np.angle(fft_spectrum2) - np.angle(fft_spectrum1)

        # Convert phase difference to microns
        wavelength_difference = 1 / 0.0000135 # should this be the diff between two points on x axis
        phase_difference_microns = phase_difference_radians * wavelength_difference / (2 * np.pi)

        return abs(phase_difference_microns)
