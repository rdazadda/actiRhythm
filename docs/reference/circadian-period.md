# Endogenous Circadian Period Estimation (Lomb-Scargle)

Estimates the dominant (endogenous) circadian PERIOD (tau) of an
activity time series using the Lomb-Scargle periodogram. Unlike the
classical Fast Fourier Transform (FFT), the Lomb-Scargle method does not
require evenly sampled data, so it correctly handles the irregular and
gappy sampling that results from non-wear periods, dropped epochs, or
mixed epoch lengths in accelerometer recordings.
