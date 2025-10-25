[![DOI](https://zenodo.org/badge/1083289004.svg)](https://doi.org/10.5281/zenodo.17444636)
User Manual 
BY: Tianxiang Ren (rentianxiang1100@163.com)
Data: 2025-10-25

0. Requirements & Conventions
MATLAB: R2022b or later (tested on R2024a).
Toolboxes: Signal Processing Toolbox (for resample, butter, filtfilt, spectrogram).

1. Program Description for Figures 2 and 3 (free_oscillation_spectrum_show.m)
Purpose.
Produce spectrum figure(fig.2. and fig.3.)—Fourier amplitude spectra of in-plane mean normal stress (σₘ) and maximum in-plane shear (τₘₐₓ, implemented as the magnitude of complex shear). All mode annotations and stylistic embellishments are removed.

Signal definitions.
σₘ(t) = 0.5·(σₓₓ + σᵧᵧ) (real)
complex shear q(t) = 0.5·(σₓₓ − σᵧᵧ) + i·τₓᵧ; τₘₐₓ spectrum uses |FFT(q)|

Segment selection: take the last len_hr hours from the tail (or use idx), with an optional t0_hr offset from the end.

Detrend: remove linear trend for σₘ and q(t).
Resampling rule: enforce Fs_eff >= 2.2 *2* bp2. If Fs_tgt is lower, resample up using resample with rational approximation (rat(...,1e-6)). If Fs_tgt empty → default 0.1 Hz before enforcement.
Bandpass (optional): 4th-order Butterworth in [bp1, bp2] Hz with zero-phase filtfilt, only if Fs_eff > 2*bp2.
Spectrum: Hann window (periodic), normalize by sum(w), single-sided with inner bins ×2.
Frequency axis: create a common axis between fmin_mhz and fmax_mhz with step = median Δf of the raw FFT (fallback Fs_eff/N * 1000). Interpolate both spectra to this axis (linear).

Minimal plotting.
Colors fixed: σₘ in red [0.85 0.33 0.10], τₘₐₓ in blue [0 0.45 0.74].
Axes: Frequency (mHz) vs Amplitude (Pa).
Only axis limits, grid, minor ticks, title, and legend. No mode labels, no bottom subpanels.

Inputs.
len_hr, testx, Fs, Fs_tgt, bp1, bp2. testx(:,4:6) = [σₓₓ σᵧᵧ τₓᵧ].
Name–value options: fmin_mhz, fmax_mhz, YLim, Title, t0_hr, idx.

Outputs.
R.f_mHz, R.A_sigma_m, R.A_tau_max for reproducibility and downstream checks.
R.Fs_eff, R.bp_Hz, R.dur_hr metadata, and figure/axes handles.

What's intentionally omitted.
Any mode labeling, binning, dead-zone logic, multi-panel layouts, fine typography, and style knobs.

Visualization-oriented filter choice.
To prevent 0S2 (~0.309 mHz) from dominating the y-axis and masking higher-frequency toroidal bands, the broadband spectrum uses a 4th-order Butterworth band-pass whose lower corner is set near 0.309 mHz. This intentionally compresses the 0S2 peak (≈ −3 dB at the corner) while leaving the band above ~0.6 mHz effectively unaffected (attenuation ≤ 0.02 dB at ≥ 2×fc for a 4th-order Butterworth). All quantitative estimates (e.g., amplitudes/Q for modes) are obtained from dedicated windows/filters placing the mode well inside the passband; this visualization setting is for display only and is documented here for clarity.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example
% Inputs (illustrative)
len_hr = 88;
Fs     = 10;                  % original sampling rate (Hz)
Fs_tgt = 0.05;                   % let the function default/enforce rule
bp1 =0.309e-3;  %         
bp2 =10e-3 ;             
% testx must be N x 6; columns 4:6 = [sigma_xx sigma_yy tau_xy] in Pa
% e.g., testx = [zeros(N,3), sigma_xx, sigma_yy, tau_xy];

R = free_oscillation_spectrum_show( ...
    len_hr, testx, Fs, Fs_tgt, bp1, bp2, ...
    'fmin_mhz',0.2, 'fmax_mhz',5.5, 'YLim',[0 1.3], ...
    'Title','Fourier Amplitude Spectrum (0.2–5.5 mHz)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

2. Program Description for Table1 
Purpose & Scope
“free_oscillation_spectrum_modalleakage.m” supports modal leakage cross-checks. It runs a linear chain (segment selection / resampling guard / 4th-order bandpass / Hann window, single-sided FFT) and performs peak picking on native FFT bins, returning only the peak frequency lists (mHz) for σm and a proxy of τmax. The parameter set below is designed for a stable, review-ready baseline consistent with prior textual analyses.

Rationale:
Full-record idx maximizes duration and frequency resolution df, stabilizing “resolvable vs. unresolvable” decisions. Global-fraction prominence (2%) keeps a consistent, energy-referenced threshold across datasets, emphasizing the most representative spectral peaks. MinPeakDistance = 0.05 mHz effectively suppresses duplicate detections at the current duration/band, focusing the table on auditable primary peaks. Native-bin picking avoids interpolation bias, preserving window response and df for traceable attribution.

Output:
R.peak_freq_mHz: {[1×K1 double], [1×K2 double]} peak frequencies for σm and τmax，with K≤40.
THEN use "modeleak_sigmam.m" and "modeleak_taumax.m" to analysis the modal leakage.

Both modeleak_sigmam.m and modeleak_taumax.m return a table T_out listing up to 40 observed peak frequencies (for σm or τmax, respectively) matched to the nearest PREM mode. Columns are: idx, observed frequency (µHz), suspected leakage source mode (leakage_source_mode), matched PREM frequency (prem_uHz), and signed deviation (*_delta_uHz); the same information is also printed to the console.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 
len_hr = 88;
Fs     = 10;                  % original sampling rate (Hz)
Fs_tgt = 0.05;                   % let the function default/enforce rule
bp1 =0.28e-3;  %         
bp2 =10e-3 ;   

r = free_oscillation_spectrum_modalleakage( ...
    len_hr, testx, Fs, Fs_tgt, bp1, bp2, ...
    'fmin_mhz', 0.2, ...
    'fmax_mhz', 5.5, ...
    'ProminenceMode','globalfrac', ...   % MinPeakProminence = 2% of band max
    'PromFrac', 0.02, ...
    'MinPeakDistance_mHz', 0.05,  ...   % suppress duplicates / pseudo doublets
    'idx', (1:size(testx,1)).' );

modeleak_sigmam(r, 'modes_catalog_all.csv', 1.5, 3);
modeleak_taumax(r, 'modes_catalog_all.csv', 1.5, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

3. Program Description for Figures 4 
"tfr_sharp_horizontal.m"
Purpose
Produce a clean time–frequency map (STFT) for ultra-low frequencies (e.g., 0.2–11 mHz), with anti-alias resampling, zero-phase band-pass filtering, robust per-frequency background removal, and a vertical colorbar outside on the right. No mode labels or decorative elements.

Inputs
x — Real (or complex) time series vector (row or column).
Fs — Original sampling rate in Hz.

Name–Value Options (defaults in brackets)
'Mode' — Analysis mode placeholder (STFT is always used) ['stft'].
'WindowType' — STFT window: 'gaussian'|'hanning'|'hann'|'hamming' ['hanning'].
'TwinHours' — Window length in hours for STFT [15].
'HopMin' — Hop size in minutes for STFT [10].
'FsTgt' — Target sampling rate after anti-alias resampling (Hz) [auto: >= 4.4*bp2].
'Alpha' — Gaussian window shape if 'gaussian' is selected [3.5].
'SigmaT', 'SigmaF' — Kept for API compatibility (not used).
'Norm' — Per-frequency normalization: 'median' or 'zscore' ['median'].
'Clip_dB' — dB floor (negative, applied after normalization) [-40 dB].
'bp1', 'bp2' — Band-pass edges in Hz [0.2e-3, 11e-3].

Processing Pipeline (summary)
Detrend and mean-remove; linearly fill NaNs if present.
Anti-alias resample to FsTgt (rational resampling).
4th-order zero-phase Butterworth band-pass [bp1, bp2].
STFT with selected window, TwinHours and HopMin.
Keep only frequencies within the band; convert axes to hours and mHz.
Per-frequency background removal (median or zscore), then convert to dB, normalize to peak 0 dB, and floor at Clip_dB.
Plot with imagesc and vertical colorbar outside (right).

Outputs
out is a struct with:
out.Fs_ds — Sampling rate after resampling (Hz).
out.T_hours — Time axis (hours).
out.F_mHz — Frequency axis (mHz, within [bp1,bp2]).
out.P_dB — Normalized, clipped dB matrix for plotting.
out.params — Echo of key parameters.

Plot Characteristics
X-axis: time in hours; Y-axis: frequency in mHz.
Colormap: parula (changeable).
Colorbar: vertical, at the right, outside the axes ('eastoutside').
Fonts: Times New Roman, size 14.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example channels (e.g., mean in-plane normal stress)
sigmam = 0.5*(testx(:,4) + testx(:,5));

out = tfr_sharp_horizontal4_min( ...
    sigmam, 10, ...
    'Mode','stft', ...
    'WindowType','hanning', ...
    'TwinHours',15, ...
    'HopMin',10, ...
    'FsTgt',0.1, ...
    'Alpha',3.5, ...          % only used if WindowType='gaussian'
    'SigmaT',6, 'SigmaF',0.6, ... % recorded, not used
    'Norm','median', ...
    'Clip_dB',-40, ...
    'bp1',0.2e-3, ...
    'bp2',11e-3 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Practical Notes
FsTgt rule: if not provided, it’s set to >= 2.2 × (2×bp2) to protect against aliasing.
Window choice: 'hanning' (or 'hann') is a robust default for spectral leakage control at these bands.
Median vs. Z-score: median is more robust to transients and outliers; zscore can highlight weaker, persistent bands.
Units: Passband edges are in Hz; axes are reported in mHz and hours.
Toolboxes: Requires Signal Processing Toolbox (resample, butter, filtfilt, spectrogram).

Troubleshooting
Blank/very dark plot: relax Clip_dB (e.g., -50) or verify the passband overlaps data.
Aliasing artifacts: increase FsTgt or narrow [bp1, bp2].
Coarse time resolution: decrease TwinHours or HopMin (mind the time–frequency trade-off).

4. Program Description for Figures 5 7 9 11
"bandpass_fft_cosine3.m" performs a zero-phase, frequency-domain bandpass (with cosine tapers), then derives the Hilbert envelope and a smoothed “slow” component. It optionally reconstructs a de-beaten time series and fits an exponential decay to estimate modal Q and e-folding time. A compact figure can be produced for publication review.

Method (summary)
Mirror padding → FFT → centered cosine-taper bandpass → iFFT (zero-phase). Trim edges affected by padding. Hilbert transform → envelope; moving-mean smoothing gives the slow component. Optional de-beating: replace amplitude by long-window mean while preserving instantaneous phase. Log-envelope linear fit on a chosen time window to estimate Q via Q=−πf0/slope, with f0 in Hz and time in seconds.

The input and output parameters can be found in the program comments.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1S2 parameters input
f0 = 0.679850;  % mHz
fband_1S2 = [f0-0.010, f0+0.010];  % [0.66485 0.69485] mHz
Fs = 10;
% fband_0S2 = [0.295 0.325];           % mHz
sm_input=((testx(:,4)+testx(:,5))/2);
t0 = datetime(2025,7,30,7,24,0,'TimeZone','UTC+8');  % 
x = sm_input(:);   
x = detrend(x,'linear'); 

[tDT, y1S2, env1S2, slow1S2, D1S2] = bandpass_fft_cosine3( ...
    x, 10, fband_1S2,'1S2', ...
    'Taper_uHz',12, 'Pad_hours',36, 'Trim_hours',12, ...
    'Smooth_hours',10, 'ReturnDatetime',true, 't0', t0, ...
    'Debeat',true, 'BeatHours',10, ...
    'DecayFitHours',[11 43], 'f0_mHz', f0, ...
    'PlotLog',true, 'MakeFigure',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


5. Program Description for Figures 6 8 10 12
"stress_petal_show.m" renders a single stress-petal plot at a specified time (or sample index) from three in-plane stress components. It is designed for manuscript figures where one snapshot is needed.

Color & sign convention
By editorial convention in this figure: blue = tensile (positive) stress, red = compressive (negative) stress.

Inputs (required)
sx, sy, sxy — column vectors of stress components (Pa)
Fs — sampling rate (Hz)
t0 — start time as datetime (TimeZone recommended)

What the plot shows
Stress petal (tension in blue, compression in red) in the local in-plane frame. A light dashed line with a small arrow marking the epicentral direction. The principal stress axis that is nearest to the epicentral direction is highlighted (light dashed). 
Three header lines at the top: timestamp, values of s1 and s2 with angle of s1 θ, and the minimal axis-angle δθ to the epicentral azimuth. Axes are labeled in Pa (Times New Roman).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = stress_petal_show(testx(:,4), testx(:,5), testx(:,6), Fs, t0, ...
    'DoBandpass', true, 'Fband_mHz', [0.295 0.325], ...
    'SliceTime', datetime(2025,7,31,0,0,0,'TimeZone','Asia/Shanghai'), ...
    'FigurePath','petal_single.png', ...
    'InnerFontSize', 12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



