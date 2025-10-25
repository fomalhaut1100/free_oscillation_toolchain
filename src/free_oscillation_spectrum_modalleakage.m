function R = free_oscillation_spectrum_modalleakage(len_hr, testx, Fs, Fs_tgt, bp1, bp2, varargin)
% FREE_OSCILLATION_SPECTRUM_MODALLEAKAGE
% Purpose (for modal leakage studies):
%   STRICTLY on the NATIVE FFT FREQUENCY BINS (no interpolation), and
%   return ONLY the peak frequency lists for sigma_m and tau_max.
%
%   Returned field:
%     R.peak_freq_mHz : { [1xK1 double], [1xK2 double] }, K<=Ntop (default 40)
%       - Channel 1: sigma_m = 0.5*(sigma_xx + sigma_yy)  (real)
%       - Channel 2: tau_max proxy via complex shear q(t) = 0.5*(sigma_xx - sigma_yy) + 1i*tau_xy
%         Peaks are found on native FFT bins within [fmin_mhz, fmax_mhz], then:
%           (i) filter by MinPeakProminence and MinPeakDistance,
%           (ii) sort by peak HEIGHT descending, take top Ntop,
%           (iii) sort ascending by frequency for output.
%
% Inputs
%   len_hr : tail length (hours) to include if 'idx' not provided
%   testx  : N x 6, columns 4:6 = [sigma_xx, sigma_yy, tau_xy] (Pa)
%   Fs     : original sampling rate (Hz)
%   Fs_tgt : target sampling rate (Hz); if empty, default 0.1; will be enforced to >= 2.2*2*bp2
%   bp1,bp2: bandpass edges (Hz). If empty/invalid, bandpass is skipped.
%
% Name-Value Options (minimal but leakage-oriented)
%   'fmin_mhz' (0.2)    : lower freq limit for peak search (mHz)
%   'fmax_mhz' (5.5)    : upper freq limit for peak search (mHz)
%   't0_hr'    (0)      : time offset (hours) from the end
%   'idx'      ([])     : explicit index vector (overrides len_hr/t0_hr)
%   'Ntop'     (40)     : keep at most Ntop peaks per channel
%   'MinPeakDistance_mHz' ([]) : if empty, defaults to max(2*df_mHz, 0.005)
%   'ProminenceMode' ('localnoise') : 'localnoise' | 'globalfrac'
%   'PromFrac'  (0.02)  : used when ProminenceMode='globalfrac' (2% of band max)
%   'NoiseK'    (6)     : k for local-noise prominence threshold
%   'NoiseWindowBins' (31): window (bins) for moving-median baseline
%
% Output
%   R.peak_freq_mHz : cell {freqs_sigma_m, freqs_tau_max} in mHz (ascending)
%
% Notes
%   - No plotting, no mode labels, no amplitude arrays are returned—this is
%     a lean extractor for peak FREQUENCIES to feed a leakage matcher.
%   - Peak picking on NATIVE bins preserves df and avoids interpolation bias.
%
% Dependencies: Signal Processing Toolbox (butter, filtfilt, resample, findpeaks).
%
% SPDX-License-Identifier: Apache-2.0
% Copyright (c) 2025 Tianxiang Ren
% Author: Tianxiang Ren (rentianxiang1100@163.com)

% -------- parse options --------
ip = inputParser; ip.KeepUnmatched=false; ip.PartialMatching=false; ip.CaseSensitive=false;
addParameter(ip,'fmin_mhz',0.2,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'fmax_mhz',5.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'t0_hr',0,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'idx',[],@(v)isnumeric(v)&&isvector(v));
addParameter(ip,'Ntop',40,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'MinPeakDistance_mHz',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(ip,'ProminenceMode','localnoise',@(s)ischar(s)||isstring(s));
addParameter(ip,'PromFrac',0.02,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'NoiseK',6,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'NoiseWindowBins',31,@(x)isnumeric(x)&&isscalar(x)&&x>=3&&mod(x,2)==1);
parse(ip,varargin{:}); P = ip.Results;
fmin = P.fmin_mhz; fmax = P.fmax_mhz; Ntop=P.Ntop;

% -------- select tail segment (v22/v10 style) --------
N = size(testx,1);
if ~isempty(P.idx)
    idx = P.idx(:);
else
    need   = round(len_hr*3600*Fs);
    tailEnd= N - round(P.t0_hr*3600*Fs);
    tailEnd= max(min(tailEnd, N), 1);
    idx1   = max(1, tailEnd - need + 1);
    idx    = (idx1:tailEnd).';
end
idx = idx(idx>=1 & idx<=N);

% -------- build two time series --------
sig_xx  = double(testx(idx,4));
sig_yy  = double(testx(idx,5));
tau_xy  = double(testx(idx,6));
m_plane = detrend(0.5*(sig_xx + sig_yy), 'linear');               % sigma_m (real)
q_cplx  = detrend(0.5*(sig_xx - sig_yy) + 1i*tau_xy, 'linear');   % complex shear (proxy for tau_max)

% -------- resample to enforce Fs_tgt >= 2.2*2*bp2 (v10 rule) --------
if isempty(Fs_tgt), Fs_tgt = 0.1; end
Fs_tgt_eff = max(Fs_tgt, 2.2*2*bp2);
Fs_eff = Fs;
if abs(Fs_tgt_eff - Fs) > 1e-9
    [p_ds,q_ds] = rat(Fs_tgt_eff/Fs, 1e-6);
    m_plane = resample(m_plane, p_ds, q_ds);
    q_cplx  = resample(q_cplx , p_ds, q_ds);
    Fs_eff  = Fs * p_ds / q_ds;
end

% -------- optional 4th-order bandpass (if feasible) --------
if ~isempty(bp1) && ~isempty(bp2) && bp1>0 && bp2>bp1 && Fs_eff > 2*bp2
    Wn = [bp1 bp2]/(Fs_eff/2);
    [b,a] = butter(4, Wn, 'bandpass');
    m_plane = safe_filtfilt(b,a,m_plane);
    qr = safe_filtfilt(b,a,real(q_cplx));
    qi = safe_filtfilt(b,a,imag(q_cplx));
    q_cplx = complex(qr, qi);
end

% -------- v10 FFT amplitude (Hann, sum(w), single-sided with inner*2) --------
[f_sig_mHz, A_sig] = v10_fft_amp(m_plane, Fs_eff);
[f_tau_mHz, A_tau] = v10_fft_amp(q_cplx , Fs_eff);

% -------- native-bin peak picking within [fmin,fmax] --------
df_mHz = median(diff(f_sig_mHz));
if ~isfinite(df_mHz) || df_mHz<=0, df_mHz = (Fs_eff/numel(m_plane))*1000; end
if isempty(P.MinPeakDistance_mHz)
    minDist = max(2*df_mHz, 0.005); % >= two bins, and at least 5 μHz
else
    minDist = P.MinPeakDistance_mHz;
end

peakCell = cell(1,2);
peakCell{1} = pick_native_peaks(f_sig_mHz, A_sig, fmin, fmax, Ntop, ...
                                P.ProminenceMode, P.PromFrac, P.NoiseK, P.NoiseWindowBins, minDist);
peakCell{2} = pick_native_peaks(f_tau_mHz, A_tau, fmin, fmax, Ntop, ...
                                P.ProminenceMode, P.PromFrac, P.NoiseK, P.NoiseWindowBins, minDist);

% -------- output (only what leakage analysis needs) --------
R = struct('peak_freq_mHz', {peakCell});

end

% ===================== helpers =====================

function locs_out = pick_native_peaks(f_mHz, A, fmin, fmax, Ntop, mode, promFrac, kNoise, winBins, minDist)
    % Restrict to band
    mask = f_mHz>=fmin & f_mHz<=fmax & isfinite(A) & ~isnan(A);
    locs_out = [];
    if nnz(mask) < 3 || max(A(mask))<=0
        return;
    end
    F = f_mHz(mask); Y = A(mask);

    % Prominence threshold
    switch lower(string(mode))
        case "globalfrac"
            prom = max(promFrac * max(Y), eps);
        otherwise % "localnoise"
            % Robust local baseline via moving median; residual > 0 approximates peaks over baseline
            try
                bg = movmedian(Y, winBins, 'omitnan');
            catch
                % Fallback if movmedian unavailable
                bg = medfilt1(Y, winBins);
            end
            resid = Y - bg;
            pos   = resid(resid>0);
            if isempty(pos)
                nf = median(abs(resid),'omitnan'); % fallback
            else
                nf = median(pos,'omitnan');        % robust "noise floor" proxy
            end
            prom = max(kNoise * nf, promFrac * max(Y)); % guard with global 2% as floor
    end

    % Peak picking on native bins
    [pks, locs] = findpeaks(Y(:), F(:), ...
                    'MinPeakProminence', prom, ...
                    'MinPeakDistance',   minDist);

    if isempty(pks)
        return;
    end

    % Sort by peak HEIGHT (desc), keep top N
    [~, ordDesc] = sort(pks, 'descend');
    keep = min(Ntop, numel(ordDesc));
    locs = locs(ordDesc(1:keep));

    % Return ascending by frequency
    locs_out = sort(locs(:)).';
end

function [f_mHz, A] = v10_fft_amp(x, Fs)
    % v10-style FFT amplitude: Hann window, sum(w) normalization, single-sided with inner bins * 2
    x = x(:); x = x - mean(x,'omitnan'); x = fillmissing(x,'linear');
    N = numel(x);
    w = hann(N,'periodic');
    X = fft(x.*w);
    Xmag = abs(X(1:floor(N/2)+1));
    A = Xmag / sum(w);
    if N>2, A(2:end-1) = 2*A(2:end-1); end
    fHz = (0:floor(N/2))' * (Fs/N);
    f_mHz = fHz * 1000;
end

function y = safe_filtfilt(b,a,x)
    x = x(:);
    nco = max(length(a),length(b));
    nfact = 3*(nco-1);
    if length(x) <= nfact
        % Short data fallback: single-pass filter (phase shift), avoids filtfilt failure.
        y = filter(b,a,x);
        return;
    end
    y = filtfilt(b,a,x);
end
