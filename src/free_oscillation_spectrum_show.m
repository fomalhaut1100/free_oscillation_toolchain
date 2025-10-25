function R = free_oscillation_spectrum_show(len_hr, testx, Fs, Fs_tgt, bp1, bp2, varargin)
% FREE_OSCILLATION_SPECTRUM_SHOW
% Processing: select tail segment -> detrend -> optional resample to ensure
%             Fs_tgt >= 2.2*2*bp2 -> 4th-order bandpass -> Hann-window FFT
%             -> single-sided spectrum (inner bins x2) -> interpolate to a
%             common frequency axis -> plot (red/blue), axis/legend only.
%
% INPUTS
%   len_hr : length (hours) from tail to include (ignored if idx is given)
%   testx  : N x 6 array; columns 4:6 = [sigma_xx, sigma_yy, tau_xy]
%   Fs     : original sampling rate (Hz)
%   Fs_tgt : target Fs (Hz). If empty, default 0.1; enforced to >= 2.2*2*bp2
%   bp1,bp2: bandpass edges in Hz (4th-order). If empty, skip filtering.
%
% NAME-VALUE OPTIONS (minimal)
%  'fmin_mhz' (0.2)   : lower frequency limit for plotting/interp (mHz)
%  'fmax_mhz' (5.5)   : upper frequency limit for plotting/interp (mHz)
%  'YLim'     ([0 1]) : y-axis limits for amplitude (Pa)
%  'Title'    ('')    : optional title string
%  't0_hr'    (0)     : time offset (hours) from the end to define tail end
%  'idx'      ([])    : explicit index vector (overrides len_hr/t0_hr)
%
% OUTPUT struct R (essential only)
%  R.f_mHz      : common frequency axis (mHz)
%  R.A_sigma_m  : amplitude of sigma_m
%  R.A_tau_max  : amplitude of tau_max (|complex shear|)
%  R.Fs_eff     : effective sampling rate after resampling
%  R.bp_Hz      : [bp1, bp2]
%  R.dur_hr     : duration (hours) used
%  R.fig, R.ax  : figure and axis handles (top panel)
%
% Dependencies: Signal Processing Toolbox (butter, filtfilt, resample).
% Font: Times New Roman (system-available).
%
% SPDX-License-Identifier: Apache-2.0
% Copyright (c) 2025 Tianxiang Ren
% Author: Tianxiang Ren (rentianxiang1100@163.com)

% -------- parse minimal options --------
ip = inputParser; ip.KeepUnmatched = false; ip.PartialMatching = false; ip.CaseSensitive = false;
addParameter(ip,'fmin_mhz',0.2,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'fmax_mhz',5.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'YLim',[0 1],@(v)isnumeric(v)&&numel(v)==2);
addParameter(ip,'Title','',@(s)ischar(s)||isstring(s));
addParameter(ip,'t0_hr',0,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'idx',[],@(v)isnumeric(v)&&isvector(v));
parse(ip,varargin{:}); P = ip.Results;
fmin = P.fmin_mhz; fmax = P.fmax_mhz;

% -------- select tail segment (as in v22/v10) --------
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
idx = idx(idx>=1 & idx<=N); dur_hr = numel(idx)/Fs/3600;

% -------- build two time series: sigma_m (real), complex shear (complex) --------
sig_xx  = double(testx(idx,4));
sig_yy  = double(testx(idx,5));
tau_xy  = double(testx(idx,6));
m_plane = detrend(0.5*(sig_xx + sig_yy),'linear');          % sigma_m
q_cplx  = detrend(0.5*(sig_xx - sig_yy) + 1i*tau_xy,'linear'); % complex shear

% -------- resample to enforce Fs_tgt >= 2.2*2*bp2 (v10 rule) --------
if isempty(Fs_tgt), Fs_tgt = 0.1; end
Fs_tgt_eff = max(Fs_tgt, 2.2*2*bp2);
Fs_eff = Fs;
if Fs_tgt_eff < Fs*0.999
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
[f_mHz_all, A_sig] = v10_fft_amp(m_plane, Fs_eff);
[~,          A_tau] = v10_fft_amp(q_cplx , Fs_eff); % magnitude for complex series

% -------- crop & interpolate to common frequency axis (v10 style) --------
df_mHz = median(diff(f_mHz_all));
if ~isfinite(df_mHz) || df_mHz<=0, df_mHz = (Fs_eff/numel(m_plane))*1000; end
f_common = (fmin:df_mHz:fmax).';
mask = f_mHz_all>=fmin & f_mHz_all<=fmax;
if nnz(mask) >= 2
    y1 = interp1(f_mHz_all(mask), A_sig(mask), f_common, 'linear', nan);
    y2 = interp1(f_mHz_all(mask), A_tau(mask), f_common, 'linear', nan);
else
    f_common = (fmin:0.01:fmax).'; y1 = nan(size(f_common)); y2 = nan(size(f_common));
end

% -------- minimal plot: top panel only (red/blue, axes, legend) --------
fig = figure('Color','w','Units','centimeters','Position',[2 2 15 8.6], 'Renderer','painters'); % 15 cm x (7.8+margin)
ax  = axes(fig); hold(ax,'on');
col1 = [0.85 0.33 0.10]; % red-ish (same as v22)
col2 = [0.00 0.45 0.74]; % blue-ish (same as v22)
plot(ax, f_common, y1, 'LineWidth',0.7, 'Color',col1, 'DisplayName','\sigma_m');
plot(ax, f_common, y2, 'LineWidth',0.7, 'Color',col2, 'DisplayName','\tau_{max}');
xlim(ax,[fmin fmax]); ylim(ax,P.YLim); grid(ax,'on'); ax.XMinorTick='on'; ax.YMinorTick='on';
set(ax,'FontName','Times New Roman','FontSize',10.5,'LineWidth',0.6,'TickLength',[0.012 0.012],'Layer','top','Box','on');
xlabel(ax,'Frequency (mHz)','FontName','Times New Roman','FontSize',11);
ylabel(ax,'Amplitude (Pa)','FontName','Times New Roman','FontSize',11);
ttl = char(P.Title);
if strlength(string(ttl))==0
    ttl = sprintf('Fourier Amplitude Spectrum (%.1f–%.1f mHz) \\approx %d hours', fmin, fmax, round(dur_hr));
end
title(ax, ttl, 'FontName','Times New Roman','FontSize',12.5,'FontWeight','bold');
legend(ax,'Location','northeast','AutoUpdate','off','Interpreter','tex');

% -------- outputs --------
R = struct('f_mHz',f_common(:),'A_sigma_m',y1(:),'A_tau_max',y2(:), ...
           'Fs_eff',Fs_eff,'bp_Hz',[bp1,bp2],'dur_hr',dur_hr, ...
           'fig',fig,'ax',ax);
end

% ---------- helpers (v10 style) ----------
function [f_mHz, A] = v10_fft_amp(x, Fs)
    x = x(:); x = x - mean(x,'omitnan'); x = fillmissing(x,'linear');
    N = numel(x); w = hann(N,'periodic');
    X = fft(x.*w); Xmag = abs(X(1:floor(N/2)+1));
    A = Xmag / sum(w); if N>2, A(2:end-1) = 2*A(2:end-1); end
    fHz = (0:floor(N/2))' * (Fs/N); f_mHz = fHz * 1000;
end

function y = safe_filtfilt(b,a,x)
    x = x(:); nco = max(length(a),length(b)); nfact = 3*(nco-1);
    if length(x) <= nfact, y = filter(b,a,x); return; end
    y = filtfilt(b,a,x);
end
