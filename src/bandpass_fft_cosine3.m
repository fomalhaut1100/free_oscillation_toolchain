function [t_out, y_bp_or_debeat, env, env_slow, decay] = ...
    bandpass_fft_cosine3(x, Fs, fband_mHz, Mode_1, varargin)
%BANDPASS_FFT_COSINE3
% Zero-phase frequency-domain bandpass with optional de-beating,
% envelope/slow-component extraction, and exponential decay fitting.
% A compact, publication-ready workflow that can also produce a single figure.
%
% INPUTS
%   x            N×1 real vector (Pa)
%   Fs           sampling rate (Hz)
%   fband_mHz    [f1 f2] in mHz (e.g., 0S2: [0.295 0.325])
%
% NAME–VALUE PAIRS (compatible with older versions)
%   'Taper_uHz'       (default 12)   Cosine taper width at each band edge (microhertz)
%   'Pad_hours'       (default 36)   Mirror padding length (hours) before FFT
%   'Trim_hours'      (default 12)   Discard at each end after iFFT (hours)
%   'Smooth_hours'    (default 8)    Envelope moving-mean window (hours) -> slow component
%   'ReturnDatetime'  (default true) Return datetime axis if t0 provided; otherwise hours
%   't0'              ([])           Start time (datetime or datenum) for plotting
%   'Debeat'          (default false)Reconstruct de-beaten series (long-window amplitude, original phase)
%   'BeatHours'       (default 12)   Smoothing window for de-beating (hours), ~ 0.5–1× beat period
%   'DecayFitHours'   (default [12 inf]) Fit window [t_start t_end] in hours (relative to segment start)
%   'f0_mHz'          ([])           Central frequency for decay model (mHz); defaults to band mid
%   'PlotLog'         (default true) Use semilog-y for envelope plot (recommended for decay)
%   'MakeFigure'      (default true) Create summary figure
%
% OUTPUTS
%   t_out            Time axis (datetime if ReturnDatetime & t0 provided; otherwise hours)
%   y_bp_or_debeat   Bandpassed (or de-beaten) time series (Pa)
%   env              Hilbert envelope (Pa)
%   env_slow         Smoothed envelope (Pa)
%   decay            Struct with fields:
%                    .Q, .e_fold_hours, .a_model (fitted envelope), .fit_hours, .f0_mHz
%
% Notes:
% - Bandpass is implemented in the frequency domain with a cosine (raised-cosine) taper.
% - Mirror padding reduces edge effects; trimming removes regions most affected by padding.
% - Decay is estimated via linear regression on log(envelope) within the specified fit window.
% SPDX-License-Identifier: Apache-2.0
% Copyright (c) 2025 Tianxiang Ren
% Author: Tianxiang Ren (rentianxiang1100@163.com)

% ---------- Parse inputs ----------
p = inputParser;
addParameter(p,'Taper_uHz',12, @(v) isscalar(v) && v>=0);
addParameter(p,'Pad_hours',36, @(v) isscalar(v) && v>=0);
addParameter(p,'Trim_hours',12, @(v) isscalar(v) && v>=0);
addParameter(p,'Smooth_hours',8, @(v) isscalar(v) && v>=0);
addParameter(p,'ReturnDatetime',true, @islogical);
addParameter(p,'t0',[], @(v) isdatetime(v) || isnumeric(v) || isempty(v));
addParameter(p,'Debeat',false, @islogical);
addParameter(p,'BeatHours',12, @(v) isscalar(v) && v>0);
addParameter(p,'DecayFitHours',[12 inf], @(v) isvector(v) && numel(v)==2);
addParameter(p,'f0_mHz',[], @(v) isempty(v) || (isscalar(v) && v>0));
addParameter(p,'PlotLog',true, @islogical);
addParameter(p,'MakeFigure',true, @islogical);
parse(p,varargin{:});

TAPER  = p.Results.Taper_uHz;
PADH   = p.Results.Pad_hours;
TRIMH  = p.Results.Trim_hours;
SMH    = p.Results.Smooth_hours;
RETDT  = p.Results.ReturnDatetime;
t0     = p.Results.t0;
DEBEAT = p.Results.Debeat;
BEATH  = p.Results.BeatHours;
FITH   = p.Results.DecayFitHours;
f0_mHz = p.Results.f0_mHz;
PLOG   = p.Results.PlotLog;
MKF    = p.Results.MakeFigure;

if isempty(f0_mHz), f0_mHz = mean(fband_mHz); end

% ---------- Mirror padding ----------
x = x(:);
N = numel(x);
padN = min(round(PADH*3600*Fs), max(0, floor(N/10)));
if padN>0
    xpad = [flipud(x(1:padN)); x; flipud(x(end-padN+1:end))];
else
    xpad = x;
end
Npad = numel(xpad);

% ---------- Frequency-domain cosine-tapered bandpass (unit passband gain) ----------
X  = fftshift(fft(xpad));
X  = X(:);
f  = ((0:Npad-1)' - floor(Npad/2)) * (Fs/Npad);  % Hz (centered)
af = abs(f);

W  = zeros(Npad,1);
f1 = fband_mHz(1)*1e-3;  % Hz
f2 = fband_mHz(2)*1e-3;  % Hz
taper = TAPER*1e-6;      % Hz

% Low-edge cosine ramp
if taper>0
    idx = af>=max(f1-taper,0) & af<f1;
    W(idx) = 0.5*(1 - cos(pi*(af(idx)-(f1-taper))/taper));
end
% Passband
idx = af>=f1 & af<=f2;
W(idx) = 1;
% High-edge cosine ramp
if taper>0
    idx = af>f2 & af<=f2+taper;
    W(idx) = 0.5*(1 + cos(pi*(af(idx)-f2)/taper));
end

y_full = ifft(ifftshift(X.*W), 'symmetric');

% Remove padding & trim ends
if padN>0
    y_full = y_full(padN+1:end-padN);
end
trimS = round(TRIMH*3600*Fs);
if 2*trimS>=numel(y_full), trimS=0; end
y_bp  = y_full(1+trimS:end-trimS);

% ---------- Time axis ----------
t_s = (0:numel(y_bp)-1)'/Fs;         % seconds from segment start
t_h = t_s/3600;                      % hours from segment start
if RETDT
    if isempty(t0)
        t_out = t_h;                 % fall back to hours
        RETDT = false;
        warning('ReturnDatetime=true but no t0 provided; returning hours axis.');
    else
        if isnumeric(t0) && ~isdatetime(t0)
            t0 = datetime(t0,'ConvertFrom','datenum');
        end
        t_out = t0 + seconds(trimS/Fs) + seconds(t_s);
    end
else
    t_out = t_h;
end

% ---------- Envelope and slow component ----------
z   = hilbert(y_bp);
env = abs(z);
if SMH>0
    win = max(1, round(SMH*3600*Fs));
    env_slow = movmean(env, win);
else
    env_slow = env;
end

% ---------- Optional de-beating (keep phase, long-window amplitude) ----------
y_bp_or_debeat = y_bp;
if DEBEAT
    phi  = unwrap(angle(z));
    wB   = max(1, round(BEATH*3600*Fs));
    a_db = movmean(env, wB);
    y_bp_or_debeat = real(a_db .* exp(1j*phi));
end

% ---------- Exponential decay fit on log-envelope ----------
fit_lo = FITH(1);
fit_hi = FITH(2);
if isinf(fit_hi), fit_hi = t_h(end); end

mask = (t_h>=fit_lo) & (t_h<=fit_hi);
aa   = env(mask);
aa(aa<=0) = min(aa(aa>0))*1e-6;      % guard against log(0)
tt   = t_s(mask);

Pfit  = polyfit(tt, log(aa), 1);
slope = Pfit(1);
ln_a0 = Pfit(2);

f0   = f0_mHz*1e-3;                  % Hz
Qhat = -pi*f0 / slope;
a_model = exp(ln_a0) * exp(-(pi*f0).*t_s / Qhat);
tau_h   = Qhat/(pi*f0)/3600;         % e-folding time (hours)

decay = struct('Q',Qhat, ...
               'e_fold_hours',tau_h, ...
               'fit_hours',[fit_lo fit_hi], ...
               'a_model',a_model, ...
               'f0_mHz',f0_mHz);

% ---------- Plot (optional) ----------
if MKF
    figure('Name','Bandpass + Envelope + Decay','Color','w', ...
           'Units','pixels','Position',[200 200 850 550]);
    fontName = 'Times New Roman';

    % Top: bandpassed (or de-beaten) trace
    subplot(2,1,1);
    plot(t_out, y_bp_or_debeat, 'b'); grid on
    ttl = [Mode_1,' bandpassed time series'];
    if DEBEAT, ttl = [ttl,' (de-beaten)']; end
    title(ttl,'FontName',fontName,'FontSize',14);
    ylabel('Pa','FontName',fontName,'FontSize',13);
    xlim([t_out(1) t_out(end)]);
    set(gca,'FontName',fontName,'FontSize',13);
    if RETDT
        datetick('x','dd-mmm HH:MM','keeplimits','keepticks');
    end

    % Bottom: envelope, slow component, decay model, and fit window
    subplot(2,1,2); hold on
    if PLOG
        h1 = semilogy(t_out, env,      'Color',[0.6 0.6 0.6], 'DisplayName','Envelope');
        h2 = semilogy(t_out, env_slow, 'r', 'LineWidth',1.2,   'DisplayName','Slow component');
        h3 = semilogy(t_out, a_model,  '--', 'Color',[1 0.5 0], 'LineWidth',1.2, 'DisplayName','Decay model');
    else
        h1 = plot(t_out, env,      'Color',[0.6 0.6 0.6], 'DisplayName','Envelope');
        h2 = plot(t_out, env_slow, 'r', 'LineWidth',1.2,   'DisplayName','Slow component');
        h3 = plot(t_out, a_model,  '--', 'Color',[1 0.5 0], 'LineWidth',1.2, 'DisplayName','Decay model');
    end

    % Fit-window highlight
    fit_lo_hr = fit_lo;
    fit_hi_hr = min(fit_hi, t_h(end));
    if RETDT
        x1 = t_out(1) + hours(fit_lo_hr);
        x2 = t_out(1) + hours(fit_hi_hr);
    else
        x1 = fit_lo_hr;
        x2 = fit_hi_hr;
    end
    yl = get(gca,'YLim');
    patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.92 0.96 1.0], ...
          'EdgeColor','none','FaceAlpha',0.9, 'DisplayName','Fit window');

    % Bring curves to front; tidy axes
    uistack(h1,'top'); uistack(h2,'top'); uistack(h3,'top');
    grid on; xlim([t_out(1) t_out(end)]);
    ylabel('Pa','FontName',fontName,'FontSize',13);
    title(sprintf('Envelope, Slow & Decay  (Q = %.0f, e-fold = %.0f h)', Qhat, tau_h), ...
          'FontName',fontName,'FontSize',14);
    box on                                % show all four axis spines
    set(gca,'LineWidth',0.8,'Layer','top'); % clearer box; grid below
    legend('Location','best','FontName',fontName,'FontSize',11);
    set(gca,'FontName',fontName,'FontSize',13);

    if RETDT
        datetick('x','dd-mmm HH:MM','keeplimits','keepticks');
    end
end
end
