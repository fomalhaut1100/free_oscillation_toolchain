function out = tfr_sharp_horizontal4_min(x, Fs, varargin)
% Pipeline: detrend -> anti-alias resample -> zero-phase bandpass -> STFT
% -> per-frequency background removal -> dB normalize & clip -> plot.
%
% Inputs:
%   x  : vector (stress time series)
%   Fs : original sampling rate (Hz)
%  Name–value (compatible with full version):
%   'Mode'       : analysis mode (kept for API; STFT is used)      ['stft']
%   'WindowType' : 'gaussian'|'hanning'|'hann'|'hamming'      ['hanning']
%   'TwinHours'  : STFT window length in hours                       [15]
%   'HopMin'     : hop size in minutes                               [10]
%   'FsTgt'      : target Fs after resampling (Hz)          [>=4.4*bp2]
%   'Alpha'      : Gaussian window shape (if WindowType='gaussian') [3.5]
%   'SigmaT'     : kept for API (unused here)                         [6]
%   'SigmaF'     : kept for API (unused here)                         [0.6]
%   'Norm'       : 'median'|'zscore'                           ['median']
%   'Clip_dB'    : floor in dB (negative)                           [-40]
%   'bp1','bp2'  : band edges in Hz                        [0.2e-3,11e-3]
%
% SPDX-License-Identifier: Apache-2.0
% Copyright (c) 2025 Tianxiang Ren
% Author: Tianxiang Ren (rentianxiang1100@163.com)

% -------- parse --------
p = inputParser;
addParameter(p,'Mode','stft');
addParameter(p,'WindowType','hanning');
addParameter(p,'TwinHours',15);
addParameter(p,'HopMin',10);
addParameter(p,'FsTgt',[]);
addParameter(p,'Alpha',3.5);
addParameter(p,'SigmaT',6);
addParameter(p,'SigmaF',0.6);
addParameter(p,'Norm','median');
addParameter(p,'Clip_dB',-40);
addParameter(p,'bp1',0.2e-3);
addParameter(p,'bp2',11e-3);
parse(p,varargin{:});
Mode       = lower(p.Results.Mode);
WindowType = lower(p.Results.WindowType);
TwinH      = p.Results.TwinHours;
HopMin     = p.Results.HopMin;
FsTgt      = p.Results.FsTgt;
Alpha      = p.Results.Alpha;
SigmaT     = p.Results.SigmaT;
SigmaF     = p.Results.SigmaF; 
NormWay    = lower(p.Results.Norm);
Clip_dB    = p.Results.Clip_dB;
bp1        = p.Results.bp1;
bp2        = p.Results.bp2;

% -------- preprocess --------
x = x(:);
if any(isnan(x)), x = fillmissing(x,'linear'); end
x = detrend(x,'linear'); x = x - mean(x);

% -------- anti-alias resample --------
if isempty(FsTgt), FsTgt = max(Fs, 4.4*bp2); end   % >= 2.2*(2*bp2)
[pP,pQ] = rat(FsTgt/Fs);
x_ds  = resample(x,pP,pQ);
Fs_ds = Fs*pP/pQ;

% -------- zero-phase bandpass --------
[b,a] = butter(4, [bp1 bp2]/(Fs_ds/2), 'bandpass');
x_bp  = filtfilt(b,a,x_ds);

% -------- STFT (Mode kept for API compatibility) --------
Nw   = max(32, round(TwinH*3600*Fs_ds));
Nh   = max(1,  round(HopMin*60*Fs_ds));
nfft = 2^nextpow2(max(Nw, 4096));
switch WindowType
    case 'gaussian', w = gausswin(Nw, Alpha);
    case {'hann','hanning'}, w = hann(Nw,'periodic');
    case 'hamming', w = hamming(Nw,'periodic');
    otherwise,      w = hann(Nw,'periodic');
end
nover = max(0, Nw - Nh);
[S,F,T] = spectrogram(x_bp, w, nover, nfft, Fs_ds, 'yaxis');  P = abs(S).^2;

% -------- keep band & convert axes --------
F_mHz  = F*1000;
keepF  = (F_mHz >= bp1*1000 & F_mHz <= bp2*1000);
F_mHz  = F_mHz(keepF);
P      = P(keepF,:);
Thours = T/3600;

% -------- per-frequency background -> dB & clip --------
switch NormWay
    case 'zscore'
        mu = mean(P,2,'omitnan'); sd = std(P,0,2,'omitnan') + eps;
        Pn = max((P - mu)./sd, 0);
    otherwise % 'median'
        medF = median(P,2,'omitnan');
        Pn = max(P - medF, 0);
end
PdB = 10*log10(Pn + eps);
PdB = PdB - max(PdB(:));
PdB = max(PdB, Clip_dB);

% -------- plot (colorbar at right OUTSIDE) --------
fig = figure('Color','w');
ax  = axes('Parent',fig);
imagesc(ax, Thours, F_mHz, PdB); axis(ax,'xy');
colormap(ax, parula); grid(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',14, ...
       'Units','normalized','PositionConstraint','outerposition');
xlabel(ax,'Time (hours)','FontName','Times New Roman');
ylabel(ax,'Frequency (mHz)','FontName','Times New Roman');
title(ax, sprintf('STFT | %.3f–%.3f mHz | Tw=%.1fh, hop=%dmin', ...
      bp1*1000, bp2*1000, TwinH, HopMin), 'FontWeight','normal');

cb = colorbar(ax,'eastoutside');      % vertical, outside on the right
ylabel(cb,'Power (dB)','FontName','Times New Roman');

% -------- output --------
out = struct('Fs_ds',Fs_ds,'T_hours',Thours,'F_mHz',F_mHz,'P_dB',PdB, ...
    'params', struct('Mode',Mode,'WindowType',WindowType,'TwinHours',TwinH, ...
                     'HopMin',HopMin,'FsTgt',FsTgt,'Alpha',Alpha, ...
                     'SigmaT',SigmaT,'SigmaF',SigmaF,'Norm',NormWay, ...
                     'Clip_dB',Clip_dB,'bp1',bp1,'bp2',bp2));
end
