function out = stress_petal_show(sx, sy, sxy, Fs, t0, varargin)
% SPDX-License-Identifier: Apache-2.0
% Copyright (c) 2025 Tianxiang Ren
% Author: Tianxiang Ren (rentianxiang1100@163.com)
% STRESS_PETAL_SHOW 
%
% ---------- parse ----------
p = inputParser; p.KeepUnmatched=true; p.PartialMatching=false; p.CaseSensitive=false;
p.addParameter('DoBandpass', true);
p.addParameter('Fband_mHz', [0.295 0.325]);
p.addParameter('f0_mHz', []);
p.addParameter('Taper_uHz', 12);
p.addParameter('Pad_hours', 36);
p.addParameter('Trim_hours', 12);
p.addParameter('Smooth_hours', 8);
p.addParameter('Debeat', true);
p.addParameter('BeatHours', 14);
p.addParameter('DecayFitHours', [24 84]);
p.addParameter('MakeFilterFigures', false);

p.addParameter('SliceTime', datetime.empty);
p.addParameter('Index', []);

p.addParameter('EventLat', 52.498);
p.addParameter('EventLon', 160.264);
p.addParameter('EventTime', datetime(2025,7,29,23,24,52,'TimeZone','UTC'));
p.addParameter('StaLat', 40.133);
p.addParameter('StaLon', 113.233);
p.addParameter('XAzimuthDeg', 25);

p.addParameter('DrawEpicentralAz', true);
p.addParameter('EpicArrowReverse', true);
p.addParameter('EpicArrowColor', [0.75 0.75 0.75]);
p.addParameter('EpicArrowHeadLengthFrac', 0.10);
p.addParameter('EpicArrowHeadHalfAngleDeg', 15);

p.addParameter('LineGray', [0.9 0.9 0.9]);
p.addParameter('LineStyleGray', '--');
p.addParameter('LineWidthGray', 1.2);

p.addParameter('DrawS1Az', true);
p.addParameter('StressAxisLineStyle', '--');
p.addParameter('StressAxisLighten', 0.65);
p.addParameter('StressAxisAlpha', 0.45);
p.addParameter('StressAxisLineWidth', 1.2);

p.addParameter('InnerFontSize', []);
p.addParameter('ShowAngleRow', true);

p.addParameter('FigurePath', 'stress_petal_show.png');
p.parse(varargin{:}); a = p.Results;
if isempty(a.f0_mHz), a.f0_mHz = mean(a.Fband_mHz); end

% ---------- data prep ----------
sx = sx(:); sy = sy(:); sxy = sxy(:);
N  = min([numel(sx) numel(sy) numel(sxy)]);
sx = sx(1:N); sy = sy(1:N); sxy = sxy(1:N);

% time axis check  ✅ 修正
if ~isa(t0,'datetime')
    error('t0 must be a valid datetime (with TimeZone recommended).');
end

% bandpass (optional)
if a.DoBandpass
    [tt,  sx_bp]  = local_bp(sx,  Fs, t0, a, a.MakeFilterFigures);
    [tt2, sy_bp]  = local_bp(sy,  Fs, t0, a, false);
    [tt3, sxy_bp] = local_bp(sxy, Fs, t0, a, false);
    if ~isequal(tt,tt2) || ~isequal(tt,tt3)
        sy_bp  = interp1(datenum(tt2), sy_bp,  datenum(tt),  'linear','extrap');
        sxy_bp = interp1(datenum(tt3), sxy_bp, datenum(tt), 'linear','extrap');
    end
else
    tt = t0 + seconds((0:N-1).'/Fs);
    sx_bp = sx; sy_bp = sy; sxy_bp = sxy;
end

% principal stresses (θ relative to station X-axis, CCW+, [0,180))
[s1_all, s2_all, theta_deg_all] = principal_stress(sx_bp, sy_bp, sxy_bp);

% choose slice index
if ~isempty(a.Index)
    j = max(1, min(numel(tt), round(a.Index)));
elseif ~isempty(a.SliceTime)
    [~, j] = min(abs(tt - a.SliceTime));
else
    [~, j] = max(abs(s1_all));  % 默认：|s1| 最大
end
t_sel = tt(j); s1 = s1_all(j); s2_val = s2_all(j); theta_deg = theta_deg_all(j);   % ✅ s2_val

% scale
abs_all = abs([s1_all; s2_all]);
smax = prctile(abs_all(isfinite(abs_all)), 99);
if ~isfinite(smax) || smax<=0, smax = max(abs([s1; s2_val; 1])); end
smax = ceil(smax);
headroom = 1.28;

% epicentral azimuth
az_epic_deg = initial_bearing_deg(a.StaLat, a.StaLon, a.EventLat, a.EventLon);
phi_epic_plot = wrapTo180(90 - az_epic_deg);

% s1 geographic & delta
theta_s1_geo_deg_all = wrapTo360(a.XAzimuthDeg - theta_deg_all(:)');
delta_all = abs(wrapTo180(theta_s1_geo_deg_all - az_epic_deg));
delta_all(delta_all>90) = 180 - delta_all(delta_all>90);
delta_sel = delta_all(j);

% ---------- petal (ellipse) ----------
M = 360; phi = linspace(0,2*pi,M);
cphi = cos(phi); sphi = sin(phi);
c2   = cphi.^2;
s2phi= sphi.^2;                 % ✅ 改名，避免遮蔽 s2 主应力

ct = cosd(theta_deg); st = sind(theta_deg);

% ✅ 用按元素乘；并使用 s2_val
r  = s1.*c2 + s2_val.*s2phi;
x0 = ( r.*cphi)*ct - ( r.*sphi)*st;
y0 = ( r.*cphi)*st + ( r.*sphi)*ct;
posMask = r >= 0; negMask = ~posMask;

% ---------- figure ----------
fig = figure('Visible','on','Color','w');
set(fig,'Units','centimeters','Position',[1 1 4 4*headroom]);
set(fig,'DefaultTextFontName','Times New Roman', ...
        'DefaultAxesFontName','Times New Roman', ...
        'DefaultTextFontSize',11,'DefaultAxesFontSize',12);

ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); pbaspect(ax,[1 1 1]);
xlim(ax,[-smax smax]); ylim(ax,[-smax smax*headroom]);
ylabel(ax,'\fontname{Times New Roman}Stress (Pa)','Interpreter','tex');
xlabel(ax,'\fontname{Times New Roman}Stress (Pa)','Interpreter','tex');

% epicentral line + arrow
if a.DrawEpicentralAz
    ce = cosd(phi_epic_plot); se = sind(phi_epic_plot);
    plot(ax, smax*[-ce ce], smax*[-se se], ...
        'LineStyle',a.LineStyleGray,'LineWidth',a.LineWidthGray,'Color',[0.90 0.90 0.90]);

    D = 1; if a.EpicArrowReverse, D = -1; end
    Larr = 0.92*smax; theta_arrow = phi_epic_plot; if a.EpicArrowReverse, theta_arrow = theta_arrow + 180; end
    HL   = a.EpicArrowHeadLengthFrac*smax; alpha = a.EpicArrowHeadHalfAngleDeg;
    x_tip = D*Larr*ce; y_tip = D*Larr*se;
    xL = x_tip + HL*cosd(theta_arrow + 180 - alpha);
    yL = y_tip + HL*sind(theta_arrow + 180 - alpha);
    xR = x_tip + HL*cosd(theta_arrow + 180 + alpha);
    yR = y_tip + HL*sind(theta_arrow + 180 + alpha);
    plot(ax,[x_tip xL],[y_tip yL],'-','LineWidth',a.LineWidthGray,'Color',a.EpicArrowColor);
    plot(ax,[x_tip xR],[y_tip yR],'-','LineWidth',a.LineWidthGray,'Color',a.EpicArrowColor);
end

% petal curves
colPos = [0 0.447 0.741]; colNeg = [0.85 0.325 0.098];
plot(ax, x0(posMask), y0(posMask), 'LineWidth',1.8, 'Color',colPos);
plot(ax, x0(negMask), y0(negMask), 'LineWidth',1.8, 'Color',colNeg);

% select axis closer to epicentral azimuth
if a.DrawS1Az
    theta1 = wrapTo180(theta_deg);
    theta2 = wrapTo180(theta1 + 90);
    d1 = axis_angle_diff(theta1, phi_epic_plot);
    d2 = axis_angle_diff(theta2, phi_epic_plot);
    if d1 <= d2
        ang_sel = theta1; val_sel = s1; which_s = 1;
    else
        ang_sel = theta2; val_sel = s2_val; which_s = 2;      % ✅ s2_val
    end
    csel = colPos; if ~(isfinite(val_sel) && val_sel>0), csel = colNeg; end
    cLite = lightenColor(csel, a.StressAxisLighten);
    cs = cosd(ang_sel); ss = sind(ang_sel);
    hsel = plot(ax, smax*[-cs cs], smax*[-ss ss], ...
        'LineStyle',a.StressAxisLineStyle,'LineWidth',a.StressAxisLineWidth,'Color',cLite);
    try, set(hsel,'Color',[cLite a.StressAxisAlpha]); end 
else
    which_s = 1; val_sel = s1;
end

% header lines
innerFS = get(ax,'FontSize'); if ~isempty(a.InnerFontSize), innerFS = a.InnerFontSize; end
tcur = t_sel; sec = round(second(tcur)); if sec==60, tcur=tcur+seconds(1); sec=0; end
line1 = sprintf('%s %02d''''', datestr(tcur,'yyyy-mm-dd HH:MM'), sec);
line2 = sprintf('s_1=%s Pa, s_2=%s Pa, \\theta=%.1f^\\circ', ...
                fmt_sig2_fixed(s1), fmt_sig2_fixed(s2_val), theta_deg);   % ✅ s2_val
line3 = sprintf('\\delta\\theta=%.1f^\\circ (s_{%d}, %s)', ...
                delta_sel, which_s, ternary(val_sel>0,'positive','negative'));
core = sprintf('\\fontname{Times New Roman}%s  |  %s', line1, line2);
labelStr = sprintf('%s\n\\fontname{Times New Roman}%s', core, line3);
text(ax, 0, smax*headroom*0.98, labelStr, ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'Interpreter','tex','FontSize',innerFS,'FontWeight','bold');

% export
figPath = a.FigurePath;
[outDir,~,~] = fileparts(figPath);
if ~isempty(outDir) && ~exist(outDir,'dir'), mkdir(outDir); end
if exist('exportgraphics','file')
    exportgraphics(fig, figPath, 'BackgroundColor','white','Resolution',300);
else
    set(fig,'InvertHardcopy','off'); print(fig, figPath,'-dpng','-r1200');
end
% close(fig);

% output
out.t = tt; out.sel_index = j; out.sel_time = t_sel;
out.sx_bp = sx_bp; out.sy_bp = sy_bp; out.sxy_bp = sxy_bp;
out.s1 = s1_all; out.s2 = s2_all; out.theta_deg = theta_deg_all;
out.az_epic_deg = az_epic_deg;
out.theta_s1_geo_deg = theta_s1_geo_deg_all(:);
out.delta_deg = delta_all(:);
out.params = a;

fprintf('Single petal saved: %s\n', figPath);
end

% ========= helpers =========
function [tDT, y] = local_bp(x, Fs, t0, a, makeFig)
[tDT, y] = bandpass_fft_cosine2( ...
    x, Fs, a.Fband_mHz, ...
    'Taper_uHz',a.Taper_uHz, 'Pad_hours',a.Pad_hours, 'Trim_hours',a.Trim_hours, ...
    'Smooth_hours',a.Smooth_hours, 'ReturnDatetime',true, 't0',t0, ...
    'Debeat',a.Debeat, 'BeatHours',a.BeatHours, ...
    'DecayFitHours',a.DecayFitHours, 'f0_mHz',a.f0_mHz, ...
    'PlotLog',makeFig, 'MakeFigure',makeFig);
end

function [s1, s2, theta_deg] = principal_stress(sx, sy, sxy)
s_avg = 0.5*(sx + sy);
R     = hypot(0.5*(sx - sy), sxy);
s1    = s_avg + R;
s2    = s_avg - R;
theta_deg = 0.5 * atan2d(2*sxy, (sx - sy));  % [-90,90]
theta_deg(theta_deg < 0) = theta_deg(theta_deg < 0) + 180; % [0,180)
end

function az = initial_bearing_deg(lat1, lon1, lat2, lon2)
phi1 = deg2rad(lat1); phi2 = deg2rad(lat2);
dlon = deg2rad(lon2 - lon1);
y = sin(dlon).*cos(phi2);
x = cos(phi1).*sin(phi2) - sin(phi1).*cos(phi2).*cos(dlon);
az = rad2deg(atan2(y,x)); az = wrapTo360(az);
end

function ang = wrapTo360(ang), ang = mod(ang,360); ang(ang<0) = ang(ang<0) + 360; end
function ang = wrapTo180(ang), ang = mod(ang + 180, 360) - 180; end

function c2 = lightenColor(c, f)
c = c(:).'; if numel(c)~=3, c = [0 0 0]; end
f = max(eps, min(1, f)); c2 = 1 - (1 - c)*f;
end

function d = axis_angle_diff(a1, a2)
d = abs(wrapTo180(a1 - a2)); if d > 90, d = 180 - d; end
end

function s = fmt_sig2_fixed(x), s = sprintf('%+.1e', x); end
function s = ternary(cond, s1, s2), if cond, s=s1; else, s=s2; end, end
