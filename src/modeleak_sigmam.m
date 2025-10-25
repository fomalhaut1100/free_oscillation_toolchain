function T_out = modeleak_sigmam(result_sm_tmax, csv_file)
% MODELEAK_SIGMAM
% Modal nearest-neighbor matching for sigma_m using the first channel
% of result_sm_tmax.peak_freq_mHz{1}. Prints the first 40 rows (no CSV file).
%
% Printed header (English):
% Index,sigma_m observed frequency (uHz),Suspected leakage source,PREM frequency (uHz),sigma_m deviation (uHz)
% SPDX-License-Identifier: Apache-2.0
% Copyright (c) 2025 Tianxiang Ren
% Author: Tianxiang Ren (rentianxiang1100@163.com)

    % ---- Read PREM catalog ----
    prem_table = readtable(csv_file);
    assert(all(ismember({'Mode','fPREM_uHz'}, prem_table.Properties.VariableNames)), ...
        'The PREM CSV must contain columns: Mode and fPREM_uHz (in microhertz).');
    prem_modes = string(prem_table.Mode);
    prem_uHz   = prem_table.fPREM_uHz(:);

    % ---- Observed peaks for sigma_m (convert mHz -> uHz) ----
    obs_sigma_uHz = result_sm_tmax.peak_freq_mHz{1}(:) * 1000;  % sigma_m
    nObs = numel(obs_sigma_uHz);
    N = min(40, nObs);

    % ---- Nearest-neighbor matching (row by row) ----
    idx_arr       = (1:N).';
    picked_mode   = strings(N,1);
    picked_premHz = nan(N,1);
    delta_sigma   = nan(N,1);

    for i = 1:N
        s = obs_sigma_uHz(i);
        [~, idx_s] = min(abs(prem_uHz - s));
        picked_mode(i)   = prem_modes(idx_s);
        picked_premHz(i) = prem_uHz(idx_s);
        delta_sigma(i)   = s - picked_premHz(i); % signed deviation
    end

    % ---- Return table with English variable names ----
    T_out = table( ...
        idx_arr, obs_sigma_uHz(1:N), picked_mode, picked_premHz, delta_sigma, ...
        'VariableNames', {'idx','sigma_obs_uHz','leakage_source_mode','prem_uHz','sigma_delta_uHz'} ...
    );

    % ---- Console print (English header + all N rows) ----
    fprintf('Index,sigma_m observed frequency (uHz),Suspected leakage source,PREM frequency (uHz),sigma_m deviation (uHz)\n');
    for i = 1:N
        fprintf('%d,%.6f,%s,%.6f,%.6f\n', ...
            idx_arr(i), ...
            obs_sigma_uHz(i), ...
            picked_mode(i), ...
            picked_premHz(i), ...
            delta_sigma(i));
    end
end
