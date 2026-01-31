%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_FigureS1.m
% Purpose: Generate Figure S1 from the Supplementary Information for Paper
%          "Decoherence in Majorana Qubits by 1/f Noise"
%
% Author: Marcus C. Goffage
% Date: 12-Aug-2025
% Affiliation: University of New South Wales
%
% Paper: "Decoherence in Majorana Qubits by 1/f Noise"
% Paper Authors: A. Alase^1, M. C. Goffage^2, M. C. Cassidy^2, 
%                S. N. Coppersmith^{2*}
% Affiliations:  ^1 University of Sydney
%                ^2 University of New South Wales
%                *  Corresponding Author
%
% -------------------------------------------------------------------------
% ABOUT THIS SCRIPT
% -------------------------------------------------------------------------
% Description:
%   - Generates Figure S1 from our paper.
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_2b and Figure_2b_inset in ../results
%   - Saves entire worskpace in the ./data directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%
fprintf('\n Running Run_FigureS1 \n');

%% Load Bespoke Quasiparticle Poisoning Library 
addpath('QPP_Library_submit')

warning('off', 'all');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Load Constants from Table 1 in the Supplementary Information
run load_constants.m
delta_0_muev = delta_muev;

% Run at the higher value of delta_mu:
mu_muev = mu_muev*5; 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Declare Simulation Input Parameters

% Take from run_script_delta_search.m and go from there. 

% Kitaev Chain Parameters. Note that our QPP_Library we use the Kitaev
% convention employed in our previous paper arXiv:2504.17485. This requires
% multiplying delta and w (as used in our latest paper) 
% by -1/2 and 1/2 respectively, as done below. 
delta_0 = 1*(-1/2); % dimensionless - as used by QPP_lib
w = 3.189*(1/2);  % dimensionless - as used by QPP_lib
BC = "OBC";
delta_0_GHz = 10^(-9)*delta_0_muev/(2*pi*hbar_mueVs);

init_state = '+'; % Initialise qubit in |+> state

% TLF Parameters
mu_low = 0; % dimensionless
mu_high = mu_muev./delta_muev; % dimensionless
mu_offset = 0; % dimensionless

gamma_min_ghz = 1; % GHz 
gamma_max_ghz = 2e3;  % GHz 
num_gamma_points = 36; 

% Simulation parameters for each gamma point:
num_trials = 50;  % number of trials at each gamma points 
t_final_ns = 1; % nanoseconds
t_final = t_final_ns*t_conversion*1e-3;

% Error Bar Parameters
% Number of trials in each group. Error bar lengths are taken as the
% standard error of the means of each group. 
group_size = 5; 

% Choose chain length
chain_length = 3e-6; % metres
N = round(chain_length/a);

% Choose delta-factors
delta_factors = [0.1, 1/3, 1, 3, 9];


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Convert Gamma and Chain length to Dimensionless Values
% Gamma-Sweep Variables
gamma_sim_to_ghz_conversion = (0.5*0.0598)/10; % 1/GHz
gamma_min = gamma_min_ghz*gamma_sim_to_ghz_conversion;
gamma_max = gamma_max_ghz*gamma_sim_to_ghz_conversion;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Run Simulations

% Preallocate storage
P_qpp_final      = cell(1, numel(delta_factors));
P_qpp_final_mat  = cell(1, numel(delta_factors));

for i = 1:length(delta_factors)
    fprintf('Running Simulation for Delta: %.0f GHz \n', delta_factors(i) * delta_0_GHz);

    [~, ~, ~, ~, ...
              ~, L_even_final_mat, ...
              L_even_final_vec, ~, gamma_vec] = ...
            run_gamma_sweep(mu_low, mu_high, mu_offset, w, ...
                            delta_factors(i)*delta_0, N, BC, ...
                            num_gamma_points, gamma_min, gamma_max, ...
                            t_final, num_trials, init_state);

    % Store results
    P_qpp_final{i}     = L_even_final_vec / 2; % /2 to get P_qpp per chain in tetron
    P_qpp_final_mat{i} = L_even_final_mat/2;

    fprintf('\n');
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%
% Convert dimensionless gamma_vec (1/2 switching rate in normalized units) to
% gamma_vec_ghz (switching rate in GHz)
gamma_vec_ghz = 2*10^(-9)*gamma_vec*delta_muev/hbar_mueVs;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Get error Bars
err_bars = cell(1, numel(P_qpp_final_mat));

for i = 1:numel(P_qpp_final_mat)
    [~, err_bars{i}] = get_gamma_sweep_error_bars_func(P_qpp_final_mat{i}, group_size);
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot Figure S1

% === Plot styling ===
marker_size = 12;
line_width = 4;
fig_ax_font_size = 30;
fig_border_LW = 2.3;
RGB = orderedcolors("gem");


col_mat = RGB([4, 3, 2, 1, 5], :);
blend_factor = 0.3;
col_mat_light = col_mat + (1 - col_mat) * blend_factor;

% === Create figure ===
currFig = figure();
currFig.Position = [205   341   823   505];
ax = gca;

% === Plot Curves ===
h = cell(1, numel(P_qpp_final));

t_final_s = t_final_ns*1e-9;
% Plot R_qpp = P_qpp/t_final_s *10^(-6) [MHz]
for i = 1:numel(P_qpp_final)
    h{i} = errorbar(gamma_vec_ghz, ...
                    P_qpp_final{i}*1e-6/(t_final_s), ...
                    err_bars{i}, ...
                    'o', ...
                    'Color',           col_mat(i, :), ...
                    'MarkerSize',      marker_size, ...
                    'LineWidth',       line_width * 0.7, ...
                    'MarkerFaceColor', col_mat(i, :), ...
                    'MarkerEdgeColor', 'w');
    hold on;
end

% === Log scale axes ===
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');


% === Axis settings ===
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'off';

xlim([gamma_min_ghz*0.9, gamma_max_ghz*1.1]);
%ylim([1e-5, 1e-1]);

% ==== Labels ====
xlabel('Two level fluctuator switching rate, $\Gamma$ (GHz)', 'Interpreter', 'latex');
%ylabel('P_{QPP, 5ns}');
%ylabel(sprintf('P_{QPP, %ins}', t_final_ns));
ylabel('R_{QPP}');

% Build legend labels dynamically
legend_labels = cell(1, numel(delta_factors));
for i = 1:numel(delta_factors)
    legend_labels{i} = sprintf('$\\Delta/h = %.0f ~\\rm{GHz}$', ...
                               delta_factors(i) * delta_0_GHz);
end

% Create legend from all handles
hLeg = legend([h{:}], legend_labels, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'off');


%%  Plot Fermi Golden Rule
hold on;
delta_vec = delta_factors*delta_0_muev;


delta_mu = (mu_high - mu_low);
delta_mu_muev = delta_mu*delta_0_muev;
gamma_vec_hz = gamma_vec_ghz*(1e9);

fermi_velocity = (2*w*delta_muev)*a/hbar_mueVs; % w*a/hbar

for ii = 1:length(delta_vec)
    delta_curr_muev = delta_vec(ii);
    
    % fix fermi velocity
    fgr = chain_length*delta_mu_muev.^2./(16*hbar_mueVs*fermi_velocity*delta_curr_muev) ...
        *gamma_vec_hz.*(1 + hbar_mueVs^2*gamma_vec_hz.^2/delta_curr_muev.^2).^(-3/2);

    fgr_no_fact = gamma_vec_hz.*(1 + hbar_mueVs^2.*gamma_vec_hz.^2./delta_curr_muev).^(-3/2);
    plot(gamma_vec_ghz, fgr*1e-6, 'color', col_mat(ii, :), 'LineWidth', 3, 'LineStyle', '--'); 
    % For P_qpp, plot: fgr*t_final_ns*1e-9
end

% Create legend from all handles
hLeg = legend([h{:}], legend_labels, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'off');

fermi_velocity = (2*w*delta_muev)*a/hbar_mueVs; % w*a/hbar
R_max = (chain_length.*delta_mu_muev.^2)./(8*3^(3/2)*hbar_mueVs.^2*fermi_velocity);
P_qpp_max = R_max*t_final_ns*1e-9;



hLeg = legend([h{:}], legend_labels, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'off');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure S1

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_S1.fig');

% Save as PNG image
saveas(gcf, '../results/figure_S1.png');
saveas(gcf, '../results/figure_S1.svg');



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Data
close all;
save('../results/figure_S1_data');
