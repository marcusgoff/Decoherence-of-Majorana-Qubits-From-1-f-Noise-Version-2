%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_Figure3cd.m
% Purpose: Generate Figure 3c and 3d from the paper
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
%   - Generates Figure 3c, 3d, and 3e from our paper.
%
% Dependencies:
%   - Results from Run_FigureS1.m must be run and saved in '../results'
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_3c, Figure_3d and Figure_3e in the ./results
%   - Saves entire worskpace under Figure_3_data.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all;   
close all;

fprintf('\n Running Run_Figure3cde \n');

%% Load results from Run_FigureS1
% load('../results/figure_S1_data');                    
load('../results_messy/figure_S1_data_5ns.mat');        
warning('off', 'all');                                  

%% Plot styling
MarkerSize     = 50;
tickLengthFrac = 0.59/0.247;

%% Fermi Golden Rule (FGR)
fermi_velocity    = (2*w*delta_0_muev)*a/hbar_mueVs;     % = w*a/hbar
R_max_fgr_tetron  = 2*(chain_length.*mu_muev.^2) ./ ...
    (8*3^(3/2)*hbar_mueVs.^2*fermi_velocity);


%% ------------------------------------------------------------------------
% Setup
% ------------------------------------------------------------------------
ylim_set = [0.000000040031643, 10^(0 + 0.5)];

delta_vec_GHz = delta_0_GHz*delta_factors;
T2_fgr_s      = 1./(R_max_fgr_tetron);                  % T2 from FGR

%% Capacitance parameters (FGR-based)
C_dot  = 0.445;      % fF
C_wire = 5*C_dot;

%% ------------------------------------------------------------------------
% Delta (flat) plot: T2* vs Δ/h
% ------------------------------------------------------------------------

% T2 time for experimental parameters (seconds)
T2_s = T2_fgr_s;

x_point = log10(delta_vec_GHz(3));    % GHz
y_point = log10(T2_s);

x_vec     = linspace(1 - 0.1, 3 + 0.1);
log_y_int = y_point;

currFig = figure();
currFig.Position = [111 330 250 317];

line_1 = log_y_int*ones(size(x_vec));
loglog(10.^(x_vec), 10.^(line_1), 'LineWidth', line_width);

xlabel('\Delta/h (GHz)');
ylabel('T_2^\ast (s)');

ax = gca;
ax.FontSize  = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir   = 'in';
ax.Box       = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

axis tight;
ylim(ylim_set);

ax.YTick = 10.^([-9, -6, -3, 0]);
ax.YAxis.MinorTickValues = 10.^([-8, -7, -5, -4, -3, -2, -1]);

ax.XTick = 10.^[1, 3];
ax.XAxis.MinorTickValues = 10.^[0:13];
ax.TickLength = tickLengthFrac*ax.TickLength;

% Draw border manually 
hold on;
ylim_current = ylim;
xlim_current = xlim;
loglog(xlim_current, [min(ylim_current) min(ylim_current)], 'k', 'LineWidth', ax.LineWidth); % bottom
loglog([min(xlim_current) min(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth); % left
loglog(xlim_current, [max(ylim_current) max(ylim_current)], 'k', 'LineWidth', ax.LineWidth); % top
loglog([max(xlim_current) max(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth); % right

% Highlight experimental point
loglog(10.^(x_point), 10.^(y_point), '.', 'Color', col_mat(3,:), ...
    'MarkerSize', MarkerSize, 'LineWidth', line_width);

%% Save Figure 3c

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_3c.fig');
% Save as PNG and SVG image
saveas(gcf, '../results/figure_3c.png');
saveas(gcf, '../results/figure_3c.svg');


%% ------------------------------------------------------------------------
% 1/S0 plot: T2* vs 1/S0
% ------------------------------------------------------------------------

% Recompute T2 for clarity (same as above, kept as you wrote it)
T2_fgr_s = 1./(R_max_fgr_tetron);
T2_s     = T2_fgr_s;

S_0_wire = S_0/25;                         % µeV
x_point  = log10((1/S_0_wire));            % 1/µeV plotted as µeV^{-1} scale label
y_point  = log10(T2_s);

S_0_min_log = 0.7;
S_0_max_log = 8.3;

x_vec     = linspace(-2, S_0_max_log);
log_y_int = y_point - 1*x_point;

currFig = figure();
currFig.Position = [472 279 343 323];

line_1 = x_vec + log_y_int;
loglog(10.^(x_vec), 10.^(line_1), 'LineWidth', line_width);

xlabel('1/S_0 (\mueV)');
ylabel('T_2^\ast (s)');

ax = gca;
ax.FontSize  = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir   = 'in';
ax.Box       = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

axis tight;

ylim(ylim_set);

ax.YTick = 10.^([-9, -6, -3, 0]);
ax.YAxis.MinorTickValues = 10.^([-8, -7, -5, -4, -3, -2, -1]);

ax.XTick = 10.^[3, 6, 9];
ax.XAxis.MinorTickValues = 10.^[0:13];
ax.TickLength = tickLengthFrac*ax.TickLength;


xlim([10^(S_0_min_log), 10^(S_0_max_log)]);

% Draw border manually (since Box is off)
hold on;
ylim_current = ylim;
xlim_current = xlim;
loglog(xlim_current, [min(ylim_current) min(ylim_current)], 'k', 'LineWidth', ax.LineWidth); % bottom
loglog([min(xlim_current) min(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth); % left
loglog(xlim_current, [max(ylim_current) max(ylim_current)], 'k', 'LineWidth', ax.LineWidth); % top
loglog([max(xlim_current) max(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth); % right

% Highlight experimental point
loglog(10.^(x_point), 10.^(y_point), '.', 'Color', col_mat(3,:), ...
    'MarkerSize', MarkerSize, 'LineWidth', line_width);

%% Save Figure 3d

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_3d.fig');
% Save as PNG and SVG image
saveas(gcf, '../results/figure_3d.png');
saveas(gcf, '../results/figure_3d.svg');


%% ------------------------------------------------------------------------
% 1/S0 versus C_wire plot
% ------------------------------------------------------------------------

x_point = log10(C_wire);        % fF
y_point = log10((1/S_0_wire));  % µeV

slope = 2;

C_min_log = -0.1;
C_max_log = 4.1;

x_vec     = linspace(C_min_log, C_max_log);
log_y_int = y_point - slope*x_point;

currFig = figure();
currFig.Position = [463 108 327 323];

line_1 = slope.*x_vec + log_y_int;
loglog(10.^(x_vec), 10.^(line_1), 'LineWidth', line_width);

xlabel('C_{wire} (fF)');
ylabel('1/S_0 (\mueV)');

ax = gca;
ax.FontSize  = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir   = 'in';
ax.Box       = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = tickLengthFrac*ax.TickLength;

% Current ticks/limits
ax.YTick = 10.^[3, 6, 9];
ax.YAxis.MinorTickValues = 10.^[0:13];

ax.XTick = 10.^[0, 3];
ax.XAxis.MinorTickValues = 10.^[0:13];

ylim([10^(S_0_min_log), 10^(S_0_max_log)]);
xlim([10^(C_min_log), 10^(C_max_log)]);

% Draw border manually 
hold on;
ylim_current = ylim;
xlim_current = xlim;
loglog(xlim_current, [min(ylim_current) min(ylim_current)], 'k', 'LineWidth', ax.LineWidth); % bottom
loglog([min(xlim_current) min(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth); % left
loglog(xlim_current, [max(ylim_current) max(ylim_current)], 'k', 'LineWidth', ax.LineWidth); % top
loglog([max(xlim_current) max(xlim_current)], ylim_current, 'k', 'LineWidth', ax.LineWidth); % right

% Highlight experimental point
loglog(10.^(x_point), 10.^(y_point), '.', 'Color', col_mat(3,:), ...
    'MarkerSize', MarkerSize, 'LineWidth', line_width);

%% Save Figure 3e

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_3e.fig');
% Save as PNG and SVG image
saveas(gcf, '../results/figure_3e.png');
saveas(gcf, '../results/figure_3e.svg');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Data
close all;
save('../results/figure_3_data');
