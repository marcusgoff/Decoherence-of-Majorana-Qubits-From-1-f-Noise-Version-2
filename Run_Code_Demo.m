%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: code_demo.m
% Purpose: Demo of the code used in our paper
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
%   - Demo code which runs in under 10 minutes. Calculates and plots P_qpp
%     as a function of time for L = 3 microns, and for a TLF rate of 200
%     GHz, for 3 nanoseconds or untill P_{qpp} = 0.1 (whichever occurs first),
%     this corresponds to the yellow line in figure 2b in the main tex). 
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_2b in the ./results
%   - Saves entire worskpace in the ./data directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


try
    run('Code/code_demo.m');
catch ME
    disp(getReport(ME));
end


