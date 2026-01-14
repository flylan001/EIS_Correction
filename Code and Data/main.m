%% ========================================================================
%  EIS Experimental Data Visualization and Correction
%  Author: LAN
%  Date: 2025-11-12
%  ------------------------------------------------------------------------
%  This script:
%   (1) Plots raw EIS experimental data (µΩ)
%   (2) Scans parameters (a, b) for DRT-based reconstruction
%   (3) Evaluates RMSE and visualizes parameter sensitivity
%  ========================================================================


clc; clear; close all;
addpath(pwd);
%% EIS Data 
run('Code/Exp1.m')

%% ECM Part
run('Code/Exp1_ECM.m')

%% DRT Part
run('Code/Exp1_DRT.m')

%% DRT_RL Part
run('Code/Exp1_DRT_RL.m')

%% DRT_RL Part for Validation Experiment

run('Code/Exp1_DRT_RL.m')
