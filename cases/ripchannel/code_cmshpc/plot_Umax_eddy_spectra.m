clear all
close all
%-----------------------------------------------------
%% Plot the centerline maximum velocity spectra:
% a) how does eddy shedding frequency vary with
%    rip-width and wave height?
% b) how does directional spread affect variability
%    (frequency content) in rip-channel?
%-----------------------------------------------------
%
%% 0) load compiled run stats:
% 
addpath(genpath('/storage/cms/grimesdj_lab/grimesdj/git/funwave/'))
%
%% Directories:
rootDIR  = '/scratch/grimesdj/ripchannel/';
figDIR   = '/scratch/grimesdj/ripchannel/figures/';
outDIR   = '/scratch/grimesdj/ripchannel/mat_data/';
%
%% Input file with compiled stats:
fileIn = [outDIR,'BulkVelocityStats.mat'];
load(fileIn)
%
%-----------------------------------------------------
%% Plot the centerline maximum velocity spectra:
% a) how does eddy shedding frequency vary with
%    rip-width and wave height?
%-----------------------------------------------------
