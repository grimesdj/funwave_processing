% Example driver for the nonlinear Bowen circulation solver.
% Requires the user's existing ripchannel functions and
% ripchannel_theory_circulation.m on the MATLAB path.

clear; close all;

runBATHY = 'spreadRip';
gridName = 'barRip1';
runWAVES = 'h10t10s10d00';

% Load the same local run structure used by the FUNWAVE preparation scripts.
info = prep_local_ripchannel_info(runBATHY,[gridName,'_',runWAVES]);
% Restore x, y and h into info if the saved run-info MAT file does not
% contain the generated bathymetry arrays.
if ~isfield(info,'h') || ~isfield(info,'x') || ~isfield(info,'y')
    info = make_ripchannel_bathy(info,0);
end

opts = struct();
opts.maxCells = 80000;
opts.gamma = 0.78;
opts.ndir = 9;
opts.nu = 2.0;             % horizontal eddy viscosity [m^2/s]
opts.bottomDrag = 0.02;    % Bowen linear drag parameter cL [m/s]
opts.hmin = 0.20;          % fixed wet/dry threshold [m]
opts.cfl = 0.30;
opts.maxIter = 60000;
opts.checkEvery = 100;
opts.tol = 1e-6;
opts.relax = 0.8;
opts.spongeFraction = 0.12;
opts.spongeStrength = 0.05;
opts.doPlot = true;
opts.verbose = true;
opts.saveFile = [gridName,'_',runWAVES,'_bowen_nonlinear.mat'];

out = ripchannel_bowen_nonlinear(info,opts);

fprintf('\nConverged: %d after %d iterations\n',out.converged,out.iterations);
if ~isempty(out.residualHistory)
    fprintf('Final normalized residual: %.3e\n',out.residualHistory(end));
end

% Save fields in a form convenient for comparison with FUNWAVE output.
x = out.x; y = out.y; h = out.h; eta = out.eta;
U = out.Ux; V = out.Uy; omega = out.omega; %#ok<NASGU>
save([gridName,'_',runWAVES,'_bowen_fields.mat'], ...
    'x','y','h','eta','U','V','omega','-v7.3');
