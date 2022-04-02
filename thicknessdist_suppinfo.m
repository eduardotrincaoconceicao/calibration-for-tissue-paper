%% Input
dFile = 'local-sheet-thickness-distribution_vs_grammage.dat';

% Load experimental data.
% T1 = readtable('exptdata2cal_thickness-tissue-ISO-12625-3_grammage-all_unpressed-unrefined.dat', 'Delimiter', ',');
T2 = readtable('fibre_morphology.dat', 'Delimiter', ',');

% Fibre morphology.
fibre_length = T2.length; % um
fibre_width = T2.width; % um
fibre_thickness = T2.thickness; % um
coarseness = T2.coarseness; % mg/m
fibre_length = round(fibre_length/fibre_width); % voxel

% Sheet data.
% grammage_mean = T1.grammage_avg; % g/m2
grammage_mean = (20:10:120)'; % g/m2
Nx = 400;

% Simulation data.
acceptanceprob = 0.1;
alpha = 1;
delxy = fibre_width;
flex = repelem(5, 11);
W = repelem(8, 11);
L = repelem(1, 11);
% flex = [];
% W = [];
% L = [];
delz = fibre_thickness./(2*W + L);

%% Compute the Local Thickness of the Simulated Sheets
rng default % for reproducibility
nobs = length(grammage_mean);

% Conversion to required units.
coarseness = coarseness * 1e-3; % g/m
delxy = delxy * 1e-6; % m

ypred = zeros(Nx*Nx, nobs);
for jj = 1:nobs
    fout = forming(...
                   fibre_length, flex(jj), Nx, Nx, ...
                   'acceptanceprob', acceptanceprob, ...
                   'alpha', alpha, ...
                   'stop_criterion', 'Grammage', ...
                   'grammage', grammage_mean(jj), ...
                   'mass', 'coarseness', ...
                   'coarseness', coarseness, ...
                   'delxy', delxy, ...
                   'wall_thick', W(jj), ...
                   'lumen_thick', L(jj));
    ypred(:, jj) = thickness(fout);
end

% Compute the physical thickness.
ypred = ypred.*delz;

%% Save Results
T = table(repelem(grammage_mean, Nx*Nx), ypred(:));
writetable(T, dFile, 'Delimiter', ',')

%% Local Functions
function sheet_thickness = thickness(x)

Nz = size(x.web, 1);
sheet_thickness = zeros(length(x.sheet_top(:)), 1);
idx = x.sheet_top(:) < Nz+1;
sheet_thickness(idx) = x.sheet_bottom(idx) - x.sheet_top(idx) + 1;

end