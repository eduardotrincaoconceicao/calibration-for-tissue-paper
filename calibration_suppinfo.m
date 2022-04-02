%% Input
% sFile = 'flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-prctile-80-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined';
% dFile = 'cal_flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-prctile-80-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined.dat';
% diary flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-prctile-80-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined.log
sFile = 'flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-box-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined';
dFile = 'cal_flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-box-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined.dat';
diary flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-box-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined.log

% Load experimental data.
T1 = readtable('exptdata2cal_thickness-tissue-ISO-12625-3_grammage-all_unpressed-unrefined.dat', 'Delimiter', ',');
T2 = readtable('fibre_morphology.dat', 'Delimiter', ',');

% Fibre morphology.
fibre_length = T2.length; % um
fibre_width = T2.width; % um
fibre_thickness = T2.thickness; % um
coarseness = T2.coarseness; % mg/m
fibre_length = round(fibre_length/fibre_width); % voxel

% Sheet data.
grammage_mean = T1.grammage_avg; % g/m2
sheetthick_mean = T1.thickness_avg; % um
sheetthick_std = T1.thickness_stddev; % um
sheetthick_intrnsc = T1.thickness_correct_avg; % um
porosity_intrnsc = T1.porosity_correct;

% Simulation data.
acceptanceprob = 0.1;
alpha = 1;
delxy = fibre_width;

%% Computer Model Calibration
options = optimset('TolX', 1, 'Display', 'iter');
set_pred = 'box';
FitMethod = 'biasinvariant';
LeaveoutCrossVal = true;
repeats = 5;
rng default % for reproducibility
nobs = length(grammage_mean);

switch set_pred
    case 'box'
        selector = 'box';
    case 'prctile-80'
        selector = 'apparent';
    otherwise
        error('Unknown method %s.', set_pred)
end
switch FitMethod
    case 'leastsquares'
        fun = @(x) norm(x);
    case 'biasinvariant'
        fun = @(x) std(x);
    otherwise
        error('Unknown method %s.', FitMethod)
end

% Compute (brute-force method) terms of the Farey series of order N.
N = 20;
[flex, thick] = rat(sort(unique(cell2mat(arrayfun(@(N) (1:N-1)/N, 2:N, 'UniformOutput', false)))));
flex = [flex 2]; thick = [thick 2]; % 1/1 -> 2/2 fibres have to be at least two voxels thick
% Distribute total thickness between fibre wall and lumen.
% Note that calculating grammage from coarseness does not depend on the
% wall-lumen ratio.
W = double(idivide(int16(thick), 2, 'floor'));
n = length(thick);
L = zeros(1, n);
L(mod(thick, 2) ~= 0) = 1;

% Compute the RVE.
A3 = 800; % um3
% Source:
%   S. Rolland du Rascoat, M. Decain, X. Thibault, C. Geindreau, J.-F. Bloch (2007).
%   Estimation of microstructural properties from synchrotron X-ray
%   microtomography and determination of the REV in paper materials.
%   Acta Materialia 55(8), 2841-2850.
%   https://doi.org/10.1016/j.actamat.2006.11.050
err = 2e-4; % 0.02%
Nx = floor(sqrt( 4*(porosity_intrnsc.*(1-porosity_intrnsc)*A3)./(repeats*err^2*sheetthick_intrnsc) )/delxy);

% Conversion to required units and transformed data.
sheetthick_mean = sheetthick_mean/fibre_thickness;
coarseness = coarseness * 1e-3; % g/m
delxy = delxy * 1e-6; % m

if LeaveoutCrossVal
    for k = 1:nobs
        ind = [1:k-1 k+1:nobs];
        objective = @(idx) mean(CritVal(fun, ...
            sheetthick_mean(ind), sheetthick_std(ind), grammage_mean(ind), ...
            flex(floor(idx)), L(floor(idx)), W(floor(idx)), acceptanceprob, alpha, ...
            Nx(ind), delxy, ...
            fibre_length, coarseness, ...
            repeats, ...
            selector));
        starttime(k) = datetime;
        tic;
        [x(k), ~, exitflag(k), output(k)] = fminbnd(objective, 1, n+1, options);
        tElapsed(k) = toc;
        finishtime(k) = datetime;

        F{k} = [int2str(flex(floor(x(k)))), '/', int2str(thick(floor(x(k))))];
    end
else
    objective = @(idx) mean(CritVal(fun, ...
        sheetthick_mean, sheetthick_std, grammage_mean, ...
        flex(floor(idx)), L(floor(idx)), W(floor(idx)), acceptanceprob, alpha, ...
        Nx, delxy, ...
        fibre_length, coarseness, ...
        repeats, ...
        selector));
    starttime = datetime;
    tic;
    [x, ~, exitflag, output] = fminbnd(objective, 1, n+1, options);
    tElapsed = toc;
    finishtime = datetime;

    F = [int2str(flex(floor(x))), '/', int2str(thick(floor(x)))];
end

diary off

% -> R script 'discrepancy+figs-suppinfo.R'

% Pick solutions.
flex = flex(floor(x)); thick = thick(floor(x)); W = W(floor(x)); L = L(floor(x));
% Position of additional properties.
npar = length(x);
thickness_after_idx = npar;
rba_after_idx = 2*npar;
coverage_after_idx = 3*npar;
numcrossings_after_idx = 4*npar;

% Compute the fitted response values.
ypred = zeros(repeats, nobs, 5*npar);
for k = 1:npar
    for jj = 1:nobs
        for ii = 1:repeats
            fout = forming(...
                           fibre_length, flex(k), Nx(jj), Nx(jj), ...
                           'acceptanceprob', acceptanceprob, ...
                           'alpha', alpha, ...
                           'stop_criterion', 'Grammage', ...
                           'grammage', grammage_mean(jj), ...
                           'mass', 'coarseness', ...
                           'coarseness', coarseness, ...
                           'delxy', delxy, ...
                           'wall_thick', W(k), ...
                           'lumen_thick', L(k));
            ypred(ii, jj, k) = porosity(fout).interfiber;
            ypred(ii, jj, thickness_after_idx+k) = thickness(fout).(selector);
            ypred(ii, jj, rba_after_idx+k) = rba(fout);
            ypred(ii, jj, coverage_after_idx+k) = coverage(fout).mean;
            ypred(ii, jj, numcrossings_after_idx+k) = mean(numcrossings(fout));
        end
    end
end

ypred_mean = squeeze(mean(ypred));

%% Voxel Height
delz = fibre_thickness./thick;

% Compute the physical thickness of the simulated sheets.
ypred_mean(:, thickness_after_idx+(1:npar)) = ypred_mean(:, thickness_after_idx+(1:npar)).*delz;

%% Save Results and Input Data
save(sFile, ...
     'F', 'x', ...
     'flex', 'thick', 'W', 'L', ...
     'exitflag', 'output', ...
     'starttime', 'finishtime', 'tElapsed', ...
     'ypred', ...
     'delz', ...
     'fibre_length', 'fibre_width', 'fibre_thickness', 'coarseness', ... # input data starts here
     'acceptanceprob', 'alpha', ...
     'N', 'A3')
T = table(grammage_mean, ypred_mean);
writetable(T, dFile, 'Delimiter', ',')

%% Local Functions
function val = CritVal(lossfun, y_expt, y_expt_std, x_expt, ...
                       flex, L, W, acceptanceprob, alpha, ...
                       Nx, delxy, ...
                       fibre_length, coarseness, ...
                       repeats, ...
                       selector)

residual = ( y_expt - sheet_thickness(x_expt, ...
                                      flex, L, W, acceptanceprob, alpha, ...
                                      Nx, delxy, ...
                                      fibre_length, coarseness, ...
                                      repeats, ...
                                      selector) )./y_expt_std;
val = zeros(repeats, 1);
for jj = 1:repeats
    val(jj) = lossfun(residual(:, jj));
end

end

function ypred = sheet_thickness(grammage, ...
                                 flex, L, W, acceptanceprob, alpha, ...
                                 Nx, delxy, ...
                                 fibre_length, coarseness, ...
                                 repeats, ...
                                 selector)

nobs = size(grammage, 1);
ypred = zeros(nobs, repeats);
for jj = 1:repeats
    for ii = 1:nobs
        fout = forming(...
                       fibre_length, flex, Nx(ii), Nx(ii), ...
                       'acceptanceprob', acceptanceprob, ...
                       'alpha', alpha, ...
                       'stop_criterion', 'Grammage', ...
                       'grammage', grammage(ii), ...
                       'mass', 'coarseness', ...
                       'coarseness', coarseness, ...
                       'delxy', delxy, ...
                       'wall_thick', W, ...
                       'lumen_thick', L);
       ypred(ii, jj) = thickness(fout).(selector)/(2*W + L);
    end
end

end