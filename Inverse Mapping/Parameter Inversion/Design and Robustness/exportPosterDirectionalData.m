function exportPosterDirectionalData(matFile)
%EXPORTPOSTERDIRECTIONALDATA  Write directional_data.csv for plot_poster_directional.py (no plotting in MATLAB).
%   exportPosterDirectionalData()                     newest Results/V2/V2_20*.mat
%   exportPosterDirectionalData('...\V2_<stamp>.mat') that run
%
%   One row per evaluated angle phi: the second-direction energy fraction f_2 (in %) at the first crossing of
%   the bias significance d_sys = 0.5, and of the coherence detection threshold, with the detection status
%   (detected | not detected | missing). Lives next to plotDirectionalForcingStudy.m.

if nargin < 1 || isempty(matFile)
    P = setupStudyPaths();
    L = dir(fullfile(P.results, 'V2', 'V2_20*.mat'));
    if isempty(L), error('exportPosterDirectionalData:none', 'No V2 results in %s', fullfile(P.results, 'V2')); end
    [~, k] = max([L.datenum]); matFile = fullfile(L(k).folder, L(k).name);
end
M = load(matFile);
X = M.results.crossings; TH = M.cfg.V2.thresholds;
ang = M.results.angleDeg(:);
fprintf('exportPosterDirectionalData: %s\n', matFile);

n = numel(ang); fB = NaN(n, 1); fD = NaN(n, 1); st = repmat({'missing'}, n, 1);
for i = 1:n
    onPhi = abs([X.phiDeg] - ang(i)) < 1e-6;
    xb = X(onPhi & strcmp({X.criterion}, 'd_sys') & abs([X.threshold] - TH.dSys(1)) < 1e-10);
    if ~isempty(xb) && xb(1).crossed, fB(i) = 100 * xb(1).f2; end
    xd = X(onPhi & strcmp({X.criterion}, 'detection'));
    if ~isempty(xd)
        if xd(1).crossed, fD(i) = 100 * xd(1).f2; st{i} = 'detected'; else, st{i} = 'not detected'; end
    end
end
T = table(ang, fB, fD, st, 'VariableNames', {'phi_deg', 'f2_bias_pct', 'f2_detect_pct', 'detect_status'});
writetable(T, 'directional_data.csv');
disp(T);
fprintf('wrote directional_data.csv; now run:  python plot_poster_directional.py\n');
end