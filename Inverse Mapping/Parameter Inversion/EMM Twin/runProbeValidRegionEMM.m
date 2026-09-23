%% runProbeValidRegionEMM.m
% Map the neighbourhood of the reference point (beta0, gamma0, R0) in which
% the five-feature map modalFeatureVector(...,'Strict',true) remains valid,
% for the anchor mode chosen in ANCHOR_MODE below:
%   'fixed'     : the frozen anchors (first run: R valid < 1%).
%   'predicted' : Jacobian-predicted anchors (predictFeatureAnchors.m),
%                 using Inverse Mapping/anchorPredictionModel.mat.
% Its output sets the bounds, the twin truth p_true and the starting
% guesses for runEMMTwinInversion.m.
%
% EXPENSIVE: one modalFeatureVector call per probe point (~47 points with
% StopAtFirstFailure = false: 1 reference + 30 axis + 8 corners per round,
% possibly several rounds). Per-point runtimes print as it goes.
%
% SETTINGS FOR THIS FIRST RUN
%   StopAtFirstFailure = false : probe every delta, so that non-monotone
%       validity (e.g. +10% invalid but +20% valid again) is detected
%       rather than assumed away. Later probes can use early stopping.
%   CaptureDetails = true      : keep modalFeatureVector's per-feature
%       diagnostics, so we can see which feature approaches its anchor-
%       bracket edge first.
%
% OUTPUT (in Parameter Inversion/Results/)
%   probeValidRegion_<ref|pred>_<timestamp>.mat : the full probe struct
%   probeValidRegion_<ref|pred>_<timestamp>.txt : console log of the whole run
%   ('ref' = fixed anchors, 'pred' = predicted anchors)
%
% NOTE on invalid points: under 'Strict',true, modalFeatureVector stops at
% the FIRST feature (in feature order) that fails its checks, so an
% invalid point names one failing feature, not necessarily the only one.
%
% Lives in Inverse Mapping/Parameter Inversion/EMM Twin/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));         % probeValidRegion.m (Parameter Inversion)
addpath(fullfile(thisDir, '..', '..'));   % modalFeatureVector.m (Inverse Mapping); it adds the Forward Model itself

clearvars -except thisDir
close all, clc

%% ---- Choose the anchor mode here ------------------------------------
ANCHOR_MODE = 'predicted';     % 'fixed' or 'predicted'
modelFile = fullfile(thisDir, '..', '..', 'anchorPredictionModel.mat');
%% ----------------------------------------------------------------------

switch ANCHOR_MODE
    case 'fixed'
        anchorArgs = {};
        tag = 'ref';
    case 'predicted'
        if ~isfile(modelFile)   % fail now, not after the diary and the reference solve
            error('runProbeValidRegionEMM:missingModel', ...
                'Prediction model not found: %s (build it from the saved eps=1e-4 Jacobian first).', modelFile);
        end
        anchorArgs = {'AnchorMode', 'predicted', 'AnchorModel', modelFile};
        tag = 'pred';
    otherwise
        error('runProbeValidRegionEMM:badMode', 'ANCHOR_MODE must be ''fixed'' or ''predicted''.');
end

resultsDir = fullfile(thisDir, '..', 'Results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
baseName = fullfile(resultsDir, ['probeValidRegion_', tag, '_', stamp]);
diary([baseName, '.txt']);
try   % so the log is closed even if the run errors (onCleanup does not fire at the end of a script)

%% Reference point (reference/test parameter set, NOT the physical floe)
beta0  = 4.6985e-5;
gamma0 = 1.4548e-3;
R0     = 0.3830;
pRef   = [beta0; gamma0; R0];

fwd = @(p) modalFeatureVector(p(1), p(2), p(3), 'Strict', true, anchorArgs{:});

fprintf('runProbeValidRegionEMM  (%s)\n', stamp);
fprintf('anchor mode: %s', ANCHOR_MODE);
if strcmp(ANCHOR_MODE, 'predicted'), fprintf('  (model: %s)', modelFile); end
fprintf('\n');
fprintf('pRef = [beta0 gamma0 R0] = [%.6g %.6g %.6g]\n\n', pRef);

%% Probe
probe = probeValidRegion(fwd, pRef, ...
    'ExpectedErrorIDs', {'modalFeatureVector:invalidFeatureFit'}, ...
    'ParamNames', {'beta', 'gamma', 'R'}, ...
    'StopAtFirstFailure', false, ...
    'CaptureDetails', true, ...
    'Verbose', true);

save([baseName, '.mat'], 'probe');
fprintf('Saved: %s.mat\n\n', baseName);

%% EMM-specific breakdowns
pts = probe.points;
isValid = strcmp({pts.status}, 'valid');

% (1) Which feature is closest to its anchor-bracket edge, at every valid
%     axis point (margin 0.5 = centred in its bracket, 0 = at the edge).
fprintf('Closest-to-edge feature at each valid AXIS point:\n');
fprintf('  %-6s %8s   %-10s %s\n', 'param', 'delta', 'feature', 'min margin');
for k = find(isValid & strcmp({pts.stage}, 'axis'))
    [mm, idx] = min(pts(k).bracketMargin);
    fprintf('  %-6s %+7.1f%%   %-10s %.3f\n', pts(k).paramName, 100 * pts(k).delta, ...
        pts(k).details(idx).label, mm);
end

% (2) Per feature: smallest margin seen anywhere (valid points, axis or
%     corner), and where. This is "which feature runs out of anchor first".
labels = {pts(1).details.label};
nF = numel(labels);
fprintf('\nSmallest anchor-bracket margin per feature (all valid points):\n');
fprintf('  %-10s %10s   %s\n', 'feature', 'margin', 'at x = [beta gamma R]/ref');
for f = 1:nF
    best = Inf; where = [];
    for k = find(isValid)
        if pts(k).bracketMargin(f) < best
            best = pts(k).bracketMargin(f);
            where = pts(k).x;
        end
    end
    fprintf('  %-10s %10.3f   %s\n', labels{f}, best, mat2str(where.', 4));
end

% (3) Invalid points: which feature failed first, and which check(s).
fprintf('\nInvalid points (first failing feature under Strict):\n');
fprintf('  %-9s %-28s %-10s %s\n', 'stage', 'x = [beta gamma R]/ref', 'feature', 'failed check(s)');
for k = find(strcmp({pts.status}, 'invalid'))
    tok = regexp(pts(k).errorMessage, '^(.*?): sub-grid fit failed its validity check\(s\) \[(.*?)\]', 'tokens', 'once');
    if isempty(tok)
        tok = {'?', pts(k).errorIdentifier};
    end
    fprintf('  %-9s %-28s %-10s %s\n', pts(k).stage, mat2str(pts(k).x.', 4), tok{1}, tok{2});
end

% (4) Bounds to carry forward.
fprintf('\nBox status: %s\n', probe.boxStatus);
if strcmp(probe.boxStatus, 'verified')
    fprintf('  scaled   lb = %s   ub = %s\n', mat2str(probe.suggestedLB.', 5), mat2str(probe.suggestedUB.', 5));
    fprintf('  physical lb = %s\n', mat2str(probe.suggestedLBp.', 6));
    fprintf('  physical ub = %s\n', mat2str(probe.suggestedUBp.', 6));
end
fprintf('\nDone: %d forward evaluations, %.1f min.\n', probe.nEvaluations, probe.totalRuntimeSec / 60);

catch runErr
    diary('off');
    rethrow(runErr);
end
diary('off');