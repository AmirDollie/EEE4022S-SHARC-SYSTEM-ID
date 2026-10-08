%% migrateTwinCacheAlpha.m
% One-off migration of the twin FRF cache for the incident-amplitude fix (2026-10-08).
%
% Every cached node FRF (HaNodes) was computed by the unpatched computeSensorFRF, which returned
% the response to an incident wave of amplitude 1/alpha. The patched computeSensorFRF multiplies
% by alpha = depth omega^2 / g. Multiplying the cached NODE values by alpha(omegaNodes) therefore
% gives exactly what a fresh solve with the patched code would store, without any EMM solve.
% (Interpolation happens later, on the migrated nodes, exactly as for a fresh solve.)
%
% SAFEGUARDS
%   * Reads from  Results/twinCache_preAlphaFix  (rename the old folder to this first).
%   * Writes to   Results/twinCache               (must be empty or absent).
%   * The original cache is never modified.
%   * Each output file carries 'alphaMigration'; an input file that already has it is refused,
%     so no FRF is scaled twice.
%   * depth and g are read from each file's key.frfOptions ('WaterDepth', 'Gravity'), falling back
%     to the computeSensorFRF defaults (1.88, 9.81). This matters for V1's depth-mismatch fields.
%   * Only source 'emm' files are migrated (function-handle FRFs are never cached).
%   * interpCheck is kept: relErrT is unchanged by the fix; relErrH is approximately unchanged
%     (alpha is smooth). Flagged in alphaMigration.note.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
src = fullfile(P.results, 'twinCache_preAlphaFix');
dst = P.cache;                                         % Results/twinCache
DEFAULT_DEPTH = 1.88; DEFAULT_G = 9.81;                % computeSensorFRF defaults

if ~exist(src, 'dir')
    error('migrate:noSource', 'Rename Results/twinCache to Results/twinCache_preAlphaFix first (%s not found).', src);
end
if exist(dst, 'dir') && ~isempty(dir(fullfile(dst, 'twinFRF_*.mat')))
    error('migrate:destNotEmpty', '%s already contains twinFRF files. Refusing to overwrite.', dst);
end
if ~exist(dst, 'dir'), mkdir(dst); end

files = dir(fullfile(src, 'twinFRF_*.mat'));
fprintf('migrateTwinCacheAlpha: %d files in %s\n', numel(files), src);
nDone = 0; nSkip = 0; depthsSeen = [];
for i = 1:numel(files)
    f = fullfile(src, files(i).name);
    S = load(f);
    if isfield(S, 'alphaMigration')
        error('migrate:alreadyMigrated', '%s already carries alphaMigration. Aborting (would scale twice).', files(i).name);
    end
    if ~isfield(S, 'key') || ~isfield(S, 'HaNodes') || ~isfield(S, 'omegaNodes')
        warning('migrate:badFile', 'Skipping %s: missing key / HaNodes / omegaNodes.', files(i).name);
        nSkip = nSkip + 1; continue
    end
    if ~strcmp(S.key.source, 'emm')
        warning('migrate:notEmm', 'Skipping %s: source ''%s'' is not an EMM FRF.', files(i).name, S.key.source);
        nSkip = nSkip + 1; continue
    end
    [depth, g] = depthGravity(S.key.frfOptions, DEFAULT_DEPTH, DEFAULT_G);
    a = depth * S.omegaNodes(:).' .^ 2 / g;                 % 1 x nW
    if size(S.HaNodes, 2) ~= numel(a)
        error('migrate:shape', '%s: HaNodes has %d columns, omegaNodes %d.', files(i).name, size(S.HaNodes, 2), numel(a));
    end
    key = S.key; omegaNodes = S.omegaNodes; interpCheck = S.interpCheck; solveSeconds = S.solveSeconds; %#ok<NASGU>
    HaNodes = a .* S.HaNodes; %#ok<NASGU>
    alphaMigration = struct('date', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'sourceFile', files(i).name, ...
        'depth', depth, 'gravity', g, 'alphaRange', [min(a) max(a)], ...
        'note', ['HaNodes multiplied by alpha(omegaNodes) once (incident amplitude 1/alpha -> 1). ' ...
                 'interpCheck carried over from the unscaled solve.']); %#ok<NASGU>
    save(fullfile(dst, files(i).name), 'key', 'omegaNodes', 'HaNodes', 'interpCheck', 'solveSeconds', ...
        'alphaMigration', '-v7');
    nDone = nDone + 1;
    depthsSeen(end + 1) = depth; %#ok<AGROW>
end
fprintf('Migrated %d files, skipped %d. Depths used: %s\n', nDone, nSkip, mat2str(unique(round(depthsSeen, 6))));
fprintf('Original cache untouched in %s\n', src);

function [depth, g] = depthGravity(opts, depth, g)
for k = 1:2:numel(opts) - 1
    switch lower(char(opts{k}))
        case 'waterdepth', depth = opts{k + 1};
        case 'gravity', g = opts{k + 1};
    end
end
end