function Anchors = loadFeatureAnchors(refFile)
%LOADFEATUREANCHORS  Freeze the fixed sub-grid interpolation anchors used
%   by modalFeatureVector.m from a modalFeatureRefinement.m results file.
%
% For each feature saved in refFile (a modalFeatureRefinementResults.mat
% produced by modalFeatureRefinement.m), this extracts the three
% fine-grid omega points bracketing that feature's already-extracted
% discrete argmax/argmin: omega(idx-1), omega(idx), omega(idx+1), where
% idx is the index at which refined(w).omega exactly equals
% refined(w).featureOmega. These three numbers per feature are the FIXED
% anchors modalFeatureVector.m re-evaluates the EMM field at for every
% trial (beta,gamma,R); they are never re-derived per trial (see the
% design note at the top of modalFeatureVector.m for why).
%
% Anchors = loadFeatureAnchors(refFile) returns a struct with fields
%   nu, M, P, N, H, g   : forward-model constants frozen from refFile
%   modeList            : the (n,j) mode-index table from refFile
%   features            : 1xK struct array, one entry per feature, each
%                          with label, row (into modeList), kind
%                          ('max'/'zero'), and omega3 (the 1x3 fixed
%                          anchor bracket, ascending).
%
% Errors (not silent fallbacks) if refFile is missing, is missing an
% expected field, a feature's featureOmega does not exactly match a grid
% point in its own omega array, a feature sits at its refinement window's
% edge (no point on one side to bracket it with), or -- see
% EXPECTEDLABELORDER below -- the features are not in the exact expected
% order.
%
% EXPECTEDLABELORDER: modalFeatureVector.m's fVec(1:5) ordering has no
% meaning of its own -- it is just "whatever order refined(:) happens to
% be in" inside refFile. That is fine as long as it never silently
% drifts, since a Jacobian row's physical meaning depends on which
% feature it is. So this function refuses to proceed unless the loaded
% labels are EXACTLY this sequence, in this order:
EXPECTED_LABEL_ORDER = {'A_{2,0}', 'A_{0,1}', 'A_{1,1}', 'A_{2,1}', 'A_{0,0}'};

if ~isfile(refFile)
    error('modalFeatureVector:missingRefinementResults', ...
        ['%s not found -- run modalFeatureRefinement.m first (it must be run ' ...
         'from the same folder as modalFeatureVector.m).'], refFile);
end
S = load(refFile);
required = {'refined', 'modeList', 'nu', 'M', 'P', 'N', 'H', 'g'};
missing = required(~isfield(S, required));
if ~isempty(missing)
    error('modalFeatureVector:malformedRefinementResults', ...
        '%s is missing expected field(s): %s. Rerun modalFeatureRefinement.m.', ...
        refFile, strjoin(missing, ', '));
end

actualLabels = {S.refined.label};
if ~isequal(actualLabels, EXPECTED_LABEL_ORDER)
    error('modalFeatureVector:unexpectedFeatureOrder', ...
        ['refined(:) labels in %s are {%s}, expected exactly {%s} in this order -- ' ...
         'modalFeatureVector.m''s fVec ordering is hardcoded to this sequence and must ' ...
         'not silently drift if modalFeatureRefinement.m''s featureSpecs is ever reordered.'], ...
        refFile, strjoin(actualLabels, ', '), strjoin(EXPECTED_LABEL_ORDER, ', '));
end

Anchors.nu = S.nu; Anchors.M = S.M; Anchors.P = S.P; Anchors.N = S.N;
Anchors.H = S.H; Anchors.g = S.g; Anchors.modeList = S.modeList;

nFeat = numel(S.refined);
Anchors.features = repmat(struct('label', '', 'row', 0, 'kind', '', 'omega3', []), 1, nFeat);
for w = 1:nFeat
    r = S.refined(w);
    omegaGrid = r.omega;
    idx = find(omegaGrid == r.featureOmega, 1);
    if isempty(idx)
        error('modalFeatureVector:featureOmegaMismatch', ...
            ['refined(%d).featureOmega does not exactly match any point in ' ...
             'refined(%d).omega in %s -- the results file may be stale or ' ...
             'corrupted; rerun modalFeatureRefinement.m.'], w, w, refFile);
    end
    if idx <= 1 || idx >= numel(omegaGrid)
        error('modalFeatureVector:featureAtWindowEdge', ...
            ['%s feature sits at its refinement window edge in %s (no fine-grid ' ...
             'point on one side to bracket it) -- widen that window in ' ...
             'modalFeatureRefinement.m and rerun before using modalFeatureVector.m.'], ...
            r.label, refFile);
    end
    Anchors.features(w).label = r.label;
    Anchors.features(w).row = r.row;
    Anchors.features(w).kind = r.kind;
    Anchors.features(w).omega3 = omegaGrid(idx-1:idx+1);
end
end