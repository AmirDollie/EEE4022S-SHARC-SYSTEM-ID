function probe = probeValidRegion(forwardFcn, pRef, varargin)
%PROBEVALIDREGION  Map out, by direct evaluation, the neighbourhood of a
%   reference point in which a (strict) forward map succeeds, and suggest
%   conservative, evidence-based bounds for solveInverse.m.
%
%   This is a DIAGNOSTIC, not an optimiser. For the EMM twin it answers:
%   within what demonstrated neighbourhood of pRef can the current
%   fixed-anchor five-feature map be trusted? That evidence, not intuition,
%   then sets the inversion bounds, the twin truth and the starting guesses.
%
% USAGE
%   probe = probeValidRegion(forwardFcn, pRef)
%   probe = probeValidRegion(forwardFcn, pRef, Name, Value, ...)
%
%   For the EMM:
%     fwd = @(p) modalFeatureVector(p(1), p(2), p(3), 'Strict', true);
%     probe = probeValidRegion(fwd, [beta0; gamma0; R0], ...
%         'ExpectedErrorIDs', {'modalFeatureVector:invalidFeatureFit'}, ...
%         'ParamNames', {'beta', 'gamma', 'R'}, 'CaptureDetails', true);
%
% WHAT IT DOES
%   1. Reference: evaluates forwardFcn(pRef). This MUST succeed; any error
%      here, including an "expected" one, is re-thrown, because a failing
%      reference point is a setup error, not an out-of-region finding.
%   2. Axis probes: for each parameter j, each sign, and each relative
%      perturbation delta (ascending), evaluates x = 1 +/- delta in
%      coordinate j only (others at 1), i.e. p_j = pRef_j*(1 +/- delta).
%      By default, once a direction fails, larger deltas in that direction
%      are skipped (recorded as 'skipped'), since only the contiguous
%      region around pRef matters for the bounds.
%   3. Candidate box: in each direction, the verified extent is the
%      largest delta such that EVERY tested delta up to it succeeded
%      (contiguous from the reference outward). The candidate box half-
%      widths are BoundFraction times that (default 0.8), so the box sits
%      strictly inside the demonstrated region rather than on its last
%      verified point.
%   4. Combined probes: all 2^n corners of the candidate box are evaluated,
%      because individually valid changes do not guarantee a valid
%      simultaneous change. If any corner fails, all half-widths are
%      multiplied by ShrinkFactor (default 0.5) and the corners re-tested,
%      up to MaxShrinks times.
%   5. Suggested bounds: the box whose corners all succeeded, as scaled
%      [suggestedLB, suggestedUB] (and in physical units).
%
% WHAT THE RESULT DOES AND DOES NOT PROVE
%   A point where forwardFcn runs but returns complex or non-finite values
%   is recorded as 'invalid' with identifier
%   'probeValidRegion:nonFiniteOutput' (distinct from a feature-fit
%   failure), rather than aborting the run. A wrongly SIZED or non-numeric
%   output is still thrown (probeValidRegion:badForwardOutput): that is a
%   bug in the forward map, not a property of the region.
%
%   Validity is only known AT THE TESTED POINTS. If +10% succeeds and +20%
%   fails, the boundary lies somewhere in between, and it has not been
%   located. Likewise the box is verified at its axis points and corners,
%   not at every interior point; for the EMM, where invalidity means a
%   feature drifting out of its fixed anchor bracket as parameters change
%   smoothly, validity between tested points is expected but not proven.
%   solveInverse.m's ExpectedErrorIDs handling remains the safety net if
%   the solver finds a gap.
%
% NAME-VALUE OPTIONS (defaults in brackets)
%   'Perturbations'      ([0.01 0.02 0.05 0.10 0.20]) relative deltas,
%                        each in (0, 1). Sorted internally.
%   'ExpectedErrorIDs'   ({}) error identifiers meaning "outside the valid
%                        region" (searched through the error's cause chain,
%                        as in solveInverse.m). Any other error is
%                        re-thrown: a bug is never recorded as a boundary.
%   'ParamNames'         ({'p1', 'p2', ...}) labels for printing/records.
%   'StopAtFirstFailure' (true) skip larger deltas in a direction after
%                        its first failure. false probes every delta (more
%                        forward calls, and shows non-monotone behaviour).
%   'BoundFraction'      (0.8) in (0, 1]; candidate half-width as a
%                        fraction of the verified extent.
%   'ShrinkFactor'       (0.5) in (0, 1); box shrink per failed corner round.
%   'MaxShrinks'         (3) shrink rounds allowed after the first.
%   'CaptureDetails'     (false) call [f, details] = forwardFcn(p) and keep
%                        details. If details is a struct array with fields
%                        omega3 and omegaStar (as modalFeatureVector
%                        returns), each feature's position inside its
%                        anchor bracket is also computed:
%                          margin = min(w* - w_lo, w_hi - w*) / (w_hi - w_lo)
%                        0.5 = centred, 0 = at a bracket edge, < 0 = outside.
%                        A handle @(p) modalFeatureVector(...) passes two
%                        outputs through, so no extra wrapper is needed.
%   'EdgeMarginWarn'     (0.1) valid points whose smallest bracket margin is
%                        below this are flagged "valid but near an anchor
%                        edge" in the summary.
%   'Verbose'            (true) print progress and a summary.
%
% OUTPUT: probe struct
%   pRef, paramNames, perturbations, settings
%   reference        : the reference-point record (see points below); its
%                      f is the baseline f0
%   points           : struct array, one record per probe (reference, axis
%                      and combined, in evaluation order):
%                        stage ('reference' | 'axis' | 'combined'),
%                        paramIndex, paramName (axis only), delta (signed,
%                        axis only), round (combined only), x, p,
%                        status ('valid' | 'invalid' | 'skipped'), f,
%                        relShift (= (f - f0)./f0), errorIdentifier,
%                        errorMessage, runtimeSec, details,
%                        bracketMargin (per feature), minBracketMargin
%   axis             : struct array, one per parameter: name,
%                      verifiedLow / verifiedHigh (largest contiguous
%                      verified delta, 0 if the first probe failed),
%                      firstFailLow / firstFailHigh (NaN if none failed),
%                      nonMonotoneLow / nonMonotoneHigh (true if a larger
%                      delta succeeded AFTER a smaller one failed; only
%                      detectable with StopAtFirstFailure = false)
%   combined         : struct array, one per corner round: round,
%                      halfWidthLow, halfWidthHigh, nCorners, nValid
%   boxStatus        : 'verified' | 'noAxisNeighbourhood' (some parameter
%                      failed at its smallest delta in BOTH directions) |
%                      'cornersFailed' (no round had all corners valid)
%   suggestedLB, suggestedUB   : scaled bounds (NaN unless 'verified')
%   suggestedLBp, suggestedUBp : the same in physical units
%   nNearEdge        : number of valid points flagged near an anchor edge
%   nEvaluations, totalRuntimeSec
%
% Error classification uses the same cause-chain logic as solveInverse.m
% (duplicated here as a local function, kept identical on purpose).

%% ---------------- Parse options --------------------------------------
s = struct( ...
    'Perturbations', [0.01 0.02 0.05 0.10 0.20], ...
    'ExpectedErrorIDs', {{}}, ...
    'ParamNames', {{}}, ...
    'StopAtFirstFailure', true, ...
    'BoundFraction', 0.8, ...
    'ShrinkFactor', 0.5, ...
    'MaxShrinks', 3, ...
    'CaptureDetails', false, ...
    'EdgeMarginWarn', 0.1, ...
    'Verbose', true);
names = fieldnames(s);

if mod(numel(varargin), 2) ~= 0
    error('probeValidRegion:badOption', 'Name-value options must come in pairs.');
end
for k = 1:2:numel(varargin)
    nm = varargin{k};
    if ~(ischar(nm) || (isstring(nm) && isscalar(nm)))
        error('probeValidRegion:badOption', 'Option names must be text.');
    end
    hit = strcmpi(char(nm), names);
    if ~any(hit)
        error('probeValidRegion:badOption', 'Unknown option ''%s''.', char(nm));
    end
    s.(names{hit}) = varargin{k+1};
end

%% ---------------- Validate ------------------------------------------
if ~isa(forwardFcn, 'function_handle')
    error('probeValidRegion:badForwardFcn', 'forwardFcn must be a function handle, f = forwardFcn(p).');
end
if ~(isnumeric(pRef) && isvector(pRef) && isreal(pRef) && all(isfinite(pRef)) && all(pRef > 0))
    error('probeValidRegion:badPRef', 'pRef must be a real vector of finite, strictly positive values.');
end
pRef = pRef(:);
n = numel(pRef);
if n > 10
    error('probeValidRegion:tooManyParameters', ...
        'The corner test evaluates 2^n points; n = %d is too many for this probe.', n);
end

d = s.Perturbations;
if ~(isnumeric(d) && isvector(d) && isreal(d) && ~isempty(d) && all(isfinite(d)) && all(d > 0 & d < 1))
    error('probeValidRegion:badOption', 'Perturbations must be a non-empty vector of values in (0, 1).');
end
d = unique(d(:)).';          % ascending, no duplicates
s.Perturbations = d;

ids = s.ExpectedErrorIDs;
if ischar(ids) || (isstring(ids) && isscalar(ids))
    ids = {char(ids)};
end
if ~(iscell(ids) && all(cellfun(@(c) ischar(c) && ~isempty(c), ids)))
    error('probeValidRegion:badOption', 'ExpectedErrorIDs must be a cell array of non-empty error identifiers.');
end
s.ExpectedErrorIDs = ids(:).';

if isempty(s.ParamNames)
    s.ParamNames = arrayfun(@(j) sprintf('p%d', j), 1:n, 'UniformOutput', false);
end
if ~(iscell(s.ParamNames) && numel(s.ParamNames) == n && all(cellfun(@ischar, s.ParamNames)))
    error('probeValidRegion:badOption', 'ParamNames must be a cell array of %d char labels.', n);
end
isScalarIn = @(v, lo, hi, incHi) isnumeric(v) && isscalar(v) && isfinite(v) && v > lo && (v < hi || (incHi && v == hi));
if ~isScalarIn(s.BoundFraction, 0, 1, true)
    error('probeValidRegion:badOption', 'BoundFraction must be a scalar in (0, 1].');
end
if ~isScalarIn(s.ShrinkFactor, 0, 1, false)
    error('probeValidRegion:badOption', 'ShrinkFactor must be a scalar in (0, 1).');
end
if ~(isnumeric(s.MaxShrinks) && isscalar(s.MaxShrinks) && s.MaxShrinks >= 0 && s.MaxShrinks == round(s.MaxShrinks))
    error('probeValidRegion:badOption', 'MaxShrinks must be a non-negative integer.');
end
if ~(isnumeric(s.EdgeMarginWarn) && isscalar(s.EdgeMarginWarn) && isfinite(s.EdgeMarginWarn))
    error('probeValidRegion:badOption', 'EdgeMarginWarn must be a finite scalar.');
end
s.StopAtFirstFailure = logical(s.StopAtFirstFailure);
s.CaptureDetails = logical(s.CaptureDetails);
s.Verbose = logical(s.Verbose);

tAll = tic;
say('probeValidRegion: %d parameters, perturbations %s\n', n, mat2str(d));

%% ---------------- 1. Reference (must succeed; nothing caught) --------
t0 = tic;
if s.CaptureDetails
    [f0, det0] = forwardFcn(pRef);
else
    f0 = forwardFcn(pRef);
    det0 = [];
end
if ~(isnumeric(f0) && isvector(f0) && isreal(f0) && all(isfinite(f0)))
    error('probeValidRegion:badForwardOutput', ...
        'forwardFcn(pRef) must return a real, finite numeric vector; got %s %s.', mat2str(size(f0)), class(f0));
end
f0 = f0(:);
m = numel(f0);
ref = newRecord('reference', 0, '', NaN, NaN, ones(n, 1));
ref.status = 'valid';
ref.f = f0;
ref.relShift = zeros(m, 1);
ref.runtimeSec = toc(t0);
ref.details = det0;
[ref.bracketMargin, ref.minBracketMargin] = bracketMargin(det0);
points = ref;
say('  reference: valid (%.1fs), f0 = %s\n', ref.runtimeSec, mat2str(f0.', 7));

%% ---------------- 2. Axis probes -------------------------------------
axisInfo = struct('name', s.ParamNames, 'verifiedLow', 0, 'verifiedHigh', 0, ...
    'firstFailLow', NaN, 'firstFailHigh', NaN, 'nonMonotoneLow', false, 'nonMonotoneHigh', false);
for j = 1:n
    for sgn = [-1, +1]
        stopped = false;
        contiguous = true;
        nonMonotone = false;
        verified = 0;
        firstFail = NaN;
        for delta = d
            x = ones(n, 1);
            x(j) = 1 + sgn * delta;
            if stopped
                rec = newRecord('axis', j, s.ParamNames{j}, sgn * delta, NaN, x);
                rec.status = 'skipped';
            else
                rec = evaluatePoint(newRecord('axis', j, s.ParamNames{j}, sgn * delta, NaN, x));
                if strcmp(rec.status, 'valid')
                    if contiguous
                        verified = delta;
                    else
                        nonMonotone = true;   % valid again beyond a failure
                    end
                else
                    contiguous = false;
                    if isnan(firstFail), firstFail = delta; end
                    stopped = s.StopAtFirstFailure;
                end
            end
            points(end+1) = rec; %#ok<AGROW>
            if ~strcmp(rec.status, 'skipped')
                say('  %-6s %+6.1f%%: %-7s (%.1fs)%s\n', s.ParamNames{j}, 100 * sgn * delta, rec.status, ...
                    rec.runtimeSec, marginNote(rec));
            end
        end
        if sgn < 0
            axisInfo(j).verifiedLow = verified;  axisInfo(j).firstFailLow = firstFail;
            axisInfo(j).nonMonotoneLow = nonMonotone;
        else
            axisInfo(j).verifiedHigh = verified; axisInfo(j).firstFailHigh = firstFail;
            axisInfo(j).nonMonotoneHigh = nonMonotone;
        end
    end
end

%% ---------------- 3-4. Candidate box and corner rounds ---------------
combined = struct('round', {}, 'halfWidthLow', {}, 'halfWidthHigh', {}, 'nCorners', {}, 'nValid', {});
suggestedLB = NaN(n, 1);
suggestedUB = NaN(n, 1);

vLow = [axisInfo.verifiedLow].';
vHigh = [axisInfo.verifiedHigh].';
if any(vLow == 0 & vHigh == 0)
    boxStatus = 'noAxisNeighbourhood';
    say('  no verified neighbourhood for: %s (smallest delta failed in both directions)\n', ...
        strjoin(s.ParamNames(vLow == 0 & vHigh == 0), ', '));
else
    boxStatus = 'cornersFailed';
    hLow = s.BoundFraction * vLow;
    hHigh = s.BoundFraction * vHigh;
    nCorners = 2^n;
    for rnd = 1:(s.MaxShrinks + 1)
        nValid = 0;
        for c = 0:(nCorners - 1)
            hi = logical(bitget(c, 1:n)).';
            x = 1 - hLow;
            x(hi) = 1 + hHigh(hi);
            rec = evaluatePoint(newRecord('combined', 0, '', NaN, rnd, x));
            points(end+1) = rec; %#ok<AGROW>
            nValid = nValid + strcmp(rec.status, 'valid');
        end
        combined(end+1) = struct('round', rnd, 'halfWidthLow', hLow, 'halfWidthHigh', hHigh, ...
            'nCorners', nCorners, 'nValid', nValid); %#ok<AGROW>
        say('  corners, round %d: box %s .. %s -> %d/%d valid\n', rnd, ...
            mat2str((1 - hLow).', 5), mat2str((1 + hHigh).', 5), nValid, nCorners);
        if nValid == nCorners
            boxStatus = 'verified';
            suggestedLB = 1 - hLow;
            suggestedUB = 1 + hHigh;
            break
        end
        hLow = s.ShrinkFactor * hLow;
        hHigh = s.ShrinkFactor * hHigh;
    end
end

%% ---------------- 5. Assemble -----------------------------------------
isValid = strcmp({points.status}, 'valid');
minMargins = [points.minBracketMargin];
nNearEdge = sum(isValid & ~isnan(minMargins) & minMargins < s.EdgeMarginWarn);

probe = struct();
probe.pRef = pRef;
probe.paramNames = s.ParamNames;
probe.perturbations = d;
probe.settings = s;
probe.reference = ref;
probe.points = points;
probe.axis = axisInfo;
probe.combined = combined;
probe.boxStatus = boxStatus;
probe.suggestedLB = suggestedLB;
probe.suggestedUB = suggestedUB;
probe.suggestedLBp = suggestedLB .* pRef;
probe.suggestedUBp = suggestedUB .* pRef;
probe.nNearEdge = nNearEdge;
probe.nEvaluations = sum(~strcmp({points.status}, 'skipped'));
probe.totalRuntimeSec = toc(tAll);

printSummary();

%% ---------------- Nested functions -----------------------------------
    function rec = evaluatePoint(rec)
        % Evaluate forwardFcn at rec.x; classify the outcome. Unexpected
        % errors are re-thrown unchanged.
        t = tic;
        try
            if s.CaptureDetails
                [f, det] = forwardFcn(rec.p);
            else
                f = forwardFcn(rec.p);
                det = [];
            end
        catch err
            if ~matchErrorChain(err, s.ExpectedErrorIDs)
                rethrow(err);
            end
            rec.status = 'invalid';
            rec.errorIdentifier = err.identifier;
            rec.errorMessage = err.message;
            rec.runtimeSec = toc(t);
            return
        end
        rec.runtimeSec = toc(t);
        if ~(isnumeric(f) && isvector(f) && numel(f) == m)
            error('probeValidRegion:badForwardOutput', ...
                'forwardFcn returned %s %s at x = %s; expected a numeric vector with %d elements (as at pRef).', ...
                mat2str(size(f)), class(f), mat2str(rec.x.', 6), m);
        end
        if ~(isreal(f) && all(isfinite(f)))
            % The forward map ran but produced unusable numbers. Recorded
            % as invalid (with its own identifier, so it is never confused
            % with a feature-fit failure) rather than aborting what may be
            % a long run; the summary counts these separately.
            rec.status = 'invalid';
            rec.f = f(:);
            rec.errorIdentifier = 'probeValidRegion:nonFiniteOutput';
            rec.errorMessage = sprintf('forwardFcn returned a complex or non-finite value at x = %s.', ...
                mat2str(rec.x.', 6));
            rec.details = det;
            return
        end
        rec.status = 'valid';
        rec.f = f(:);
        rec.relShift = (f(:) - f0) ./ f0;
        rec.details = det;
        [rec.bracketMargin, rec.minBracketMargin] = bracketMargin(det);
    end

    function rec = newRecord(stage, j, name, delta, rnd, x)
        rec = struct('stage', stage, 'paramIndex', j, 'paramName', name, 'delta', delta, ...
            'round', rnd, 'x', x(:), 'p', x(:) .* pRef, 'status', '', 'f', [], 'relShift', [], ...
            'errorIdentifier', '', 'errorMessage', '', 'runtimeSec', NaN, 'details', [], ...
            'bracketMargin', [], 'minBracketMargin', NaN);
    end

    function txt = marginNote(rec)
        txt = '';
        if strcmp(rec.status, 'invalid')
            txt = ['  [', rec.errorIdentifier, ']'];
        elseif ~isnan(rec.minBracketMargin)
            txt = sprintf('  min bracket margin %.2f', rec.minBracketMargin);
            if rec.minBracketMargin < s.EdgeMarginWarn
                txt = [txt, '  <-- near anchor edge'];
            end
        end
    end

    function say(varargin)
        if s.Verbose
            fprintf(varargin{:});
        end
    end

    function printSummary()
        if ~s.Verbose, return; end
        fprintf('\n  Axis summary (validity known only at tested points):\n');
        for jj = 1:n
            a = axisInfo(jj);
            fprintf('    %-6s verified %6.1f%% .. %+6.1f%%   first failure: %s / %s\n', a.name, ...
                -100 * a.verifiedLow, 100 * a.verifiedHigh, pctOrNone(-a.firstFailLow), pctOrNone(a.firstFailHigh));
            if a.nonMonotoneLow || a.nonMonotoneHigh
                sides = {'-', '+'};
                fprintf('           NON-MONOTONE: valid again beyond a failure on the %s side\n', ...
                    strjoin(sides(logical([a.nonMonotoneLow, a.nonMonotoneHigh])), ' and '));
            end
        end
        nNonFinite = sum(strcmp({points.errorIdentifier}, 'probeValidRegion:nonFiniteOutput'));
        if nNonFinite > 0
            fprintf('  %d point(s) returned complex/non-finite output (recorded as invalid).\n', nNonFinite);
        end
        fprintf('  Box status: %s\n', boxStatus);
        if strcmp(boxStatus, 'verified')
            fprintf('    suggestedLB = %s\n    suggestedUB = %s   (scaled units)\n', ...
                mat2str(suggestedLB.', 6), mat2str(suggestedUB.', 6));
        end
        if nNearEdge > 0
            fprintf('  %d valid point(s) have a feature within %.2f of an anchor-bracket edge.\n', ...
                nNearEdge, s.EdgeMarginWarn);
        end
        fprintf('  %d forward evaluations, %.1fs total\n\n', probe.nEvaluations, probe.totalRuntimeSec);
    end
end

%% ------------------------------------------------------------------------
function txt = pctOrNone(v)
if isnan(v)
    txt = 'none';
else
    txt = sprintf('%+.1f%%', 100 * v);
end
end

function [margins, minMargin] = bracketMargin(det)
% Position of each feature's vertex inside its anchor bracket, if the
% details look like modalFeatureVector's (fields omega3, omegaStar).
margins = [];
minMargin = NaN;
if ~(isstruct(det) && ~isempty(det) && isfield(det, 'omega3') && isfield(det, 'omegaStar'))
    return
end
margins = NaN(numel(det), 1);
for k = 1:numel(det)
    w = sort(det(k).omega3(:));
    width = w(end) - w(1);
    if width > 0
        margins(k) = min(det(k).omegaStar - w(1), w(end) - det(k).omegaStar) / width;
    end
end
minMargin = min(margins);
end

function [tf, allIDs] = matchErrorChain(err, ids)
% True if err, or anything in its cause chain, has an identifier in ids.
% Identical to the local function in solveInverse.m.
allIDs = {};
queue = {err};
while ~isempty(queue)
    e = queue{1};
    queue(1) = [];
    allIDs{end+1} = e.identifier; %#ok<AGROW>
    try
        c = e.cause;
    catch
        c = {};
    end
    if iscell(c)
        queue = [queue, c(:).']; %#ok<AGROW>
    end
end
tf = any(ismember(allIDs, ids));
end