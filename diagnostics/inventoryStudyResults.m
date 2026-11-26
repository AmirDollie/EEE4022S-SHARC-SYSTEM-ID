%% inventoryStudyResults.m  --  F10 step 0: what did each thesis study save?   (2026-11-26)
%
% WHY. The F10 refactor (classify(raw, rule)) and the reclassification under new rules are meant to be
% post-processing of saved results. That only works where a study saved the RAW engine numbers for every
% point it evaluated, not just the classes or the boundaries it found. This script reports, per study:
%   1. the structure of the saved 'results' struct (fields, types, sizes, three levels deep);
%   2. every place a class is stored, how many points it covers, and which raw fields sit next to it:
%        sigma, dthetaTot, dTot, condF, singular  (needed to reproduce the class: the regression test)
%        dthetaW, dW, dthetaSys, dSys, corr, F, sMax, f0, bW (wanted for the new rules and for F2)
%   3. a verdict: can the thesis classes be recomputed from what was saved?
% Read-only: loads the latest Results/<study>/<study>_*.mat for each study and changes nothing.
%
% RUNTIME  Under a minute.   WRITES  results/inventoryStudyResults.log
% HOW      run('diagnostics/inventoryStudyResults.m')   (from any folder)

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
root = fileparts(here);
resRoot = fullfile(root, 'Inverse Mapping', 'Parameter Inversion', 'Results');
outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = fullfile(outDir, 'inventoryStudyResults.log');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile); diary('on');
cleanupDiary = onCleanup(@() diary('off'));

studies = {'A1', 'D1', 'D2', 'D3', 'R1', 'R2', 'R3', 'R4', 'V1', 'V2', 'M1'};
classNames = {'class', 'finalClass', 'inverseClass', 'classHi', 'classLo', 'classes'};
needRaw = {'sigma', 'dthetaTot', 'dTot', 'condF', 'singular'};
wantRaw = {'dthetaW', 'dW', 'dthetaSys', 'dSys', 'corr', 'F', 'sMax', 'f0', 'bW'};

fprintf('==== inventoryStudyResults  (%s) ====\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Results folder: %s\n', resRoot);
fprintf('Needed per point to recompute a class: %s\n', strjoin(needRaw, ', '));
fprintf('Wanted for the new rules and F2:       %s\n', strjoin(wantRaw, ', '));

summary = cell(numel(studies), 4);
for s = 1:numel(studies)
    id = studies{s};
    fprintf('\n==================== %s ====================\n', id);
    files = dir(fullfile(resRoot, id, sprintf('%s_*.mat', id)));
    if isempty(files)
        fprintf('  no saved results found\n');
        summary(s, :) = {id, '-', 0, 'no file'};
        continue
    end
    [~, k] = max([files.datenum]);
    f = fullfile(files(k).folder, files(k).name);
    fprintf('  file: %s  (%.1f MB; %d run(s) on disk)\n', files(k).name, files(k).bytes / 1e6, numel(files));
    S = load(f);
    if ~isfield(S, 'results')
        fprintf('  no ''results'' variable; variables: %s\n', strjoin(fieldnames(S), ', '));
        summary(s, :) = {id, files(k).name, 0, 'no results variable'};
        continue
    end

    fprintf('\n  -- structure (3 levels) --\n');
    printTree('results', S.results, 0, 3);

    fprintf('\n  -- where classes are stored --\n');
    hits = findClassHolders(S.results, 'results', classNames, 0, 8);
    if isempty(hits)
        fprintf('  no class field found anywhere\n');
        summary(s, :) = {id, files(k).name, 0, 'no class stored'};
        continue
    end
    best = 0; bestOk = false;
    for h = 1:numel(hits)
        H = hits(h);
        haveNeed = needRaw(ismember(needRaw, H.fields));
        haveWant = wantRaw(ismember(wantRaw, H.fields));
        missNeed = setdiff(needRaw, haveNeed, 'stable');
        fprintf('  %-60s %5d point(s)   class field: %s\n', H.path, H.n, H.classField);
        fprintf('      needed present: %-40s missing: %s\n', listOr(haveNeed), listOr(missNeed));
        fprintf('      wanted present: %s\n', listOr(haveWant));
        ok = isempty(missNeed);
        if H.n > best || (H.n == best && ok), best = H.n; bestOk = ok; end
    end
    if bestOk, verdict = 'raw fields saved with the classes: regression test possible';
    else, verdict = 'classes saved WITHOUT all raw fields: needs fuller saving in the rerun'; end
    summary(s, :) = {id, files(k).name, best, verdict};
end

fprintf('\n==================== SUMMARY ====================\n');
fprintf('  %-4s %-26s %8s   %s\n', 'id', 'file', 'points', 'verdict');
for s = 1:numel(studies)
    fprintf('  %-4s %-26s %8d   %s\n', summary{s, 1}, summary{s, 2}, summary{s, 3}, summary{s, 4});
end
fprintf(['\n  Note: a study can also store boundaries (D2 windows, R1 noise limits, tolerances) found by a search\n' ...
    '  under the old rule. Check the structure listing for the sweep grids themselves; the new rules need\n' ...
    '  the raw numbers at every grid point to redo those searches by interpolation.\n']);
fprintf('==== done ====\n');
diary('off');

%% ==== local functions ================================================================================================
function printTree(name, v, depth, maxDepth)
% One line per node: name, class, size; structs and cells are expanded down to maxDepth.
    pad = repmat('  ', 1, depth + 1);
    sz = mat2str(size(v));
    if isstruct(v)
        fn = fieldnames(v);
        fprintf('%s%s  struct %s, %d field(s)\n', pad, name, sz, numel(fn));
        if depth < maxDepth && ~isempty(v)
            for i = 1:numel(fn)
                printTree(fn{i}, v(1).(fn{i}), depth + 1, maxDepth);
            end
            if numel(v) > 1, fprintf('%s  (fields of element 1 of %d)\n', pad, numel(v)); end
        end
    elseif iscell(v)
        fprintf('%s%s  cell %s', pad, name, sz);
        if ~isempty(v), fprintf(', first element %s %s', class(v{1}), mat2str(size(v{1}))); end
        fprintf('\n');
        if depth < maxDepth && ~isempty(v) && isstruct(v{1})
            printTree([name '{1}'], v{1}, depth + 1, maxDepth);
        end
    elseif ischar(v) || isstring(v)
        txt = char(v); if numel(txt) > 40, txt = [txt(1:40) '...']; end
        fprintf('%s%s  %s ''%s''\n', pad, name, class(v), txt);
    else
        fprintf('%s%s  %s %s\n', pad, name, class(v), sz);
    end
end

function hits = findClassHolders(v, path, classNames, depth, maxDepth)
% Every struct (or struct array) that has a class-like field, with the number of points it covers and
% the names of all its fields (so the raw fields stored next to the class can be listed).
    hits = struct('path', {}, 'n', {}, 'classField', {}, 'fields', {});
    if depth > maxDepth, return; end
    if isstruct(v) && ~isempty(v)
        fn = fieldnames(v);
        cf = fn(ismember(fn, classNames));
        if ~isempty(cf)
            c = v(1).(cf{1});
            if numel(v) > 1, n = numel(v);
            elseif iscell(c) || isstring(c), n = numel(c);
            else, n = 1; end
            hits(end + 1) = struct('path', path, 'n', n, 'classField', cf{1}, 'fields', {fn}); %#ok<AGROW>
        end
        for i = 1:numel(fn)
            sub = v(1).(fn{i});
            if isstruct(sub) || iscell(sub)
                p = path; if numel(v) > 1, p = [p '(i)']; end
                hits = [hits, findClassHolders(sub, [p '.' fn{i}], classNames, depth + 1, maxDepth)]; %#ok<AGROW>
            end
        end
    elseif iscell(v) && ~isempty(v)
        k = find(cellfun(@isstruct, v), 1);
        if ~isempty(k)
            sub = findClassHolders(v{k}, [path '{i}'], classNames, depth + 1, maxDepth);
            for j = 1:numel(sub)
                if sub(j).n == 1, sub(j).n = sum(cellfun(@isstruct, v(:))); end
            end
            hits = [hits, sub];
        end
    end
end

function s = listOr(c)
    if isempty(c), s = '(none)'; else, s = strjoin(c, ', '); end
end
