function run = startStudy(id, titleText, resultsRoot)
%STARTSTUDY  Open a study run: results folder, timestamped base name, and a log of the command window.
%   run = startStudy('R4', 'Uncertainty in the known gamma', resultsRoot)
%
%   Creates <resultsRoot>/<id>/ and starts diary(<id>_<stamp>.log), so everything printed during the
%   run is kept. Close it with finishStudy(run, results, cfg, summaryLines), which saves the .mat and the
%   human-readable summary. Standard layout of every study folder:
%       <ID>_<stamp>.log           full command-window output
%       <ID>_<stamp>.mat           results struct, cfg, run metadata
%       <ID>_<stamp>_summary.txt   headline numbers and deliverables (what goes into the thesis)
%       <ID>_<table>.csv           tables, one per deliverable (written with writeStudyCsv)
%       figures/<ID>_figN_<slug>.pdf / .png
%   finishStudy records every table and figure written during the run in runInfo.manifest.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
run.id = id;
run.title = titleText;
run.stamp = datestr(now, 'yyyymmdd_HHMMSS');
run.dir = fullfile(resultsRoot, id);
if ~exist(run.dir, 'dir'), mkdir(run.dir); end
run.base = fullfile(run.dir, sprintf('%s_%s', id, run.stamp));
run.started = datestr(now, 'yyyy-mm-dd HH:MM:SS');
run.startNum = now - 2 / 86400;           % files written after this belong to the run (manifest)
run.tic = tic;
run.commit = gitCommit(resultsRoot);
run.platform = platformName();
diary([run.base '.log']);
diary('on');
fprintf('\n==== %s: %s ====\n', id, titleText);
fprintf('  started %s on %s, code commit %s\n', run.started, run.platform, run.commit);
fprintf('  results in %s\n\n', run.dir);
end

function c = gitCommit(d)
c = 'unknown';
try
    [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', d));
    if st == 0, c = strtrim(out); end
catch
end
end

function p = platformName()
if exist('OCTAVE_VERSION', 'builtin'), p = ['Octave ' OCTAVE_VERSION]; else, p = ['MATLAB ' version]; end
end
