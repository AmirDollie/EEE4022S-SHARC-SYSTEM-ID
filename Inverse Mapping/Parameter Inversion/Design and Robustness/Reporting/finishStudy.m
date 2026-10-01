function run = finishStudy(run, results, cfg, summaryLines)
%FINISHSTUDY  Close a study run: save results and cfg (.mat), write the summary (.txt), stop the log.
%   run = finishStudy(run, results, cfg, summaryLines)
%
%   summaryLines  cell array of text lines (the headline numbers and deliverables). They are
%                 printed to the command window as well, so the log ends with the same summary.
%   The manifest (every CSV in <Results>/<ID>/ and every figure in figures/ written since startStudy)
%   is listed at the end of the summary and saved in runInfo.manifest, so the .mat records exactly
%   which tables and figures this run produced.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
run.seconds = toc(run.tic);
run.manifest = runManifest(run);           % tables and figures written during this run
run.finished = datestr(now, 'yyyy-mm-dd HH:MM:SS');
header = {sprintf('%s: %s', run.id, run.title), ...
    sprintf('run %s, %s, code commit %s, %.0f s', run.stamp, run.platform, run.commit, run.seconds), ''};
lines = [header, summaryLines(:).', {''}, {'Files written by this run:'}, ...
    cellfun(@(f) ['  ' f], run.manifest, 'UniformOutput', false)];
fprintf('\n');
fprintf('%s\n', lines{:});
fid = fopen([run.base '_summary.txt'], 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
runInfo = rmfield(run, 'tic'); %#ok<NASGU>
save([run.base '.mat'], 'results', 'cfg', 'runInfo', '-v7');
fprintf('\n  saved %s.mat and _summary.txt\n', run.base);
diary('off');
end

function m = runManifest(run)
m = {};
L = [dir(fullfile(run.dir, '*.csv')); dir(fullfile(run.dir, 'figures', '*.pdf')); ...
    dir(fullfile(run.dir, 'figures', '*.png'))];
for i = 1:numel(L)
    if L(i).datenum >= run.startNum
        rel = L(i).name;
        if ~strcmpi(rel(end - 3:end), '.csv'), rel = ['figures/' rel]; end
        m{end + 1} = rel; %#ok<AGROW>
    end
end
m = sort(m);
end
