%% testReporting.m
% Smoke test for the Reporting suite (seconds, no EMM). Builds a throw-away study in a temporary folder:
% single-panel and 2 x 2 figures, series and sensor series, thresholds, ticks, panel labels, an in-panel
% and a shared legend, a CSV with awkward text, the PDF and PNG exports, the .mat, the summary and the
% manifest; then checks that the files exist and the figure sizes follow sharcStyle.
% It catches broken paths, MATLAB/Octave differences and accidental API changes, not visual taste.
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Reporting'));
S = sharcStyle();
nPass = 0; nFail = 0;
tmp = tempname; mkdir(tmp);

run = startStudy('ZZ', 'Reporting smoke test', tmp);

% 1  single panel: two series, a threshold, ticks, in-panel legend
[f1, a1] = sharcFigure('single');
x = 1:10;
h1 = sharcSeries(a1, x, x / 10, 1, 'DisplayName', 'slot 1');
h2 = sharcSeries(a1, x, x / 20, 2, 'Marker', 'none', 'DisplayName', 'slot 2');
ylim(a1, [0 1.2]);
sharcThreshold(a1, 'y', 0.5, 'd = 0.5');
sharcTicks(a1, 'x', [1 5 10]);
sharcLegend(a1, [h1 h2], {'slot 1', 'slot 2'}, 'northwest');
sharcPanel(a1, 'a', 'Single');
ok = isequal(get(h1, 'Color'), S.series(1, :)) && strcmp(get(h1, 'Marker'), S.markers{1}) && ...
    strcmp(get(h2, 'Marker'), 'none') && abs(get(h1, 'LineWidth') - S.lineWidth) < 1e-12;
[nPass, nFail] = check(ok, 'sharcSeries: slot colour, marker and line width from sharcStyle; overrides apply', ...
    nPass, nFail);
set(f1, 'Units', 'centimeters'); p = get(f1, 'Position');
[nPass, nFail] = check(abs(p(3) - S.width.single) < 0.05, sprintf('single figure is %.2f cm wide (style %.2f)', ...
    p(3), S.width.single), nPass, nFail);
saveSharcFigure(f1, run, 1, 'single');

% 2  2 x 2: sensor series on a log axis, thresholds as ticks, shared legend
[f2, a2] = sharcFigure('full', [2 2]);
hs = [];
for k = 1:numel(S.sensorNames)
    hs(end + 1) = sharcSensorSeries(a2(1), x, k * x, S.sensorNames{k}); %#ok<AGROW>
end
set(a2(1), 'YScale', 'log'); ylim(a2(1), [0.5 100]);
sharcTicks(a2(1), 'y', [1 10 100]);
sharcThreshold(a2(1), 'y', 2, '', 'tick');
for k = 2:4, sharcSeries(a2(k), x, sin(x), k); sharcPanel(a2(k), char('a' + k - 1)); end
sharcPanel(a2(1), 'a', 'Sensors');
kD = find(strcmp(S.sensorNames, 'lsm6dsv16x'));
ok = isequal(get(hs(kD), 'Color'), S.sensorColour('lsm6dsv16x')) && ...
    strcmp(get(hs(kD), 'DisplayName'), S.sensorLabel('lsm6dsv16x')) && ...
    any(abs(get(a2(1), 'YTick') - 2) < 1e-12);
[nPass, nFail] = check(ok, 'sharcSensorSeries: fixed sensor colour and label; threshold added as a tick', ...
    nPass, nFail);
set(f2, 'Units', 'centimeters'); h0 = get(f2, 'Position');
sharcSharedLegend(f2, a2, hs, cellfun(S.sensorLabel, S.sensorNames, 'UniformOutput', false), 3);
h1p = get(f2, 'Position');
[nPass, nFail] = check(abs(h1p(3) - S.width.full) < 0.05 && h1p(4) > h0(4), ...
    sprintf('2 x 2 figure %.2f cm wide; shared legend adds a %.2f cm band below', h1p(3), h1p(4) - h0(4)), ...
    nPass, nFail);
saveSharcFigure(f2, run, 2, 'grid');

% 3  CSV with commas, quotes and NaN
writeStudyCsv(run, 'table', {'name', 'value', 'note'}, ...
    {'plain', 1.5, 'ok'; 'with, comma', NaN, 'said "hi"'; 'flag', true, ''});
txt = fileread(fullfile(run.dir, 'ZZ_table.csv'));
ok = ~isempty(strfind(txt, '"with, comma"')) && ~isempty(strfind(txt, '"said ""hi"""')) && ...
    ~isempty(strfind(txt, 'flag,1,'));
[nPass, nFail] = check(ok, 'CSV: text with commas and quotes is quoted (RFC 4180), NaN empty, logical 0/1', ...
    nPass, nFail);

% 4  finish: .mat, summary, manifest
close(f1); close(f2);
run = finishStudy(run, struct('x', x), struct('note', 'test cfg'), {'headline line'});
expected = {'ZZ_table.csv', 'figures/ZZ_fig1_single.pdf', 'figures/ZZ_fig1_single.png', ...
    'figures/ZZ_fig2_grid.pdf', 'figures/ZZ_fig2_grid.png'};
filesOk = all(cellfun(@(f) exist(fullfile(run.dir, f), 'file') == 2, expected)) && ...
    exist([run.base '.mat'], 'file') == 2 && exist([run.base '_summary.txt'], 'file') == 2 && ...
    exist([run.base '.log'], 'file') == 2;
M = load([run.base '.mat']);
ok = filesOk && isequal(sort(M.runInfo.manifest(:)), sort(expected(:))) && isfield(M, 'cfg') && ...
    ~isempty(strfind(fileread([run.base '_summary.txt']), 'figures/ZZ_fig2_grid.pdf'));
[nPass, nFail] = check(ok, sprintf('files written; manifest in .mat and summary lists all %d deliverables', ...
    numel(expected)), nPass, nFail);

fprintf('\n%d passed, %d failed\n', nPass, nFail);
rmdirQuiet(tmp);
if nFail > 0, error('testReporting:failed', '%d check(s) failed.', nFail); end

function rmdirQuiet(d)
try, rmdir(d, 's'); catch, end
end

function [nPass, nFail] = check(ok, label, nPass, nFail)
if ok
    fprintf('  PASS  %s\n', label); nPass = nPass + 1;
else
    fprintf('  FAIL  %s\n', label); nFail = nFail + 1;
end
end