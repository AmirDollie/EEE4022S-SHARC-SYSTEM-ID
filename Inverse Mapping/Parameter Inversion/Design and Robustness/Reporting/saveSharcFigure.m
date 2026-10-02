function files = saveSharcFigure(fig, run, number, slug)
%SAVESHARCFIGURE  Save a study figure as vector PDF (for LaTeX) and 300 dpi PNG, with standard names.
%   files = saveSharcFigure(fig, run, 1, 'gamma_tolerance')
%
%   Writes  <Results>/<ID>/figures/<ID>_fig<number>_<slug>.pdf  and  .png
%   (run from startStudy). The names carry no timestamp, so a rerun replaces the figure that the
%   thesis includes; the timestamped .mat and .txt of the run record which run produced it.
%   finishStudy lists every figure written during the run in its manifest (runInfo.manifest).
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Reporting/.
s = sharcStyle();
figDir = fullfile(run.dir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
base = fullfile(figDir, sprintf('%s_fig%d_%s', run.id, number, slug));
files = {[base '.pdf'], [base '.png']};
if exist('exportgraphics', 'file') && ~exist('OCTAVE_VERSION', 'builtin')
    axs = findall(fig, 'Type', 'axes');           % no interactive axes toolbar in the exported image
    for k = 1:numel(axs), try, axs(k).Toolbar.Visible = 'off'; catch, end, end
    exportgraphics(fig, files{1}, 'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(fig, files{2}, 'Resolution', s.dpi, 'BackgroundColor', 'white');
else
    print(fig, files{1}, '-dpdf');
    print(fig, files{2}, '-dpng', sprintf('-r%d', s.dpi));
end
fprintf('  figure %d saved: %s (.pdf, .png)\n', number, base);
end