function exportPosterLayoutData(d3File, d1File)
%EXPORTPOSTERLAYOUTDATA  Write layouts_data.csv for plot_poster_layouts.py (no plotting in MATLAB).
%   exportPosterLayoutData()                 newest Results/D3/D3_20*.mat and Results/D1/D1_*.mat
%   exportPosterLayoutData(d3File, d1File)   those runs
%
%   Rows: panel (known | unknown), sensor name, r/R, theta (deg, thesis convention), is_ref, sigma (%).
%     known   : the D1 A-optimal layout for n_s = 2 and its reference; sigma_lnbeta at the known heading
%               (results.best(1, k, 1).A.sB, headline noise case)
%     unknown : the D3a robust (worst-case) optimum for n_s = 4; worst-case sigma_lnbeta over the headings
%   Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/.

P = setupStudyPaths();
if nargin < 1 || isempty(d3File), d3File = newest(fullfile(P.results, 'D3'), 'D3_20*.mat'); end
if nargin < 2 || isempty(d1File), d1File = newest(fullfile(P.results, 'D1'), 'D1_*.mat'); end
M3 = load(d3File); M1 = load(d1File);
fprintf('exportPosterLayoutData:\n  D3 %s\n  D1 %s\n', d3File, d1File);

% D1: n_s = 2 optimum at the known heading
R1 = M1.results; G1 = R1.geometry;
k = find(R1.ns == 2, 1);
b = R1.best(1, k, 1).A;
if ~isfield(b, 'sB') || ~isfinite(b.sB), error('exportPosterLayoutData:d1', 'D1 sigma for n_s = 2 not found'); end
rows = {};
for i = double(b.layout(:).')
    rows(end + 1, :) = rowFor('known', G1, i, any(double(b.ref) == i), 100 * b.sB); %#ok<AGROW>
end

% D3a: n_s = 4 robust optimum, worst case over headings
R3 = M3.results; A = R3.d3a(1); G3 = R3.geom;
pick = @(k, name) A.perNs(k).sel(strcmp({A.perNs(k).sel.name}, name));
k4 = find(R3.nsList == 4, 1);
rw = pick(k4, 'robustWorst');
ref = []; if isfield(rw, 'ref'), ref = double(rw.ref); end
for i = double(rw.layout(:).')
    rows(end + 1, :) = rowFor('unknown', G3, i, any(ref == i), 100 * rw.sBworstN); %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', {'panel', 'name', 'r_over_R', 'theta_deg', 'is_ref', 'sigma_pct'});
writetable(T, 'layouts_data.csv');
disp(T);
fprintf('wrote layouts_data.csv; now run:  python plot_poster_layouts.py\n');
end

function r = rowFor(panel, G, i, isRef, sig)
xy = G.xy(i, :) / G.Rphys;
r = {panel, G.names{i}, round(hypot(xy(1), xy(2)), 6), round(mod(atan2d(xy(2), xy(1)), 360), 4), ...
    double(isRef), round(sig, 4)};
end

function f = newest(folder, pattern)
L = dir(fullfile(folder, pattern));
if isempty(L), error('exportPosterLayoutData:none', 'No %s in %s', pattern, folder); end
[~, k] = max([L.datenum]); f = fullfile(L(k).folder, L(k).name);
end