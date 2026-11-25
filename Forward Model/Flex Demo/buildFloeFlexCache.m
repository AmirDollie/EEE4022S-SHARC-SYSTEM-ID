function buildFloeFlexCache(varargin)
%BUILDFLOEFLEXCACHE  Precompute every floe the flexing demo can show, so the demo never solves the EMM live.
%
%   buildFloeFlexCache            % builds (or resumes) floeFlexCache.mat next to this file
%   buildFloeFlexCache('Force', true)
%
%   FLOES
%     grid     beta / beta0 in {0.25 0.5 1 2 4} x R / R0 in {0.5 0.75 1 1.5 2}, gamma = gamma0, nu = 0.3
%              (the A1 one-at-a-time and corner cases are all on this grid)
%     presets  physical floes converted exactly as A1 does (D = E h^3 / 12(1 - nu^2), beta = D / (rho_w g H^4),
%              gamma = rho_i h / (rho_w H), R = R_phys / H):
%                1 m PP disc (5 mm), SHARC model floe (4 cm PP plate, E = 1.5 GPa),
%                MIZ floes R = 10 m and 50 m (1 m ice, E = 6 GPa), Froude-scaled from H = 122 m
%   FREQUENCIES  29 nodes from 2.98 to 8.50 rad/s (the working band, about 0.2 rad/s spacing, as the engine)
%   TRUNCATION   production [50 10 10]
%
%   COST. 34 floes x 29 EMM solves. The log prints the time of the first floe and a projection. Progress is
%   saved after every floe, so an interrupted build resumes where it stopped.
%
%   Lives in Forward Model/Flex Demo/.

ip = inputParser; ip.addParameter('Force', false); ip.parse(varargin{:});
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..')); addpath(fullfile(here, '..', 'Animation'));
file = fullfile(here, 'floeFlexCache.mat');

H = 1.88; g = 9.81; rhoB = 1000;
p0 = [4.6985e-5, 1.4548e-3, 0.3830];
omega = linspace(2.98, 8.50, 29);
betaMult = [0.25 0.5 1 2 4];
RMult = [0.5 0.75 1 1.5 2];

% ---- floe list ----
floes = struct('name', {}, 'kind', {}, 'iB', {}, 'iR', {}, 'p', {}, 'nu', {}, 'Hdisp', {}, 'note', {});
for iB = 1:numel(betaMult)
    for iR = 1:numel(RMult)
        nm = sprintf('beta x %g, R x %g', betaMult(iB), RMult(iR));
        if betaMult(iB) == 1 && RMult(iR) == 1, nm = 'p0 basin disc (R = 0.72 m)'; end
        floes(end + 1) = struct('name', nm, 'kind', 'grid', 'iB', iB, 'iR', iR, ...
            'p', [p0(1) * betaMult(iB), p0(2), p0(3) * RMult(iR)], 'nu', 0.3, 'Hdisp', H, ...
            'note', 'basin scale'); %#ok<AGROW>
    end
end
floes(end + 1) = physFloe('1 m PP disc, 5 mm', 1.5e9, 0.005, 0.5, 0.4, 910, rhoB, H, H, 'basin scale');
floes(end + 1) = physFloe('SHARC model floe', 1.5e9, 0.04, 0.22, 0.4, 910, rhoB, H, H, ...
    'basin scale; 4 cm PP plate, E = 1.5 GPa');
floes(end + 1) = physFloe('MIZ floe, R = 10 m', 6e9, 1.0, 10, 0.3, 917, 1025, 122, 122, ...
    'Froude-scaled from 122 m depth');
floes(end + 1) = physFloe('MIZ floe, R = 50 m', 6e9, 1.0, 50, 0.3, 917, 1025, 122, 122, ...
    'Froude-scaled from 122 m depth');

% ---- build or resume ----
cache = struct('omega', omega, 'H', H, 'g', g, 'p0', p0, 'betaMult', betaMult, 'RMult', RMult, ...
    'sensorNames', {{'s2', 's6', 's14', 's26'}}, 'ref', 4, 'floes', floes, 'pat', {cell(1, numel(floes))}, ...
    'built', datestr(now)); %#ok<TNOW1,DATST>
if exist(file, 'file') && ~ip.Results.Force
    old = load(file, 'cache');
    if isequal(old.cache.omega, omega) && numel(old.cache.pat) == numel(floes) && ...
            isequal({old.cache.floes.name}, {floes.name})
        cache.pat = old.cache.pat;
        fprintf('Resuming: %d of %d floes already built.\n', nnz(~cellfun(@isempty, cache.pat)), numel(floes));
    end
end
tAll = tic; nDone = 0;
for k = 1:numel(floes)
    if ~isempty(cache.pat{k}), continue; end
    f = floes(k);
    cache.pat{k} = floeFlexPatterns(f.p, f.nu, omega, 'Depth', H);
    nDone = nDone + 1;
    left = nnz(cellfun(@isempty, cache.pat));
    fprintf('  %2d/%d  %-32s R/l_f = %5.2f   %5.1f s   (about %.0f min left)\n', k, numel(floes), f.name, ...
        f.p(3) / f.p(1)^0.25, cache.pat{k}.seconds, left * toc(tAll) / nDone / 60);
    save(file, 'cache', '-v7.3');
end
save(file, 'cache', '-v7.3');
fprintf('Cache complete: %s (%.1f min this session)\n', file, toc(tAll) / 60);
end

function f = physFloe(name, E, h, Rphys, nu, rhoI, rhoW, Hdepth, Hdisp, note)
% A1 conversion: nondimensionalise by the depth of the floe's own water (Froude similarity for the MIZ floes)
g = 9.81;
D = E * h^3 / (12 * (1 - nu^2));
p = [D / (rhoW * g * Hdepth^4), rhoI * h / (rhoW * Hdepth), Rphys / Hdepth];
f = struct('name', name, 'kind', 'preset', 'iB', NaN, 'iR', NaN, 'p', p, 'nu', nu, 'Hdisp', Hdisp, 'note', note);
end