%% computeSensorFRFTester.m
% Unit tests for computeSensorFRF.m (Inverse Mapping/).
%
% Uses REAL EMM solves (one per frequency), so it is kept small: about a
% dozen solves in total (a few seconds each). Every test except 7 reuses
% the same three in-range frequencies.
%
% Groups
%   1  'displacement' equals direct precomputeDeflectionData + evaluateDeflection (exact)
%   2  'acceleration' equals -omega.^2 .* displacement; info.Hdisplacement consistent
%   3  Sensor ordering preserved (permuted sensors -> permuted rows, exact)
%   4  Single-frequency calls stacked equal the vector call (exact)
%   5  Repeated calls are deterministic (exact)
%   6  Physics sanity: mirror symmetry about the incidence axis,
%      eta(r, theta) = eta(r, -theta), for waves incident along theta = 0
%   7  Validated-alpha flag and warning fire only outside [1.7 13.85]
%   8  info contents (alpha, truncation, sensors, p, runtimes)
%   9  Error handling (no EMM solves)
%  10  Truncation option honoured; convergence at mid-band REPORTED, not asserted
%
% Lives in Inverse Mapping/Inverse Mapping Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));                                   % computeSensorFRF
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
clearvars -except thisDir
close all, clc

nPass = 0;
maxAll = @(X) max(abs(X(:)));

%% Reference floe, baseline sensors (Level 2C kappa-long [2 6 14 26], explicit)
beta = 4.6985e-5; gamma = 1.4548e-3; R = 0.3830; nu = 0.3;
depth = 1.88; g = 9.81;
p = [beta gamma R];
sensors = [0.3*R, 0;  0.3*R, pi;  0.5*R, pi;  0.9*R, 0];
omega = [4.5 6.0 7.8];                                               % alpha = 3.88, 6.90, 11.66
tStart = tic;

%% 1  displacement == direct Forward Model calls
fprintf('=== 1: displacement vs direct Forward Model ===\n');
[Hd, infoD] = computeSensorFRF(omega, sensors, p, 'Output', 'displacement');
Hdirect = complex(zeros(size(sensors, 1), numel(omega)));
for k = 1:numel(omega)
    data = precomputeDeflectionData(depth * omega(k)^2 / g, beta, gamma, R, nu, 50, 10, 10);
    for j = 1:size(sensors, 1)
        Hdirect(j, k) = evaluateDeflection(data, sensors(j, 1), sensors(j, 2));
    end
end
d1 = maxAll(Hd - Hdirect) / maxAll(Hdirect);
assert(isequal(size(Hd), [4 3]), '1: wrong size');
assert(d1 < 1e-12, sprintf('1: differs from direct calls (rel %.2e)', d1));
fprintf('PASS: max relative difference %.1e\n', d1); nPass = nPass + 1;

%% 2  acceleration == -omega^2 * displacement
fprintf('\n=== 2: acceleration = -omega^2 * displacement ===\n');
[Ha, infoA] = computeSensorFRF(omega, sensors, p);                   % default Output
assert(strcmp(infoA.output, 'acceleration'), '2: default Output should be acceleration');
d2 = maxAll(Ha - (-(omega.^2) .* Hd)) / maxAll(Ha);
assert(d2 < 1e-14, sprintf('2: acceleration relation off (rel %.2e)', d2));
assert(isequal(infoA.Hdisplacement, Hd), '2: info.Hdisplacement differs from displacement output');
assert(isequal(infoD.Hdisplacement, Hd), '2: displacement info.Hdisplacement differs from H');
% normalised spatial vector unchanged up to sign (-1)
vD = Hd(:, 2) / norm(Hd(:, 2)); vA = Ha(:, 2) / norm(Ha(:, 2));
assert(maxAll(vA + vD) < 1e-14, '2: normalised spatial vector should flip sign only');
fprintf('PASS: rel %.1e; spatial vector H/||H|| identical up to sign\n', d2); nPass = nPass + 1;

%% 3  sensor ordering
fprintf('\n=== 3: sensor ordering ===\n');
perm = [3 1 4 2];
Hp = computeSensorFRF(omega, sensors(perm, :), p);
assert(isequal(Hp, Ha(perm, :)), '3: permuted sensors do not give permuted rows');
fprintf('PASS\n'); nPass = nPass + 1;

%% 4  single-frequency calls == vector call
fprintf('\n=== 4: single vs vector frequency calls ===\n');
Hs = complex(zeros(size(Ha)));
for k = 1:numel(omega)
    Hs(:, k) = computeSensorFRF(omega(k), sensors, p);
end
assert(isequal(Hs, Ha), '4: single-frequency results differ from vector call');
Hcol = computeSensorFRF(omega([1 3]).', sensors, p);                % column omega accepted
assert(isequal(Hcol, Ha(:, [1 3])), '4: column OMEGA gives a different result');
fprintf('PASS\n'); nPass = nPass + 1;

%% 5  determinism
fprintf('\n=== 5: repeated calls ===\n');
Ha2 = computeSensorFRF(omega(2), sensors, p);                        % one solve is enough
assert(isequal(Ha2, Ha(:, 2)), '5: repeated call differs');
fprintf('PASS\n'); nPass = nPass + 1;

%% 6  mirror symmetry about the incidence axis
fprintf('\n=== 6: eta(r, theta) = eta(r, -theta) ===\n');
mir = [0.5*R, pi/3; 0.5*R, -pi/3; 0.8*R, 2.2; 0.8*R, -2.2];
Hm = computeSensorFRF(6.0, mir, p, 'Output', 'displacement');
d6 = max(abs(Hm(1) - Hm(2)), abs(Hm(3) - Hm(4))) / max(abs(Hm));
assert(d6 < 1e-10, sprintf('6: mirror symmetry broken (rel %.2e)', d6));
% and theta = 0 / pi sensors must be distinct (not a trivially symmetric field)
assert(abs(Ha(1, 2) - Ha(2, 2)) > 1e-3 * abs(Ha(1, 2)), '6: (0.3R,0) and (0.3R,pi) unexpectedly equal');
fprintf('PASS: rel %.1e\n', d6); nPass = nPass + 1;

%% 7  validated-alpha flag and warning
fprintf('\n=== 7: validated-alpha warning ===\n');
lastwarn('', '');
[~, i7a] = computeSensorFRF(omega, sensors, p);
[~, id] = lastwarn;
assert(isempty(id), sprintf('7: warning fired in range (%s)', id));
assert(~any(i7a.outsideValidated), '7: in-range frequencies flagged');
lastwarn('', '');
wHigh = sqrt(14.2 * g / depth);                                      % alpha = 14.2 > 13.85
fprintf('(the warning printed next is EXPECTED)\n');            % a disabled warning would not set lastwarn
[~, i7b] = computeSensorFRF([6.0 wHigh], sensors(1, :), p);
[~, id] = lastwarn;
assert(strcmp(id, 'computeSensorFRF:outsideValidatedAlpha'), '7: warning did not fire above the ceiling');
assert(isequal(i7b.outsideValidated, [false true]), '7: wrong frequencies flagged');
% the flag follows ValidatedAlphaRange, not a hard-coded value
lastwarn('', '');
[~, i7c] = computeSensorFRF(wHigh, sensors(1, :), p, 'ValidatedAlphaRange', [1.7 15]);
[~, id] = lastwarn;
assert(isempty(id) && ~i7c.outsideValidated, '7: custom ValidatedAlphaRange not honoured');
fprintf('PASS: omega = %.3f rad/s (alpha 14.2) flagged, in-range calls silent\n', wHigh); nPass = nPass + 1;

%% 8  info contents
fprintf('\n=== 8: info ===\n');
assert(maxAll(infoA.alpha - depth * omega.^2 / g) < 1e-14, '8: alpha wrong');
assert(isequal(infoA.truncation, [50 10 10]), '8: default truncation not recorded');
assert(isequal(infoA.sensors, sensors) && isequal(infoA.p, p), '8: sensors/p not recorded');
assert(isequal(infoA.omega, omega), '8: omega not recorded');
assert(numel(infoA.runtimeSec) == numel(omega) && all(infoA.runtimeSec > 0), '8: runtimes missing');
assert(infoA.depth == depth && infoA.gravity == g && infoA.nu == nu, '8: constants not recorded');
fprintf('PASS: mean %.2f s per EMM solve\n', mean(infoA.runtimeSec)); nPass = nPass + 1;

%% 9  error handling (all rejected before any EMM solve)
fprintf('\n=== 9: error handling ===\n');
cases = {
    'badOmega',          {[], sensors, p}
    'badOmega',          {[6 -1], sensors, p}
    'badOmega',          {[6 NaN], sensors, p}
    'badOmega',          {[6 1i], sensors, p}
    'badSensors',        {6, sensors(:, 1), p}
    'badSensors',        {6, [-0.1 0], p}
    'badSensors',        {6, [0.1 NaN], p}
    'sensorOutsideFloe', {6, [1.01*R 0], p}
    'badP',              {6, sensors, [beta gamma]}
    'badP',              {6, sensors, [-beta gamma R]}
    'badP',              {6, sensors, [beta gamma 0]}
    'badOption',         {6, sensors, p, 'Output'}
    'badOption',         {6, sensors, p, 'Output', 'velocity'}
    'badOption',         {6, sensors, p, 'Truncation', [50 10]}
    'badOption',         {6, sensors, p, 'Truncation', [50 10 0]}
    'badOption',         {6, sensors, p, 'ValidatedAlphaRange', [5 2]}
    'badOption',         {6, sensors, p, 'WaterDepth', -1}
    'badOption',         {6, sensors, p, 'Nope', 1}
    };
for c = 1:size(cases, 1)
    id = ['computeSensorFRF:', cases{c, 1}];
    try
        computeSensorFRF(cases{c, 2}{:});
        error('tester:noError', '9.%d: expected %s but no error was thrown', c, id);
    catch err
        assert(strcmp(err.identifier, id), sprintf('9.%d: expected %s, got %s (%s)', ...
            c, id, err.identifier, err.message));
    end
end
fprintf('PASS: %d error cases\n', size(cases, 1)); nPass = nPass + 1;

%% 10  truncation option honoured; convergence reported
fprintf('\n=== 10: truncation option (convergence reported, not asserted) ===\n');
[Ht, i10] = computeSensorFRF(6.0, sensors, p, 'Truncation', [60 12 12]);
assert(isequal(i10.truncation, [60 12 12]), '10: truncation not recorded');
d10 = maxAll(Ht - Ha(:, 2)) / maxAll(Ha(:, 2));
assert(d10 > 0, '10: higher truncation gave a bit-identical answer; option not reaching the Forward Model?');
fprintf('PASS: [60 12 12] vs [50 10 10] at 6 rad/s, max relative change %.2e (information only)\n', d10);
nPass = nPass + 1;

fprintf('\nAll %d groups passed in %.1f s.\n', nPass, toc(tStart));