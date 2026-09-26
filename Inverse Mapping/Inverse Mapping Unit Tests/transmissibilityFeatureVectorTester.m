%% transmissibilityFeatureVectorTester.m
% Unit tests for transmissibilityFeatureVector.m (Inverse Mapping/) and the
% stacking pair stackTransmissibility / unstackTransmissibility
% (Transmissibility/). Real EMM solves, kept to about 20 (about 20 s in MATLAB).
%
% Groups
%   1  Stacking convention: explicit index check, round trip, NaN pass-through
%   2  Direct agreement with computeSensorFRF ratios H_j / H_r
%   3  Parameter handling: ThetaIdx = [1 3] with gamma from PBase; ThetaIdx = 1:3
%   4  Reference convention: T_j^(q) = T_j^(r) / T_q^(r) (algebraic consistency
%      between references; the raw vectors are NOT expected to match)
%   5  Agreement with the saved file 5b H0 (if a sensorFRFSensitivity_*.mat exists)
%   6  Welch grid from extractTransmissibility accepted directly
%   7  Error handling
%
% Lives in Inverse Mapping/Inverse Mapping Unit Tests/.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));                                    % transmissibilityFeatureVector, computeSensorFRF
addpath(fullfile(thisDir, '..', '..', 'Transmissibility'));          % stack / unstack / extractTransmissibility
addpath(fullfile(thisDir, '..', '..', 'Forward Model'));
addpath(fullfile(thisDir, '..', '..', 'Forward Model', 'Animation'));
clearvars -except thisDir
close all, clc

nPass = 0;
maxAll = @(X) max(abs(X(:)));
p0 = [4.6985e-5, 1.4548e-3, 0.3830]; R0 = p0(3);
sensors = [0.3*R0, 0; 0.3*R0, pi; 0.5*R0, pi; 0.9*R0, 0];
w = [4.0 5.75 7.5];
tStart = tic;

%% 1  stacking convention
fprintf('=== 1: stacking convention ===\n');
T = [1+2i, 3+4i; 5+6i, 7+8i; 9+10i, 11+12i];      % m = 3, nF = 2
f = stackTransmissibility(T);
assert(isequal(f, [1; 5; 9; 2; 6; 10; 3; 7; 11; 4; 8; 12]), '1: layout is not frequency-major');
m = 3; nF = 2;
for j = 1:m
    for k = 1:nF
        assert(f((k - 1) * 2 * m + j) == real(T(j, k)) && f((k - 1) * 2 * m + m + j) == imag(T(j, k)), '1: index formula');
    end
end
% each 2m-block is one frequency: blkdiag of per-frequency blocks needs no permutation
blk = reshape(f, 2 * m, nF);
assert(isequal(blk(:, 2), [real(T(:, 2)); imag(T(:, 2))]), '1: frequency block is not contiguous');
Tr = complex(randn(3, 45), randn(3, 45));
assert(isequal(unstackTransmissibility(stackTransmissibility(Tr), 3), Tr), '1: round trip failed');
Tn = Tr; Tn(2, 7) = NaN;                           % real NaN: imag stays 0 in MATLAB
fn = stackTransmissibility(Tn);
assert(sum(isnan(fn)) == 2 && all(isnan(fn((7 - 1) * 6 + [2 5]))), '1: NaN must mark both Re and Im of that entry');
fprintf('PASS\n'); nPass = nPass + 1;

%% 2  direct agreement with computeSensorFRF
fprintf('\n=== 2: agreement with computeSensorFRF ratios ===\n');
[f0, o0] = transmissibilityFeatureVector(log(p0([1 3])), w, sensors);
H = computeSensorFRF(w, sensors, p0);
Tdirect = H(1:3, :) ./ H(4, :);
d2 = maxAll(o0.T - Tdirect) / maxAll(Tdirect);
assert(d2 < 1e-12, sprintf('2: differs from direct ratio (rel %.2e)', d2));
assert(isequal(f0, stackTransmissibility(o0.T)) && numel(f0) == 2 * 3 * numel(w), '2: stacking or size wrong');
assert(o0.reference == 4 && isequal(o0.others, 1:3), '2: default reference / others wrong');
fprintf('PASS: rel %.1e, %d features\n', d2, numel(f0)); nPass = nPass + 1;

%% 3  parameter handling
fprintf('\n=== 3: ThetaIdx and fixed gamma ===\n');
p1 = [p0(1) * 1.10, p0(2), p0(3) * 0.97];
[~, o1] = transmissibilityFeatureVector(log(p1([1 3])), w(2), sensors);
assert(max(abs(o1.p - p1) ./ p1) < 1e-14, '3: p not assembled from theta + PBase');
gOther = 2 * p0(2);
[~, o2] = transmissibilityFeatureVector(log(p1([1 3])), w(2), sensors, 'PBase', [1 gOther 1]);
assert(abs(o2.p(2) - gOther) < 1e-18 && abs(o2.p(1) - p1(1)) < 1e-18, '3: fixed gamma not taken from PBase');
[f3, o3] = transmissibilityFeatureVector(log(p1), w(2), sensors, 'ThetaIdx', 1:3);
[f1, ~] = transmissibilityFeatureVector(log(p1([1 3])), w(2), sensors);
% exp(log(gamma0)) differs from gamma0 in the last bit, so compare with a tolerance
assert(max(abs(o3.p - o1.p) ./ o1.p) < 1e-14 && maxAll(f3 - f1) < 1e-10 * maxAll(f1), '3: three-parameter path differs');
assert(maxAll(o2.T - o1.T) > 0, '3: changing gamma had no effect at all (suspicious)');
fprintf('PASS\n'); nPass = nPass + 1;

%% 4  reference convention (algebraic consistency)
fprintf('\n=== 4: reference convention ===\n');
[~, oq] = transmissibilityFeatureVector(log(p0([1 3])), w, sensors, 'Reference', 2);
assert(isequal(oq.others, [1 3 4]), '4: others for reference 2 wrong');
% T^(2)_j = T^(4)_j / T^(4)_2 for j = 1, 3; and T^(2)_4 = 1 / T^(4)_2
T4 = o0.T;                                      % rows: sensors 1, 2, 3 relative to 4
pred = [T4(1, :) ./ T4(2, :); T4(3, :) ./ T4(2, :); 1 ./ T4(2, :)];
d4 = maxAll(oq.T - pred) / maxAll(pred);
assert(d4 < 1e-12, sprintf('4: T^(q) != T^(r)/T_q^(r) (rel %.2e)', d4));
fprintf('PASS: rel %.1e\n', d4); nPass = nPass + 1;

%% 5  agreement with saved file 5b data
fprintf('\n=== 5: agreement with saved file 5b H0 ===\n');
resDir = fullfile(thisDir, '..', 'Parameter Inversion', 'Results');
d = dir(fullfile(resDir, 'sensorFRFSensitivity_*.mat'));
if isempty(d)
    fprintf('SKIPPED: no sensorFRFSensitivity_*.mat in %s\n', resDir);
else
    [~, newest] = max([d.datenum]);
    S = load(fullfile(resDir, d(newest).name));
    kk = [1 round(numel(S.sens.omega) / 2) numel(S.sens.omega)];
    [~, o5] = transmissibilityFeatureVector(log(S.sens.p0([1 3])), S.sens.omega(kk), S.sens.sensors, ...
        'PBase', S.sens.p0);
    T5 = S.sens.H0(1:3, kk) ./ S.sens.H0(4, kk);
    d5 = maxAll(o5.T - T5) / maxAll(T5);
    assert(d5 < 1e-12, sprintf('5: differs from saved 5b data (rel %.2e)', d5));
    fprintf('PASS: rel %.1e against %s\n', d5, d(newest).name);
end
nPass = nPass + 1;

%% 6  Welch grid accepted directly
fprintf('\n=== 6: Welch grid from extractTransmissibility ===\n');
est = extractTransmissibility(randn(4, 5000), 0.1, 'Reference', 4, 'SegmentLength', 512, 'Band', [3 8.5]);
assert(numel(est.omega) == 45, '6: expected 45 Welch bins in [3, 8.5] at L = 512, dt = 0.1');
[f6, o6] = transmissibilityFeatureVector(log(p0([1 3])), est.omega(1:2), sensors, 'Reference', est.reference);
assert(isequal(o6.omega, est.omega(1:2)) && isequal(o6.others, est.others), '6: grid / channel order mismatch');
assert(isequal(size(stackTransmissibility(est.T(:, 1:2))), size(f6)), '6: data and model stacks differ in size');
fprintf('PASS: model and data share omega, channel order and stacking\n'); nPass = nPass + 1;

%% 7  errors
fprintf('\n=== 7: error handling ===\n');
th = log(p0([1 3]));
cases = {
    'badTheta',   {th(1), w, sensors}
    'badTheta',   {[th(1) NaN], w, sensors}
    'badOmega',   {th, [], sensors}
    'badOmega',   {th, [5 -1], sensors}
    'badSensors', {th, w, sensors(1, :)}
    'badOption',  {th, w, sensors, 'ThetaIdx', [1 4]}
    'badOption',  {th, w, sensors, 'ThetaIdx', [1 1]}
    'badOption',  {th, w, sensors, 'PBase', [1 2]}
    'badOption',  {th, w, sensors, 'Reference', 5}
    'badOption',  {th, w, sensors, 'FRFOptions', 'Truncation'}
    'badOption',  {th, w, sensors, 'MinRefFraction', NaN}
    'badOption',  {th, w, sensors, 'MinRefFraction', -1}
    'badTheta',   {[th(1) 1000], w, sensors}
    'badOption',  {th, w, sensors, 'Nope', 1}
    };
for c = 1:size(cases, 1)
    id = ['transmissibilityFeatureVector:', cases{c, 1}];
    try
        transmissibilityFeatureVector(cases{c, 2}{:});
        error('tester:noError', '7.%d: expected %s but no error was thrown', c, id);
    catch err
        assert(strcmp(err.identifier, id), sprintf('7.%d: expected %s, got %s (%s)', c, id, err.identifier, err.message));
    end
end
% computeSensorFRF errors propagate unchanged (sensor outside a shrunken floe)
try
    transmissibilityFeatureVector([th(1), log(0.3 * R0)], w(1), sensors);
    error('tester:noError', '7: expected computeSensorFRF:sensorOutsideFloe');
catch err
    assert(strcmp(err.identifier, 'computeSensorFRF:sensorOutsideFloe'), '7: forward error not propagated');
end
fprintf('PASS: %d cases + propagated forward error\n', size(cases, 1)); nPass = nPass + 1;

fprintf('\nAll %d groups passed in %.1f s.\n', nPass, toc(tStart));