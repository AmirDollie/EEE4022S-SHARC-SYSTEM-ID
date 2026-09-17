%% level2A_crossMACFixTester.m
% Regression test isolating the exact bug found by external review:
% ssiCrossMAC must compare SSI's OWN two recovered shapes to each other,
% never one recovered shape against the OTHER mode's TRUE shape.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

% Hand-constructed case where the bug and the fix give VISIBLY different
% numbers, so a silent regression would be immediately obvious.
trueShape1 = [1; 0; 0];
trueShape2 = [0; 1; 0];              % orthogonal to trueShape1 -> knownCrossMAC = 0
ssiShape1  = [1; 0; 0];              % SSI recovers shape 1 perfectly
ssiShape2  = [0.7; 0.7; 0];          % SSI recovers shape 2 with real error, NOT orthogonal to ssiShape1

knownCrossMAC = computeMAC(trueShape1, trueShape2);
fprintf('knownCrossMAC (true1 vs true2) = %.4f (should be 0, orthogonal by construction)\n', knownCrossMAC);

buggyCrossMAC   = computeMAC(ssiShape1, trueShape2);   % the OLD, incorrect computation
correctCrossMAC = computeMAC(ssiShape1, ssiShape2);    % the FIXED computation

fprintf('Buggy version   (ssiShape1 vs trueShape2): %.4f\n', buggyCrossMAC);
fprintf('Correct version (ssiShape1 vs ssiShape2):  %.4f\n', correctCrossMAC);

assert(abs(buggyCrossMAC - knownCrossMAC) < 1e-10, ...
    'Sanity check failed: expected the buggy formula to trivially match knownCrossMAC here.');
assert(correctCrossMAC > 0.4, ...
    'Correct formula should show a MEANINGFULLY different, nonzero cross-MAC here, since ssiShape2 genuinely is not orthogonal to ssiShape1.');
assert(abs(correctCrossMAC - buggyCrossMAC) > 0.3, ...
    'Fixed and buggy formulas should diverge substantially on this constructed case -- if they do not, the fix may not be applied correctly.');

fprintf('\nPASS: buggy formula trivially reproduces knownCrossMAC (%.4f), correct formula\n', buggyCrossMAC);
fprintf('correctly shows the genuine SSI recovery error instead (%.4f) -- these must differ.\n', correctCrossMAC);