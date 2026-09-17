%% computeMAC_crossFormula_regressionTest.m
% Isolated regression test for ONE formula: ssiCrossMAC must compare
% SSI's OWN two recovered shapes to each other, never one recovered
% shape against the OTHER mode's TRUE shape. This test does NOT validate
% level2A_frequencyRegimeStudy.m as a whole -- it only proves the
% corrected formula computes what it claims to on a hand-verified case.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));
clear all, close all, clc

trueShape1 = [1; 0; 0];
trueShape2 = [0; 1; 0];              % orthogonal to trueShape1 -> knownCrossMAC = 0 exactly
ssiShape1  = [1; 0; 0];              % SSI recovers shape 1 perfectly
ssiShape2  = [0.7; 0.7; 0];          % SSI recovers shape 2 with real error

knownCrossMAC = computeMAC(trueShape1, trueShape2);
fprintf('knownCrossMAC (true1 vs true2) = %.4f (exact: 0, orthogonal by construction)\n', knownCrossMAC);

buggyCrossMAC   = computeMAC(ssiShape1, trueShape2);
correctCrossMAC = computeMAC(ssiShape1, ssiShape2);

fprintf('Buggy version   (ssiShape1 vs trueShape2): %.4f (exact expected: 0.0000)\n', buggyCrossMAC);
fprintf('Correct version (ssiShape1 vs ssiShape2):  %.4f (exact expected: 0.5000)\n', correctCrossMAC);

assert(abs(buggyCrossMAC - 0) < 1e-10, ...
    'Sanity check failed: buggy formula should trivially reproduce knownCrossMAC=0 exactly here.');
assert(abs(correctCrossMAC - 0.5) < 1e-10, ...
    'Correct formula should give exactly 0.5 for this hand-constructed case.');

fprintf('\nPASS: both formulas match their exact, hand-derived expected values -- the divergence\n');
fprintf('(0.0000 vs 0.5000) is not incidental, it is the precise, provable signature of the bug.\n');