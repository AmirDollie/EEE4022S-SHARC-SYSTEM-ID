function tr = studyTruncation()
%STUDYTRUNCATION  The production EMM truncation of every study: [M P N Mo].   (2026-11-26, paper plan F1)
%   The single source of truth. defineStudyScenarios sets cfg.emm.truncation from it, and emmField injects it
%   whenever a caller's FRFOptions carry no explicit 'Truncation', so every engine field, cache key and log
%   uses the same model.
%
%   [800 10 25 20000], agreed at Meeting 1 (26 Nov 2026):
%     M  = 800    plate tower; d_trunc <= 0.09 against the converged reference at p0 and three R = 2 R0 floes
%                 (Results/F1/truncationF1_20261125_124222 and _130313)
%     P  = 10     Gegenbauer order; P = 30 changes T by 3e-5 (converged)
%     N  = 25     angular orders, a RESOLUTION setting fixed from the largest radius in the search range
%                 (2 R0 x 1.3; angularOrderFor), not from the floe being estimated
%     Mo = 20000  open-water tower (split solver); Mo = 40000 changes d by 0.003-0.004
%   Before this date the studies used computeSensorFRF's default [50 10 10] (the thesis model).
    tr = [800 10 25 20000];
end
