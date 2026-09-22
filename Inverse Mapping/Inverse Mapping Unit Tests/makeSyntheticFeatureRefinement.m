function [refined, modeList, nu, M, P, N, H, g] = makeSyntheticFeatureRefinement(varargin)
%MAKESYNTHETICFEATUREREFINEMENT  Build a synthetic modalFeatureRefinement.m
%   results struct for modalFeatureVectorTester.m -- TEST-ONLY, not part
%   of the production pipeline.
%
% Mirrors the real project's label/row/kind assignment exactly (same five
% labels in the same order, same modeList row indices, same 'max'/'zero'
% kinds as modalFeatureRefinement.m's featureSpecs), but with small,
% cheap, arbitrary omega windows -- NOT the real refined feature
% locations -- so tests built on this run fast and do not depend on ever
% having run the real (expensive, nominal-truncation) pipeline.
%
% [refined, modeList, nu, M, P, N, H, g] = makeSyntheticFeatureRefinement()
%   uses default cheap sanity forward-model constants (nu=0.3, M=10, P=6,
%   N=10, H=1.88, g=9.81) and default well-separated omega centres
%   [2.0 2.5 3.0 3.5 4.0] rad/s at a 0.1 rad/s fixed grid step.
%
% Optional name-value overrides: 'Centers' (1x5), 'dOmega' (scalar),
% 'nu','M','P','N','H','g'. Useful for building two distinct variants
% (e.g. to check that a custom results file is never read stale).

p = inputParser;
addParameter(p, 'Centers', [2.0, 2.5, 3.0, 3.5, 4.0], @(x) isnumeric(x) && numel(x) == 5);
addParameter(p, 'dOmega', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'nu', 0.3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'M', 10, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'P', 6,  @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'N', 10, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'H', 1.88, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'g', 9.81, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});
o = p.Results;

labels = {'A_{2,0}', 'A_{0,1}', 'A_{1,1}', 'A_{2,1}', 'A_{0,0}'};
rows   = [3, 6, 7, 8, 1];              % same rows as the real featureSpecs
kinds  = {'max', 'max', 'max', 'zero', 'zero'};

refined = struct('label', {}, 'row', {}, 'kind', {}, 'omega', {}, 'A', {}, ...
    'featureOmega', {}, 'featureValue', {}, 'atEdge', {}, 'phaseFlipDeg', {});
for k = 1:5
    refined(k).label = labels{k};
    refined(k).row = rows(k);
    refined(k).kind = kinds{k};
    refined(k).omega = o.Centers(k) + o.dOmega * (-2:2);  % 5-point grid, centre at idx 3
    refined(k).A = complex(zeros(1, 5));
    refined(k).featureOmega = o.Centers(k);                % idx = 3, safely interior
    refined(k).featureValue = NaN;
    refined(k).atEdge = false;
    refined(k).phaseFlipDeg = NaN;
end

modeList = [0 0; 1 0; 2 0; 3 0; 4 0; 0 1; 1 1; 2 1];  % same 8-row table as the real project
nu = o.nu; M = o.M; P = o.P; N = o.N; H = o.H; g = o.g;
end