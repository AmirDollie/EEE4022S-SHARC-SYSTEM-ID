function zeta = evaluateSpectralDeflection(specData, r, theta, t)
%EVALUATESPECTRALDEFLECTION Physical, time-varying plate deflection under
%JONSWAP forcing, from precomputed spectral data. DIMENSIONALISED!!!!

%   SPECDATA: struct from precomputeSpectralData.m
%   R, THETA: position on the plate (must satisfy R <= specData.R)
%   T:         real PHYSICAL time (seconds), scalar or vector
%
%   Computes:
%       zeta(r,theta,t) = sum_i a_i * Re{ eta_i(r,theta) * exp(i*(omega_i*t + epsilon_i)) }
%   i.e. each frequency bin's own complex deflection eta_i (from
%   evaluateDeflection.m), weighted by that bin's physical
%   wave amplitude a_i, given its own random phase epsilon_i, oscillated
%   at its own real angular frequency omega_i, then summed across every
%   bin. This is exact (not an approximation) because the underlying
%   model is linear -> see precomputeSpectralData.m for my same note.
%
%   T may be a row vector to get a whole time series in one call
%  R, THETA must be scalars (evaluate at one spatial point per call).

    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(thisDir, '..'));
    addpath(fullfile(thisDir, '..', 'Animation'));
    
    %ensure we never evaluate a point outside of the circular floe
    if r > specData.R
        warning('evaluateSpectralDeflection:outsideDisc', ...
            'r=%.4f is outside the plate radius R=%.4f -- result is not physically meaningful (see the earlier r>R blow-up discussion).', ...
            r, specData.R);
    end

    numBins = length(specData.binData);
    t = t(:).';   % ensure row vector, so output is a time series if T has >1 element

    zeta = zeros(size(t));
    for i = 1:numBins
        eta_i = evaluateDeflection(specData.binData{i}, r, theta);
        zeta = zeta + specData.a(i) * real(eta_i * exp(1i*(specData.omega(i)*t + specData.epsilon(i))));
    end
end