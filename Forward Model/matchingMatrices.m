%% Create the matching matrices M_n^(In), M_n^(0), M_n
%These linear algebra formulas are shown on pg 289 of his thesis (in
%appendix B.1)

function [MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg)

%used to eliminate A^(0)_n and A_n in favour of A^(In)_n and U^(n).
%usage: [MIN, M0, MMAT] = MATCHINGMATRICES(EIPRIME, E0PRIME, ENPRIME, D0, D, D0G, DG)

% From appendix:
%   MIN  = -(E0prime)^{-1} * EIprime         (M+1)x(M+1)
%   M0   = (D0*E0prime)^{-1} * D0g            (M+1)x(P+1)
%   MMAT = (D*Enprime)^{-1} * Dg              (M+1)x(P+1)

% The (D0*E0prime)^{-1} and (D*Enprime)^{-1} notation means "solve the
% square (M+1)x(M+1) system for the (M+1)x(P+1) right-hand side", done
% via mldivide (\), not a literal matrix inverse -> D0g and Dg are not
% square, so "inverting" them directly is not what's meant.
% 
% This will be used downstream as: A^(0)_n = MIN*A^(In)_n + M0*U^(n),
%   and:  A_n     = MMAT*U^(n).
%
% Inputs are all things already built and independently validated:
%  EIPRIME, E0PRIME from openWaterRadialDerivatives; ENPRIME from
%  radialMatrixDerivative; D0, D from innerProductMatrix; D0G, DG from
%  gegenbauerCrossMatrix -> this function does no new physics, purely
%  assembles them.


% What this outputs:

% MIn is (M+1) x (M+1) and tells us "If I know incident amplitudes
    %A^(In)_n, and nothing else abt the interface
    % how much of that shows up in the open-water scattered amplitudes
    % A^(0)_n?
% M0 is (M+1) x (P+1) and tells us how the open water scattered amps
% respond to the interface unknown U^(n) (the Gegenbauer-basis 
% representation of the normal velocity at the plate's submerged edge,
% still unknown at this point).
% Mmat is (M+1) x (P+1) and is the disc covered analogue A_n = Mmat * U^(n)



    rc1 = rcond(E0prime);
    if rc1 < 1e-12
        warning('matchingMatrices:illConditioned', ...
            'E0prime is nearly singular (rcond=%.2e) -- MIn may be unreliable.', rc1);
    end
    MIn = -(E0prime \ EIprime);

    A0 = D0*E0prime;
    rc2 = rcond(A0);
    if rc2 < 1e-12
        warning('matchingMatrices:illConditioned', ...
            'D0*E0prime is nearly singular (rcond=%.2e) -- M0 may be unreliable.', rc2);
    end
    M0 = A0 \ D0g;

    A = D*Enprime;
    rc3 = rcond(A);
    if rc3 < 1e-12
        warning('matchingMatrices:illConditioned', ...
            'D*Enprime is nearly singular (rcond=%.2e) -- Mmat may be unreliable.', rc3);
    end
    Mmat = A \ Dg;
end