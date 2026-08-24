%% Builds the actual open water radial derivatives 
function [EIprime, E0prime] = openWaterRadialDerivatives(n, R, k)
%computes d/dr[E^(I)_n(r)] and d/dr[E^(0)_n(r)] at r=R
%eq. 2.23 ->the open-water diagonal radial matrices.
%
% usgae: [EIPRIME, E0PRIME] = OPENWATERRADIALDERIVATIVES(N, R, K)
%
%K is the open-water root vector [k0;...;kM] from dispersionRoots
%ie called with beta=0, gamma=0. Returns two (M+1)x(M+1) diagonal
%matrices: EIPRIME uses Ihat (bounded, incident-side), E0PRIME uses
%Khat (radiating, scattered-side) -> matching E^(I)_n, E^(0)_n in the summary.

    M1 = length(k);
    EIdiag = zeros(M1,1);
    E0diag = zeros(M1,1);

    for m = 1:M1
        EIdiag(m) = scaledBesselDerivative('I', n, k(m), R, R);
        E0diag(m) = scaledBesselDerivative('K', n, k(m), R, R);
    end

    EIprime = diag(EIdiag);
    E0prime = diag(E0diag);
end