%% Inner product matrices D and D_0
function Dmat = innerProductMatrix(kappaTower, gamma);

%These are derived INCORRECTLY in appendix B.1

%I have rederived the closed form relations

%   [D]_ij = integral_{-1}^{-gamma} phi_i(z)*phi_j(z) dz, closed form:

    %diagonal:
%       [D]_ii = (1/2)*[tan(kappa_i*L)/kappa_i + L*sec(kappa_i*L)^2] 

    %off diagonal:
%       [D]_ij = [kappa_i*tan(kappa_i*L) - kappa_j*tan(kappa_j*L)] / (kappa_i^2 - kappa_j^2)
%   where L = 1-gamma.


%Note that the D and D_0 are BOTH derived from this function
% If you pass in gamma = 0 you get D_0 (which should be diagonal!)
% Note this holds if V0 is orthogonal... will test this!
% Pass the actual gamma with the disc covered tower
% [kappa0;...;kappaM] (NOT including the damped modes;those are already eliminated via Ln by this point)
% which will recover D (NOT diagonal as V is not orthogonal!)

    M1 = length(kappaTower);
    L = 1-gamma;
    Dmat = zeros(M1, M1);

    for i = 1:M1
        ki = kappaTower(i);
        for j = 1:M1
            if i == j
                Dmat(i,j) = 0.5*( tan(ki*L)/ki + L/cos(ki*L)^2 );
            else
                kj = kappaTower(j);
                Dmat(i,j) = (ki*tan(ki*L) - kj*tan(kj*L)) / (ki^2 - kj^2);
            end
        end
    
    end    

end