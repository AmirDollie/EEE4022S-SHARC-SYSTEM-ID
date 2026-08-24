function Dg = gegenbauerCrossMatrix(rootVector, gammaDraught, P, normScale)
    %Dg^(0) contains <V0, Vg>
    % Dg is that between <V, Vg>
    % See appendix B.1 for more details
    % These form the cross matrices between the vertical mode basis and the
    % Gegenbauer basis.
    % This is necessary, due to:
    % 1. There exists non-zero draught, hence the vertical domains 
        % of omega and omega_0 are different
    % 2. There is a singularity of order -1/3 at the sharp edge.


    % TThis function has 2 uses, D_g^(0) and D_g

    % D_g^(0): pass rootvector = open water tower [k0,...,kM]'
        %normscale = 1
    %D_g: pass rootvector = disc covered [kappa0,...,kappaM]'
        %normscale = 1-gammaDraught (disc normalises at z = -gamma!

    % returns an (M+1) x (P+1) matrix, delta = 1/6 fixed throughout
    
% equ: (just above B.2) Dg(i,j) = Gamma(delta) * (2/(root_i*(1-gamma)))^delta ...
%                  * (2*(j-1)+delta) * (-1)^(j-1) ...
%                  * besselj(2*(j-1)+delta, root_i*(1-gamma)) / cos(root_i*normScale)
%NAMING COLLISION WARNING: every function in this codebase uses
%"gamma" as the draught parameter, so a bare call to gamma(1/6) here
%would NOT call MATLAB's Gamma function -> it would try to index the
%local variable named gamma using 1/6, and error. Use
%builtin('gamma', 1/6) explicitly to sidestep this. If you ever add
%another function that needs the Gamma special function, watch for
%this same trap.

%INTEGRATION-BOUNDS NOTE: the thesis states D_g^(0) as an integral
%over the full open-water domain [-1,0], even though the Gegenbauer
%basis is only defined on [-1,-gamma] (its weight requires z*<=1).
%This is consistent, not a typo: the Gegenbauer functions are taken
%to vanish on (-gamma,0] by the same extension convention used for
%u_n(z) elsewhere, so integrating over [-1,0] and [-1,-gamma] give
%the same value. Relevant if we ever re-verify this via quadrature
% ->integrate over [-1,-gamma], not [-1,0] literally.
    delta = 1/6;
    GammaDelta = builtin('gamma', delta);

    M1 = length(rootVector);
    Dg = zeros(M1, P+1);

    for i = 1:M1
        ri = rootVector(i);
        besselArg = ri*(1-gammaDraught);
        denom = cos(ri*normScale);    

        for jIdx = 0:P
            order = 2*jIdx + delta;
            Dg(i, jIdx+1) = GammaDelta * (2/besselArg)^delta * (2*jIdx+delta) ...
                            * (-1)^jIdx * besselj(order, besselArg) / denom;
        end
    end

end