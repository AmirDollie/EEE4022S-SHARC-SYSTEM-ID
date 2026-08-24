function [MU, Sn, Sn0] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat)
%SCATTERINGMATRICES M_n^(U) and the scattering matrices S_n, S_n^(0)
%Appendix B.1, eq. B.1c onward) -> the payoff of the entire matching system

%usage: [MU, SN, SN0] = SCATTERINGMATRICES(D0G, DG, EN_R, MIN, M0, MMAT)

%   MU  = [Dg^T*En_R*Mmat - D0g^T*M0]^{-1} * D0g^T*(I + MIn)   (P+1)x(M+1)
%   SN  = Mmat * MU                                             (M+1)x(M+1)
%   SN0 = MIn + M0 * MU                                         (M+1)x(M+1)
%
%Usage downstream: An = Sn*AIn, A0n = Sn0*AIn!!!!!

    M1 = size(MIn, 1);

    Aeq = Dg.' * En_R * Mmat - D0g.' * M0;
    rc = rcond(Aeq);
    if rc < 1e-10
        warning('scatteringMatrices:illConditioned', ...
            'Aeq nearly singular (rcond=%.2e) -- results may be unreliable.', rc);
    end

    MU = Aeq \ (D0g.' * (eye(M1) + MIn));
    Sn = Mmat * MU;
    Sn0 = MIn + M0 * MU;
end