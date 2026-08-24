%% Outputs the corrected radial matrix (used in 2.24 and derived in A.5)
function En = radialMatrix(n,r,R,kappa,V,Ln)
    
    % in keeping with Montiel's notation, I will use:
    % En is the corrected result
    % Escript_n is the diag interior amps matrix m>0
    % EscriptTilde_n are for the damped modes m = {-2, -1}

    % Thus Escript_n(r) is (M+1) x (M+1)
    % with diag elements Ihat_n(kappa_m*r) ofr m = 0, ...,M

    % EscriptTilde(r) is 2 x 2
    % with diag elements Ihat_n(kappa_{-j}*r) corr to
    %[kappa_{-2}, kappa_{-1}] to match V's column order
    %kappa is obv (M+3)x1

    % V is the (M+1) x 2 matrix from bennettsCoefficients.m

    % Ln is the 2 x (M+1) from edgeMap.m

    kappaTower = kappa(3:end);
    M1 = length(kappaTower);
    kappaDamped = [kappa(1), kappa(2)]; %ie kappa_{-2} and {-1}
    
    EscriptDiag = zeros(M1,1);
    for m = 1:M1
        EscriptDiag(m) = scaledBesselI(n, kappaTower(m), r, R);
    end
    Escript = diag(EscriptDiag);

    EscriptTildeDiag = zeros(2,1);
    for j = 1:2
        EscriptTildeDiag(j) = scaledBesselI(n, kappaDamped(j), r, R);
    end
    EscriptTilde = diag(EscriptTildeDiag);

    En = Escript + V*EscriptTilde*Ln;


end

function val = scaledBesselI(n, K, r, R)
    %Ihat_n(Kr) = I_n(Kr)/I_n(KR) -> eq. (6). Same overflow caveat as
    %besselILogDerivative in edgeMap.m
    val = besseli(n, K*r) / besseli(n, K*R);
end