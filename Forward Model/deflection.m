%% This is the plate deflection output given a specific incident wave
function eta = deflection(r, theta, alpha, beta, gamma, R, nu, M, P, N)

% Physical plate deflection (equ A.2)
% along with with phi'_m(0) corrected to phi'_m(-gamma) 
% Ties together the entire pipeline:

% roots -> Bennetts -> edge map -> matching -> scattering -> forcing 
% -> deflection sum

% usage: %   ETA = DEFLECTION(R, THETA, ALPHA, BETA, GAMMA, R_PLATE, NU, M, P, N)

% A lot of things don't depend on the angular order n
% So i will compute them outside of the loop:

    kappa = dispersionRoots(alpha, beta, gamma, M);
    k = dispersionRoots(alpha, 0, 0, M);
    V = bennettsCoefficients(alpha, beta, gamma, kappa);
    k0 = k(1);
    
    D0 = innerProductMatrix(k, 0);
    D = innerProductMatrix(kappa(3:end), gamma);
    D0g = gegenbauerCrossMatrix(k, gamma, P, 1);
    Dg = gegenbauerCrossMatrix(kappa(3:end), gamma, P, 1-gamma);

    kappaFull = kappa;   % [kappa_{-2}; kappa_{-1}; kappa0; ...; kappaM]


    etaSum = 0;
    for n = -N:N
        %compute edge map, radial derivatives, matching, scattering
        Ln = edgeMap(n, alpha, beta, gamma, R, nu, kappa);
        [EIprime, E0prime] = openWaterRadialDerivatives(n, R, k);
        Enprime = radialMatrixDerivative(n, R, kappa, V, Ln);
        [MIn, M0, Mmat] = matchingMatrices(EIprime, E0prime, Enprime, D0, D, D0g, Dg);
        En_R = radialMatrix(n, R, R, kappa, V, Ln);
        [~, Sn, ~] = scatteringMatrices(D0g, Dg, En_R, MIn, M0, Mmat);
        
        %forward propogation and actual physical forcing:
        AIn = incidentAmplitude(alpha, n, k0, R, M);
        An = Sn*AIn;
        Aen = Ln*An;
        AfullVec = [Aen(1); Aen(2); An];

        innerSum = 0;
        for idx = 1:length(kappaFull)
            km = kappaFull(idx);
            dphi = edgeDerivative(km, gamma);
            IhatVal = besseli(n, km*r) / besseli(n, km*R);
            innerSum = innerSum + dphi * AfullVec(idx) * IhatVal;
        end
        %equ A.2
        etaSum = etaSum + innerSum * exp(1i*n*theta);
    end

    eta = etaSum / (1i*alpha);        

end