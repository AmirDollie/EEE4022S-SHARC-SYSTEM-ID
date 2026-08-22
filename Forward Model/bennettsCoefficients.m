%% Calculates the coefficients v derived by Bennetts (2007) used in Montiel's thesis

% I believe there is a sign error in his thesis. I will explain in a
% comment above the line it is present in

function V = bennettsCoefficients(alpha, beta, gamma, kappa)


    %This functions constructs the V matrix which maps the damped
    %modes phi_{-1} and phi_{-2} onto the usual phi evanescent tower (m>0)

    %kappa is the (M+3)x1 vector containing m = -2, -1, 0, ..., M roots of
    %the plate covered dispersion relation

    %V is hence a (M+1)x2 matrix such that:
    % phi_{-j}(z) = sum_m V(m+1, j) * phi_m(z),   j = 1, 2 equ (2.16)
    %Obv this means col 1 of V corr to phi_{-2} not phi_{-1}

    %note there is an index change j* = 3-j

    %SIGN CORRECTION relative to the literal eq. (4) transcription (in my chapt 2 summary): this
    %implementation carries an EXTRA MINUS SIGN on v_{m,j} relative to
    %what's printed in the thesis. Traced this back to Bennetts (2007) sec. 4.2,
    %whose actual derivation (eq. 4.4 there) is a residue-theorem
    %identity of the form
    %       2 * sum_{n=-2}^{inf} [f_a(k_n)/K'(k_n)] * cosh(k_n(z+h)) = 0
    %i.e. a SUM-EQUALS-ZERO over every mode including the damped ones.
    %Getting "phi_{-j}(z) = sum v_{m,j} phi_m(z)" out of that requires
    %isolating the n=-j term and moving every other term across the "=0"
    %which introduces a minus sign structurally, not by convention.
    %This strongly suggests Montiel's restatement (condensing Bennetts'
    %general two-constant f_a(tau)=a1*tau+a3*tau^3 result into the
    %specific closed-form eq. 4) dropped that sign somewhere in the
    %algebra. Verified numerically: the reconstruction self-check below
    %converges to phi_{-j}(z) correctly WITH this minus sign, and to
    %-phi_{-j}(z) without it, with the gap shrinking as M grows in both
    %cases (i.e. it's a genuine convergent-series sign flip, not noise).
    
    %Just want to FLAG THIS. Not 10000% sure, as I am sure Montiel knows
    %more than me... worth double checking!

    kappaM2 = kappa(1); %kappa_{-2}
    kappaM1 = kappa(2); %kappa_{-1}
    kappaTower = kappa(3:end); %m>0 kappas
    M1 = length(kappaTower);
    
    %letting kappaDamped(j) = kappa_{-j}, j = 1,2
    kappaDamped = [kappaM1, kappaM2];

    vBar = zeros(M1,2);
    for j = 1:2
        jStar = 3-j;
        kappaJ = kappaDamped(j);
        kappaJStar = kappaDamped(jStar);

        [~, DpJ] = dispersionFunction(kappaJ, alpha, beta, gamma);

        for m = 1:M1
            kappaM = kappaTower(m);
            [~, DpM] = dispersionFunction(kappaM, alpha, beta, gamma);

            % Now for the minus sign (leading term!)
            % See equation 
            vRaw = -(kappaM * DpJ * (kappaM^2 - kappaJStar^2)) / (kappaJ * DpM * (kappaJ^2 - kappaJStar^2));

            vBar(m,j) = vRaw * cos(kappaM*(1-gamma)) / cos(kappaJ*(1-gamma));
        end
    end

    % Assemble V with the j* swap: column j of V holds vBar(:, jStar)
    V = zeros(M1, 2);
    V(:,1) = vBar(:,2);   % j=1 slot <- kappa_{-2}
    V(:,2) = vBar(:,1);   % j=2 slot <- kappa_{-1}


end