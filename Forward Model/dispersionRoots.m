% This function will solve the 2 dispersion relations (2.14) and (2.15) in
% their relevant domains (open water and disc covered)

function xi = dispersionRoots(alpha, beta, gamma, M)

% Note that the disc covered dispersion relation collapses to the open water
%one when beta = gamma = 0. 

% Hence pass BETA = 0 = GAMMA for open water 
% returning [k0, k1, ..., kM]' ie (M+1)x1
% with k0 being purely imaginary and km, m>0 ascending (evanescent)

% Disc covered relation:
% returning [kappa{-2}, kappa{-1}, kappa0, ...,kappaM]' ie (M+3)x1
% includes the complex conj damped pair ahead of travelling root kappa0

% If you multiply the dispersion relations by cos(K), you remove the poles
% introduced by the tan(K) term, so for open water:
% Ktan(K) = -alpha becomes Ksin(K) = -alpha*cos(K)
% aka D(K) = Ksin(K) + alpha*cos(K) whose roots are km

% And disc covered:
% D(K) = (1 - alpha*gamma + beta*K^4)*K*sin(K*(1-gamma)) + alpha*cos(K*(1-gamma)
% whose roots are kappa_p
% and again as a sanity check u can see is the same as open with
% gamma and beta = 0!

%Note that for now, this works for REAL beta, not viscoelastic yet!

% D(xi) and D'(xi) now live in their own file, dispersionFunction.m,
% since bennettsCoefficients.m (Step 1) also needs D' and duplicating
% the formula in two places would break the "only one place to fix the
% sign discrepancy" guarantee. Call dispersionFunction(xi,alpha,beta,gamma)
% wherever entireFn used to be called locally.
    
    %Binary variable to differentiate open or disc region:
    isOpenWater = (beta == 0) && (gamma == 0);
    
    kappa0 = solveTravellingRoot(alpha, beta, gamma);
    kappaEvanescent = solveEvanescentRoots(alpha, beta, gamma, M);

    if isOpenWater
        xi = [kappa0; kappaEvanescent];
    else
        [kappaM2, kappaM1] = solveDampedPair(alpha, beta, gamma);
        xi = [kappaM2; kappaM1; kappa0; kappaEvanescent];
    end
end

% ============ Now for my function definitions ===========

% ---------- travelling(m = 0) root finder -----------
function xi0 = solveTravellingRoot(alpha, beta, gamma)
    %note xi0 = i*p0, found by zeroing the Re{D(ip)} over p>0
    %This follows Sturm-Liouville theory:
    %exactly 1 root is on the positive imag axis!

    %Messed up if we consider viscoelasticity!!!!
    if ~isreal(beta)
        warning('dispersionRoots:complexBeta', ...
            ['beta is complex: this assumes the travelling root stays ' ...
             'on the imaginary axis, which is only verified for real ' ...
             'beta. Treat the result as a starting guess only.']);
    end
    
    residual = @(p) real(dispersionFunction(1i*p, alpha, beta, gamma));
    %Arb upper bound
    pMax = 10*alpha + 20;
    %scan p for near 0 (1e-8) to upper bound
    brackets = findSignChangeBrackets(residual, 1e-8, pMax, 2000);

    if isempty(brackets)
        error('dispersionRoots:noTravellingRoot', ...
              'No sign change found for the travelling root on (0, %.4g).', pMax);
    end
    %modified Matlab's fzero function for sign-change interval
    p0 = bracketAndBisect(residual, brackets(1,1), brackets(1,2));
    xi0 = 1i * p0;
end

% ---------- Evanescent tower ----------

function xi = solveEvanescentRoots(alpha, beta, gamma, M)

    % find the M real, ascending roots xi(m) ~ m*pi/(1-gamma) for big m)
    residual = @(K) dispersionFunction(K, alpha, beta, gamma);
    xMax = (M + 0.5) * pi/(1 - gamma);
    brackets = findSignChangeBrackets(residual, 1e-8, xMax, 4000*(M+1));

    if size(brackets, 1) < M
        error('dispersionRoots:tooFewEvanescentRoots', ...
              'Found only %d sign changes on (0, %.4g), need %d.', ...
              size(brackets,1), xMax, M);
    end
    
    xi = zeros(M,1);
    for m = 1:M
        xi(m) = bracketAndBisect(residual, brackets(m,1), brackets(m,2));
    end
end

% ---------- Damped Pair (requires newtonComplex) -----------

%This is the demon function. Caused me much pain today wrt convergence...
%which is why you will see there are many catch error sections.
%Needs proper validation.
function [xiM2, xiM1] = solveDampedPair(alpha, beta, gamma, seed)
    residual = @(xi) dispersionFunction(xi, alpha, beta, gamma);

    if nargin < 4 || isempty(seed)
        % Dense grid over the first quadrant: for small beta the true
        % damped root's basin of attraction is narrow (real-axis roots
        % dominate), so a handful of hand-picked seeds isn't reliable --
        % verified this empirically across several (alpha,beta,gamma).
        reRange = linspace(0.5, 30, 15);
        imRange = linspace(0.5, 20, 15);
        candidateSeeds = reRange(:) + 1i*imRange(:).';
        candidateSeeds = candidateSeeds(:);
    else
        candidateSeeds = seed;
    end

    xiM2 = [];
    bestNorm = Inf;
    for k = 1:numel(candidateSeeds)
        [root, converged] = newtonComplex(residual, candidateSeeds(k));
        if converged && abs(imag(root)) > 1e-6 && real(root) > 0.3
            if abs(residual(root)) < 1e-6 && abs(root) < bestNorm
                xiM2 = root;
                bestNorm = abs(root);
            end
        end
    end

    if isempty(xiM2)
        error('dispersionRoots:noDampedRoot', ...
              'No seed in the grid converged to a valid damped root.');
    end

    if isreal(beta)
        xiM1 = conj(xiM2);
    else
        warning('dispersionRoots:complexBeta', ...
            'beta is complex: conjugate symmetry does not hold; xiM1 is a placeholder.');
        xiM1 = conj(xiM2);
    end
end
% ---------- bracketAndBisect -----------

function root = bracketAndBisect(fun, xLow, xHigh)
    fLow = fun(xLow); fHigh = fun(xHigh);
    %checks if one is + and other - thus there exists zero crossing between
    if sign(fLow) == sign(fHigh) && fLow ~= 0 && fHigh ~= 0
        error('dispersionRoots:badBracket', ...
              'fun(%.6g)=%.6g and fun(%.6g)=%.6g do not bracket a root.', ...
              xLow, fLow, xHigh, fHigh);
    end
    %if they have opposite signs uses matlabs fzero cause im too lazy
    root = fzero(fun, [xLow, xHigh]);
end

% ---------- findSignChangeBrackets ----------
function brackets = findSignChangeBrackets(fun, xMin, xMax, nSamples)
    %here provides the actual finding of the brackets 
    %that the above wrapper acts on
    %pretty obvious how this works 
    %This was the reason i got rid of the poles so a blunt search wont
    %return an error!!
    x = linspace(xMin, xMax, nSamples);
    f = arrayfun(fun, x);
    idx = find(sign(f(1:end-1)) ~= sign(f(2:end)));
    brackets = [x(idx)', x(idx+1)'];
end

% ---------- newtonComplex -----------
%Literally Newton's method but on complex plane
%recall OG formula is x_new = x_old - f(x_old)/f'(x_old)
%here x is complex
%I use central difference as no closed form derivative is supplied
%Can make better later, maybe a genetic algorithm?
%This is currently my biggest form of error that I can think of
function [root, converged] = newtonComplex(fun, x0)
    tol = 1e-10; maxIter = 50; step = 1e-6;
    x = x0; converged = false;
    for iter = 1:maxIter
        fx = fun(x);
        if abs(fx) < tol
            converged = true;
            break
        end
        dfx = (fun(x+step) - fun(x-step)) / (2*step);
        if dfx == 0
            error('dispersionRoots:zeroDerivative', 'Zero derivative at iteration %d.', iter);
        end
        x = x - fx/dfx;
    end
    root = x;
end