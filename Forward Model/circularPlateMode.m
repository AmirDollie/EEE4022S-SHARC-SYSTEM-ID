function mode = circularPlateMode(n, j, R, nu, varargin)
%CIRCULARPLATEMODE  Dry, free-edge circular-plate eigenmode (in-vacuo basis).
%
%   mode = circularPlateMode(n, j, R, nu) returns the j-th (0-indexed)
%   radial eigenmode of angular order n for a FREE-EDGE circular thin
%   plate of radius R and Poisson ratio nu -- the dry biharmonic plate
%   problem
%
%       nabla^4 W - lambda^4/R^4 * W = 0,   lambda = k*R (nondimensional)
%
%   separated as W(r,theta) = R_{n,j}(r) * {cos(n*theta) or sin(n*theta)},
%   with the radial part built from ordinary + modified Bessel functions
%   (NOT the modified-only I_n form used in the EMM's disc-covered
%   region -- the dry plate has no evanescent-only structure, it is a
%   genuinely oscillatory eigenproblem):
%
%       R_{n,j}(r) = J_n(lambda*r/R) + C*I_n(lambda*r/R)
%
%   This file solves ONLY the dry, in-vacuo plate eigenproblem. There is
%   no water, no EMM, and no fluid coupling anywhere in this function.
%
%   *** IMPORTANT -- three concepts that must not be conflated ***
%   (1) a "dry plate eigenfrequency" (a lambda_{n,j} returned here) is
%       NOT (2) a "wet hydroelastic modal-response peak" (a resonance of
%       the projected modal RAO |A_{n,j}(omega)| once the EMM/fluid is
%       included), which is in turn NOT (3) an "SSI-identified pole"
%       (a resonance recovered from a simulated/measured time series by
%       stochastic subspace identification). They may end up close, but
%       nothing here establishes that, and the inverse-problem cost
%       function must not assume it.
%
%   RIGID-BODY MODES are not roots of the free-edge determinant below --
%   the boundary-condition equations degenerate identically at lambda=0.
%   They are constructed directly instead:
%       n=0, j=0:  R(r) = 1        (heave,        lambda = 0)
%       n=1, j=0:  R(r) = r        (tilt,         lambda = 0)
%   (n=1 has two independent rigid tilt directions -- cos(theta) and
%   sin(theta) -- but they share this one radial part/eigenvalue, so
%   only one j=0 entry is returned here; the angular partner is chosen
%   by the caller exactly as for the flexible n>0 modes.)
%
%   All other (n,j) pairs are the (0-indexed) j-th positive root of
%
%       det B_n(lambda) = 0
%
%   built from imposing the two free-edge conditions (zero radial
%   bending moment, zero Kirchhoff effective shear) on R_{n,j}(r) at
%   r=R. Both conditions are the standard closed-form circular-plate
%   free-edge conditions (see e.g. Leissa, "Vibration of Plates", NASA
%   SP-160); they were re-derived from the plate biharmonic equation for
%   this file rather than transcribed from Itao & Crandall (1979),
%   because Meylan & Squire (1996) explicitly flag that paper's
%   published normalisation constant as erroneous (corrected in Meylan
%   1995) -- rather than chase that correction down, this file sidesteps
%   the whole question: the eigenvalue equation used here does not
%   depend on that normalisation at all, and the mode is normalised
%   numerically (see below), independently of any literature constant.
%
% Inputs:
%   n   angular (nodal-diameter) order, integer >= 0
%   j   radial (nodal-circle) index, integer >= 0, 0-indexed. j=0 is the
%       lowest eigenvalue at that n -- for n=0 or n=1 this IS the
%       rigid-body mode; for n>=2 it is the lowest flexible mode.
%   R   plate radius
%   nu  Poisson ratio
%
% Name-value options:
%   'LambdaMax'   upper search limit for lambda            (default 40)
%   'NBrackets'   number of scan points used to bracket     (default 4000)
%                 roots of det B_n(lambda) over (0, LambdaMax]
%
% Output (struct mode):
%   .n, .j          as given
%   .lambda         nondimensional eigenvalue (k*R); 0 for rigid-body modes
%   .C              coefficient of the I_n branch (NaN for rigid-body modes)
%   .isRigidBody    true for the n=0,j=0 heave and n=1,j=0 tilt modes
%   .radial         function handle R_{n,j}(r), r in [0,R], numerically
%                   normalised so that
%                       int_0^{2pi} int_0^R |w_{n,j}(r,theta)|^2 r dr dtheta = 1
%                   with w_{n,j} = R_{n,j}(r)*cos(n*theta) (n=0 uses the
%                   same cos(0*theta)=1 convention). This is a direct
%                   numerical quadrature of the definition, not a
%                   transcription of any published closed-form constant.
%   .angularKind    always 'cos' here; the sine partner for n>0 shares
%                   this radial part and lambda exactly -- multiply by
%                   sin(n*theta) instead of cos(n*theta) to get it. It is
%                   not returned separately to avoid silently duplicating
%                   (lambda, C, radial shape) data that is identical
%                   between the two.
%
% See circularPlateModeTester.m for the validation this basis needs
% before anything gets projected onto it: free-edge residual checks,
% numerical orthogonality between distinct j at fixed n, rigid-body
% modes identified separately from the root search, and eigenvalue
% convergence with root-search resolution.

    p = inputParser;
    p.addParameter('LambdaMax', 40);
    p.addParameter('NBrackets', 4000);
    p.parse(varargin{:});
    lambdaMax = p.Results.LambdaMax;
    nBrackets = p.Results.NBrackets;

    if n < 0 || j < 0 || mod(n,1) ~= 0 || mod(j,1) ~= 0
        error('circularPlateMode:badIndex', 'n and j must be non-negative integers.');
    end

    mode.n = n;
    mode.j = j;
    mode.angularKind = 'cos';

    isRigidCandidate = (n == 0 && j == 0) || (n == 1 && j == 0);

    if isRigidCandidate
        mode.lambda = 0;
        mode.C = NaN;
        mode.isRigidBody = true;
        if n == 0
            rawRadial = @(r) ones(size(r));
        else % n == 1
            rawRadial = @(r) r;
        end
    else
        mode.isRigidBody = false;

        % Flexible-root index within the det(B_n)=0 root ladder: for
        % n=0/1, j=0 was consumed by the rigid-body mode above, so the
        % j-th requested mode is the (j-1)-th flexible root; for n>=2
        % there is no rigid-body mode at this n, so j itself indexes the
        % flexible root directly.
        if n <= 1
            flexIndex = j - 1; % 0-indexed among flexible roots
        else
            flexIndex = j;
        end

        [lambda, C] = findFreeEdgeRoot(n, nu, flexIndex, lambdaMax, nBrackets);
        mode.lambda = lambda;
        mode.C = C;
        rawRadial = @(r) besselj(n, lambda*r/R) + C*besseli(n, lambda*r/R);
    end

    % --- numerical normalisation (bypasses the historical Itao/Crandall
    % formula entirely -- see docstring) ---
    if n == 0
        angularFactor = 2*pi; % int_0^2pi cos(0)^2 dtheta
    else
        angularFactor = pi;   % int_0^2pi cos(n theta)^2 dtheta
    end
    radialIntegral = integral(@(r) rawRadial(r).^2 .* r, 0, R);
    A = 1/sqrt(angularFactor * radialIntegral);

    mode.radial = @(r) A*rawRadial(r);
end

% =====================================================================
function [lambda, C] = findFreeEdgeRoot(n, nu, flexIndex, lambdaMax, nBrackets)
% Finds the flexIndex-th (0-indexed) positive root of det B_n(lambda)=0,
% where B_n is the 2x2 matrix of [moment; shear] free-edge conditions
% evaluated on the [J_n branch, I_n branch] basis, then recovers the
% branch coefficient C from the (necessarily consistent, since det=0)
% homogeneous system.

    if flexIndex < 0
        error('circularPlateMode:badIndex', ...
            'j=0 for n=0 or n=1 is the rigid-body mode; there is no flexible root at that index.');
    end

    detFun = @(lam) detB(lam, n, nu);

    % scan for sign changes, starting just above 0 to avoid the
    % removable-singularity point of the BC expressions at lambda=0
    lamGrid = linspace(1e-3, lambdaMax, nBrackets);
    dvals = arrayfun(detFun, lamGrid);

    rootsFound = [];
    for k = 1:numel(lamGrid)-1
        a = lamGrid(k); b = lamGrid(k+1);
        fa = dvals(k); fb = dvals(k+1);
        if isfinite(fa) && isfinite(fb) && sign(fa) ~= sign(fb) && fa ~= 0
            r = bisectRoot(detFun, a, b);
            rootsFound(end+1) = r; %#ok<AGROW>
        end
    end

    if numel(rootsFound) <= flexIndex
        error('circularPlateMode:notEnoughRoots', ...
            ['Only found %d flexible root(s) for n=%d below LambdaMax=%.1f; ', ...
             'requested flexible index %d. Increase LambdaMax or NBrackets.'], ...
            numel(rootsFound), n, lambdaMax, flexIndex);
    end

    rootsFound = sort(rootsFound);
    lambda = rootsFound(flexIndex+1);

    % Newton polish: bisection is only guaranteed to 1e-13 in lambda, and
    % at higher (n,lambda) I_n(lambda) grows exponentially while C shrinks
    % to compensate, so the residual is far more sensitive to lambda error
    % there than at low modes. A few Newton steps on det B_n (central
    % finite-difference derivative) squeeze this further at negligible
    % cost. With the corrected BC2/Bessel-derivative formulas (see the
    % history note in besselJderivs below), this now leaves modal
    % orthogonality at machine precision at every (n,j) tested, not just
    % approximately small -- an earlier version of this file, with two
    % separate algebra bugs, could not get better than ~1e-3 to ~1e-5
    % however much this polish step was tightened, which was the tell
    % that the bugs were in the formulas, not the root-finding precision.
    h = 1e-6;
    for iter = 1:6
        d0 = detFun(lambda);
        dprime = (detFun(lambda+h) - detFun(lambda-h)) / (2*h);
        if dprime == 0 || ~isfinite(dprime)
            break;
        end
        step = d0/dprime;
        lambda = lambda - step;
        if abs(step) < 1e-14
            break;
        end
    end

    % Recover C from the moment-condition row (BC1), which is exact at a
    % true root of det B_n since the 2x2 system is then singular and the
    % same [1;C] also satisfies the shear condition (BC2) -- this is
    % checked explicitly in circularPlateModeTester.m, not just assumed.
    [bc1J, bc1I, ~, ~] = edgeConditionValues(lambda, n, nu);
    C = -bc1J/bc1I;
end

function r = bisectRoot(f, a, b)
    fa = f(a); fb = f(b);
    for iter = 1:80
        m = 0.5*(a+b);
        fm = f(m);
        if fm == 0 || (b-a) < 1e-13
            r = m; return;
        end
        if sign(fm) == sign(fa)
            a = m; fa = fm;
        else
            b = m; fb = fm;
        end
    end
    r = 0.5*(a+b);
end

function d = detB(lambda, n, nu)
    [bc1J, bc1I, bc2J, bc2I] = edgeConditionValues(lambda, n, nu);
    d = bc1J*bc2I - bc1I*bc2J;
end

function [bc1J, bc1I, bc2J, bc2I] = edgeConditionValues(lambda, n, nu)
% Nondimensional free-edge boundary condition values (moment = BC1,
% Kirchhoff effective shear = BC2) at x = lambda, for each of the two
% branches g = J_n and g = I_n, re-derived directly from
%
%   M_r/D  = R'' + nu*( R'/r - n^2 R/r^2 )                        = 0
%   V_r/D  = R''' + R''/r - (1/r^2)*[1+n^2*(2-nu)]*R'
%                  + (1/r^3)*n^2*(3-nu)*R                          = 0
%
% (standard Kirchhoff free-edge conditions for a circular plate, mode
% cos(n*theta)/sin(n*theta)) nondimensionalised with x = k*r, lambda = k*R:
%
%   BC1(lambda) = lambda^2*g''(lambda) + nu*lambda*g'(lambda) - nu*n^2*g(lambda)
%   BC2(lambda) = lambda^3*g'''(lambda) + lambda^2*g''(lambda)
%                 - lambda*[1+n^2*(2-nu)]*g'(lambda) + n^2*(3-nu)*g(lambda)
%
% NOTE: an earlier version of this file dropped the R''/r term above when
% assembling BC2 -- it produced roots whose BC1/BC2 residuals looked
% machine-precision-consistent (because det B_n=0 was still solved
% self-consistently for the WRONG B_n), but which were not true
% eigenvalues: n=0 cross-mode orthogonality came out at the ~5e-3 level
% instead of machine precision. Caught only by actually running the
% orthogonality tester below, not by the residual check alone -- residual
% self-consistency of a root-finder is not the same as testing that the
% equation it solved is the correct one. Left here as a flag for anyone
% re-deriving these formulas: recompute d/dr(nabla^2 W) in full and check
% every term survives into the final nondimensional form.

    [Jv, Jp, Jpp, Jppp] = besselJderivs(n, lambda);
    [Iv, Ip, Ipp, Ippp] = besselIderivs(n, lambda);

    bc1J = lambda^2*Jpp + nu*lambda*Jp - nu*n^2*Jv;
    bc1I = lambda^2*Ipp + nu*lambda*Ip - nu*n^2*Iv;

    bc2J = lambda^3*Jppp + lambda^2*Jpp - lambda*(1+n^2*(2-nu))*Jp + n^2*(3-nu)*Jv;
    bc2I = lambda^3*Ippp + lambda^2*Ipp - lambda*(1+n^2*(2-nu))*Ip + n^2*(3-nu)*Iv;
end

function [Jv, Jp, Jpp, Jppp] = besselJderivs(n, x)
% J_n and its first three derivatives at x, via the recurrence relation
% for J_n' and the Bessel ODE for J_n'', differentiated once more (in
% closed form, using the ODE again) for J_n'''. MATLAB's besselj handles
% the negative-order case (n-1 = -1 when n=0) natively via
% J_{-m} = (-1)^m J_m, so no special-casing is needed for n=0.
%
% NOTE: an earlier version of this file had the WRONG coefficient of Jv
% here (+n^2/x^3 instead of -3*n^2/x^3), inherited from a hand-algebra
% slip in re-differentiating the Bessel ODE. It agreed with the true
% third derivative at n=0 (where the n^2 term vanishes and the bug is
% invisible) but was wrong for every n>=1 -- exactly why the n=0 basis
% validated perfectly while n=1 (rigid tilt vs. its own flexible modes)
% showed a real, non-numerical ~5e-3 orthogonality violation. Re-derived
% and confirmed against sympy's built-in besselj third derivative to
% machine precision (see project scratch notes) before fixing:
%   J''' = J'*[(n^2+2)/x^2 - 1] + J*[1/x - 3*n^2/x^3]
% This is the second of two independent algebra errors this file has had
% (the first was a dropped R''/r term in the plate shear BC). Both were
% invisible to the self-consistency residual check alone -- that check
% only verifies a root satisfies MY OWN BC formulas, not that those
% formulas are the physically correct ones. Only the orthogonality test
% catches an internally-self-consistent-but-wrong derivation. Do not
% trust a residual-only check on this file again.
    Jv  = besselj(n, x);
    Jp  = besselj(n-1, x) - (n/x)*Jv;
    Jpp = -Jp/x - (1 - n^2/x^2)*Jv;
    Jppp = Jp*((n^2+2)/x^2 - 1) + Jv*(1/x - 3*n^2/x^3);
end

function [Iv, Ip, Ipp, Ippp] = besselIderivs(n, x)
% I_n and its first three derivatives at x, same approach as
% besselJderivs but for the modified Bessel ODE
% (x^2 I'' + x I' - (x^2+n^2) I = 0). besseli handles negative integer
% order natively via I_{-m} = I_m.
%   I''' = I'*[(n^2+2)/x^2 + 1] + I*[-1/x - 3*n^2/x^3]
% (see besselJderivs for the history of why this needed re-deriving and
% sympy-checking rather than trusting the original hand algebra).
    Iv  = besseli(n, x);
    Ip  = besseli(n-1, x) - (n/x)*Iv;
    Ipp = -Ip/x + (1 + n^2/x^2)*Iv;
    Ippp = Ip*((n^2+2)/x^2 + 1) + Iv*(-1/x - 3*n^2/x^3);
end