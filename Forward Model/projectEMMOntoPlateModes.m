function A = projectEMMOntoPlateModes(wFun, modeList, R, nu, varargin)
%PROJECTEMMONTOPLATEMODES  Project a deflection field onto the dry,
%free-edge circular-plate basis (circularPlateMode.m).
%
%   A = PROJECTEMMONTOPLATEMODES(wFun, modeList, R, nu) returns the
%   modal amplitude coefficients
%
%       A_{n,j} = <w, w_{n,j}> = int_0^{2pi} int_0^R w(r,theta) *
%                                 R_{n,j}(r) * cos(n*theta) * r dr dtheta
%
%   for each (n,j) pair in modeList, where w_{n,j}(r,theta) =
%   R_{n,j}(r)*cos(n*theta) is the ORTHONORMAL dry mode returned by
%   circularPlateMode(n,j,R,nu) (orthonormality verified to ~1e-16 in
%   circularPlateModeTester.m -- this function's correctness rests on
%   that).
%
%   This is Step 2/4 of the dry-mode inverse-mapping plan
%   (circularPlateMode.m -> circularPlateModeTester.m -> THIS FUNCTION
%   -> modalRAOReconnaissance.m): it turns a 2D deflection field into a
%   finite set of modal amplitudes, w_EMM(r,theta;alpha) -> A_{n,j}(alpha).
%
% Inputs:
%   wFun      function handle, wFun(r, theta) -> complex deflection.
%             MUST be ELEMENTWISE VECTORIZED: given r and theta as
%             arrays of the SAME size (any shape), it must return an
%             array of that same size, evaluating the field at each
%             (r(k), theta(k)) pair independently -- i.e. the ordinary
%             MATLAB broadcasting convention (besselj, exp, etc. all
%             work like this natively), NOT a meshgrid/outer-product
%             convention. r=0 must be handled (the plate centre).
%
%             This is deliberately the same calling convention as this
%             project's own Forward Model/deflection.m(r, theta, alpha,
%             beta, gamma, R, nu, M, P, N), which is itself elementwise
%             over matching-size r/theta arrays (each term is built from
%             besseli(n, km*r) and exp(1i*n*theta), both elementwise).
%             So wiring in the real EMM output later is expected to be
%             no more than
%                 wFun = @(r,theta) deflection(r, theta, alpha, beta, ...
%                     gamma, R, nu, M, P, N);
%             This has NOT been done or tested here -- deflection.m has
%             its own long dependency chain (dispersionRoots ->
%             bennettsCoefficients -> edgeMap -> matchingMatrices ->
%             scatteringMatrices -> ...) that is out of scope for this
%             file and hasn't been independently validated in this
%             session. See projectEMMOntoPlateModesTester.m, which
%             validates the PROJECTION MATH ONLY, using synthetic
%             fields, not a real EMM field.
%
%   modeList  Kx2 matrix of [n j] pairs to project onto, one row per
%             requested mode (matches circularPlateMode's (n,j)
%             convention exactly, including the rigid-body j=0 entries
%             at n=0,1). Duplicate rows are allowed but wasteful.
%
%   R, nu     plate radius, Poisson ratio -- passed straight through to
%             circularPlateMode, so must match whatever the field wFun
%             was itself computed for.
%
% Name-value options:
%   'NTheta'       number of angular quadrature points used to collapse
%                  the theta integral at each n (default 720, i.e. 0.5
%                  degree spacing). Uses the PERIODIC TRAPEZOIDAL RULE
%                  (equivalently, a uniform-weight rectangle sum) over a
%                  full period -- for a smooth 2*pi-periodic integrand
%                  (which cos(n*theta)*w(r,theta) always is, for
%                  integer n and any smooth w) this converges
%                  spectrally (faster than any fixed polynomial order),
%                  so 720 points is already far more than needed for
%                  the low n this project uses; it is not the accuracy
%                  bottleneck.
%   'AngularKind'  'cos' (default) or 'sin' -- which angular partner to
%                  project onto. Montiel's own analysis only needs
%                  'cos': for a single disc forced by a straight-crested
%                  wave, the response is symmetric about theta=0,pi (no
%                  roll component -- see Montiel_Thesis_6.pdf Sec 6.1),
%                  so only n>=0 with the cos(n*theta) family is
%                  physically populated. 'sin' is offered for
%                  completeness/future use but is NOT exercised by
%                  anything in this project yet, and errors for n=0
%                  (sin(0*theta) is identically zero, undefined to
%                  project onto). Mixed cos+sin projection in a single
%                  call is not offered -- call this function twice (with
%                  different modeList/AngularKind) if that's ever needed.
%
% Output:
%   A   Kx1 complex column vector, A(k) is the amplitude for
%       modeList(k,:). If wFun truly lies in the span of the requested
%       modes, reconstructing
%           w_recon(r,theta) = sum_k A(k) * R_{n_k,j_k}(r) * cos(n_k*theta)
%       reproduces wFun exactly (to quadrature precision); otherwise
%       w_recon is the best least-squares approximation of wFun using
%       exactly those modes (a standard orthogonal-projection property,
%       since the basis is orthonormal -- see
%       projectEMMOntoPlateModesTester.m Test 2 for the convergence
%       check this implies: reconstruction error can only decrease, or
%       stay the same, as more modes from a nested sequence are added).
%
% Method: a naive implementation would run a full 2D quadrature
% (integral2) once per requested mode. Instead, this exploits the
% angular-Fourier orthogonality of cos(n*theta) to collapse the theta
% integral ONCE PER DISTINCT n present in modeList (not once per mode),
% producing a radial function
%
%   C_n(r) = int_0^{2pi} w(r,theta) * cos(n*theta) dtheta
%
% which every j at that n then only needs a cheap 1D radial integral
% against:
%
%   A_{n,j} = int_0^R R_{n,j}(r) * C_n(r) * r dr
%
% This is the reduction described in the original 4-step plan for this
% function. Note this reuses the SAME construction (same theta grid, one
% wFun evaluation pattern) across every j at a given n -- it does not
% literally cache numeric C_n(r) values across calls, since the radial
% quadrature (integral, adaptive) may sample different r nodes for
% different j; the saving is architectural (1D collapsed integrand
% instead of a fresh 2D quadrature per mode), not a literal value cache.

    p = inputParser;
    p.addParameter('NTheta', 720);
    p.addParameter('AngularKind', 'cos');
    p.parse(varargin{:});
    nTheta = p.Results.NTheta;
    angularKind = lower(p.Results.AngularKind);

    if ~ismember(angularKind, {'cos', 'sin'})
        error('projectEMMOntoPlateModes:badAngularKind', ...
            "AngularKind must be 'cos' or 'sin', got '%s'.", angularKind);
    end
    if any(modeList(:,1) == 0) && strcmp(angularKind, 'sin')
        error('projectEMMOntoPlateModes:zeroSin', ...
            'sin(0*theta) is identically zero -- AngularKind=''sin'' is undefined/meaningless for n=0.');
    end

    nModes = size(modeList, 1);
    A = zeros(nModes, 1);

    % Periodic trapezoidal rule = uniform rectangle rule over one full
    % period; deliberately NOT including theta=2*pi as a separate node
    % (that would double-count the periodic endpoint against theta=0).
    thetaGrid = (0:nTheta-1) * (2*pi/nTheta);
    dtheta = 2*pi/nTheta;

    uniqueN = unique(modeList(:,1));
    for kn = 1:numel(uniqueN)
        n = uniqueN(kn);
        if strcmp(angularKind, 'cos')
            angFun = @(th) cos(n*th);
        else
            angFun = @(th) sin(n*th);
        end
        CnHandle = @(r) angularFourierCoeff(wFun, r, thetaGrid, dtheta, angFun);

        rowsAtN = find(modeList(:,1) == n);
        for kk = 1:numel(rowsAtN)
            k = rowsAtN(kk);
            j = modeList(k, 2);
            mode = circularPlateMode(n, j, R, nu);
            A(k) = integral(@(r) mode.radial(r) .* CnHandle(r) .* r, 0, R);
        end
    end
end

function val = angularFourierCoeff(wFun, r, thetaGrid, dtheta, angFun)
% C_n(r) = int_0^{2pi} w(r,theta) * angFun(theta) dtheta, evaluated at
% whatever r (scalar or vector, any shape) the caller (the radial
% adaptive quadrature in the main function) asks for. Vectorized over r:
% builds an (numel(r) x nTheta) grid via ndgrid, evaluates wFun ONCE on
% the whole grid, and reduces over theta with a single weighted sum --
% avoids an explicit loop over r, which matters since `integral` may
% query many r-nodes per call.
    rShape = size(r);
    rCol = r(:);
    [Rg, Tg] = ndgrid(rCol, thetaGrid); % numel(r) x nTheta
    Wg = wFun(Rg, Tg);
    val = sum(Wg .* angFun(Tg), 2) * dtheta;
    val = reshape(val, rShape);
end