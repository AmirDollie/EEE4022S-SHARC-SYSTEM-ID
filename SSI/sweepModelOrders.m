function results = sweepModelOrders(P, l, dt, nMax)
%SWEEPMODELORDERS Runs the already-validated extractStateSpace ->
%modalParameters chain at every model order from 1 to nMax, deduplicating
%complex-conjugate pole pairs so each oscillatory mode appears once.
%
%   RESULTS = SWEEPMODELORDERS(P, L, DT, NMAX)
%
%   P    - the Hankel projection matrix from hankelProjection.m
%   L    - number of output channels (matches buildHankelMatrix's l)
%   DT   - sampling period (seconds)
%   NMAX - highest model order to sweep to
%
%   Returns RESULTS, a 1 x nMax struct array. RESULTS(n) has fields:
%   frequency  - k x 1 natural frequencies (rad/s), k = number of
%                DISTINCT oscillatory modes at this order (i.e. after
%                keeping only one member of each conjugate pair)
%   damping    - k x 1 damping ratios, same ordering as frequency
%   modeShapes - l x k mode shapes, columns matching frequency/damping
%   poles      - k x 1 the continuous-time pole itself (complex), the
%                representative (Im>0) member of each conjugate pair
%
%   DEDUPLICATION: modalParameters.m returns every eigenvalue of the
%   identified A, including both members of each complex-conjugate pair
%   (this is correct and expected -- see its own docstring). A genuine
%   oscillatory mode corresponds to exactly ONE such pair, so plotting or
%   comparing both members separately would double-count every real mode
%   as two artificially-distinct poles. This function keeps only the
%   member with strictly positive imaginary part as that mode's single
%   representative, applied consistently at every order so later
%   order-to-order comparisons are comparing "one entry per mode" at
%   both ends.
%
%   NON-OSCILLATORY POLES: uses modalParameters' own isOscillatory flag
%   to decide what to drop, rather than testing imag(lambda_c) directly.
%   This matters because a real, NEGATIVE discrete eigenvalue produces a
%   spurious +i*pi/dt term after log() (MATLAB's complex branch cut),
%   which looks exactly like a genuine oscillatory pole with no
%   conjugate partner if you only look at lambda_c. Confirmed this was a
%   real failure mode empirically: it caused the mode count to exceed
%   the hard floor(n/2) bound on distinct conjugate pairs at odd orders
%   before this fix. Testing isOscillatory (which checks imag(lambda),
%   the discrete eigenvalue, before the log/branch-cut step) avoids
%   being fooled by this artifact.

    results(nMax) = struct('frequency', [], 'damping', [], ...
                            'modeShapes', [], 'poles', []);

    for n = 1:nMax
        [A, C, ~] = extractStateSpace(P, l, n);
        [omega, zeta, Phi, lambda_c, isOscillatory] = modalParameters(A, C, dt);

        keepMask = isOscillatory & (imag(lambda_c) > 0);

        results(n).frequency  = omega(keepMask);
        results(n).damping    = zeta(keepMask);
        results(n).modeShapes = Phi(:, keepMask);
        results(n).poles      = lambda_c(keepMask);
    end
end