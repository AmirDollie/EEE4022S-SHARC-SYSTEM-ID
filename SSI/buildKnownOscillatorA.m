function A = buildKnownOscillatorA(omega1, omega2, dt)
%BUILDKNOWNOSCILLATORA Constructs a real, exactly-known discrete-time
%state matrix representing two undamped forced tones at omega1 and
%omega2, for use as ground truth in Level 2C's sensor-placement study.
%
%   A = BUILDKNOWNOSCILLATORA(OMEGA1, OMEGA2, DT)
%
%   OMEGA1, OMEGA2 - the two known angular frequencies (rad/s)
%   DT             - sampling period (seconds)
%
%   Returns a real 4x4 block-diagonal rotation matrix,
%       A = blkdiag(R(omega1*dt), R(omega2*dt)),
%       R(theta) = [cos(theta) -sin(theta); sin(theta) cos(theta)],
%   with state x = [cos(omega1*t); sin(omega1*t); cos(omega2*t); sin(omega2*t)].
%
%   WHY THIS FORM:
%   Level 1 confirmed a forced tone has zeta~0 (eigenvalues
%   on the unit circle) essentially exactly. A real 2x2 rotation block is
%   the standard way to represent a complex-conjugate eigenvalue pair
%   e^{+-i*theta} using only real arithmetic -> eig(R(theta)) gives
%   exactly e^{+-i*theta}, confirmed in the tester below.
%
%   NO SSI IS USED to construct this matrix -> it is built directly from
%   the KNOWN frequencies, specifically to avoid the circularity of
%   needing sensor placement to run SSI, while needing SSI's output to
%   choose sensor placement. This A must be paired with C rows built by
%   buildCandidateCrows.m, which expects EXACTLY this state ordering
%   (cos1, sin1, cos2, sin2) -> changing the order here without updating
%   that function would silently misalign state and measurement.

    R1 = [cos(omega1*dt), -sin(omega1*dt); sin(omega1*dt), cos(omega1*dt)];
    R2 = [cos(omega2*dt), -sin(omega2*dt); sin(omega2*dt), cos(omega2*dt)];
    A = blkdiag(R1, R2);
end