function C_j = buildCandidateCrows(phi1_tilde, phi2_tilde)
%BUILDCANDIDATECROWS Constructs one candidate sensor's row of C, given
%that sensor's (already normalised) modal response at two known
%frequencies, matching the state ordering of buildKnownOscillatorA.m.
%
%   C_J = BUILDCANDIDATECROWS(PHI1_TILDE, PHI2_TILDE)
%
%   PHI1_TILDE, PHI2_TILDE - scalar complex values: this candidate
%       sensor's (already-normalised) spatial response to mode 1 and
%       mode 2 respectively, e.g. from evaluateDeflection.
%
%   Returns C_J, a 1x4 real row vector, matching the state ordering
%   x = [cos(omega1*t); sin(omega1*t); cos(omega2*t); sin(omega2*t)]
%   from buildKnownOscillatorA.m.
%
%   DERIVATION:
% 
% 
%  the physical signal at this sensor is
%   Re{phi*e^{i*omega*t}} = Re(phi)*cos(omega*t) - Im(phi)*sin(omega*t).
%   Summing this over both modes and reading off the coefficient of each
%   state component gives
%       C_j = [Re(phi1), -Im(phi1), Re(phi2), -Im(phi2)].
%
%   THIS MAPPING WAS ITSELF THE SUBJECT OF A CAUGHT BUG during the design
%   of this experiment: an earlier draft left C_j as the raw 2-column
%   complex vector [phi1, phi2], which does not match the 4-dimensional
%   REAL state buildKnownOscillatorA.m actually produces. Getting this
%   wrong would make the Gramian describe a state-space model that does
%   not generate the actual signal syntheticSensorData.m would produce
%   for the same configuration.
%
%   NORMALISATION IS THE CALLER'S RESPONSIBILITY, deliberately: this
%   function does not normalise phi1_tilde/phi2_tilde itself, so the
%   caller can make an explicit, visible choice (e.g. normalise each
%   mode's field independently across the whole candidate grid before
%   calling this function, to isolate geometric distinguishability from
%   amplitude imbalance -> see the Level 2C design notes).

%   Truly outdone myself this time by typing a kajillion word explanation
%   for a single line, love to see it!

    C_j = [real(phi1_tilde), -imag(phi1_tilde), real(phi2_tilde), -imag(phi2_tilde)];
end