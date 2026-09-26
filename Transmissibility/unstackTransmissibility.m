function T = unstackTransmissibility(f, m)
%UNSTACKTRANSMISSIBILITY Inverse of stackTransmissibility (frequency-major).
%
%   T = UNSTACKTRANSMISSIBILITY(F, M)
%
%   F - (2 m nF) real vector in the frequency-major layout of
%       stackTransmissibility; M - number of transmissibilities (rows of T).
%   T - m x nF complex.
%
%   Lives in Transmissibility/.
if ~(isnumeric(f) && isreal(f) && isvector(f))
    error('unstackTransmissibility:badF', 'F must be a real vector.');
end
if ~(isnumeric(m) && isscalar(m) && m >= 1 && m == round(m)) || mod(numel(f), 2 * m) ~= 0
    error('unstackTransmissibility:badM', 'numel(F) must be a multiple of 2*M.');
end
nF = numel(f) / (2 * m);
Z = reshape(f(:), 2 * m, nF);
T = complex(Z(1:m, :), Z(m + 1:end, :));
end