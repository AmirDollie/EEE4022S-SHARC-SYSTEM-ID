function f = stackTransmissibility(T)
%STACKTRANSMISSIBILITY The single, frozen stacking convention for
%transmissibility feature vectors. Used by BOTH the model side
%(transmissibilityFeatureVector) and the data side (estimates from
%extractTransmissibility), and by every covariance / residual built on
%them, so model and data can never be stacked differently.
%
%   F = STACKTRANSMISSIBILITY(T)
%
%   T - m x nF complex (rows = non-reference channels in ascending channel
%       order, columns = frequencies), e.g. est.T from extractTransmissibility.
%   F - (2 m nF) x 1 real column, FREQUENCY-MAJOR: one block of 2m entries
%       per frequency,
%         [ Re T_1(w_1); ...; Re T_m(w_1); Im T_1(w_1); ...; Im T_m(w_1);
%           Re T_1(w_2); ... ]
%       i.e. F((k-1)*2m + j) = Re T(j,k) and F((k-1)*2m + m + j) = Im T(j,k).
%       Each 2m-block is exactly one frequency, so a per-frequency covariance
%       Sigma_k (2m x 2m) assembles as Sigma = blkdiag(Sigma_1, ..., Sigma_nF)
%       with no index permutation.
%   Inverse: unstackTransmissibility(F, m).
%
%   A NaN in either part of T(j,k) (e.g. a flagged bin) makes BOTH its Re
%   and Im entries NaN in F (assigning a plain NaN to a complex array leaves
%   the imaginary part 0, which would otherwise enter a residual as data).
%
%   Lives in Transmissibility/.
if ~(isnumeric(T) && ismatrix(T) && ~isempty(T))
    error('stackTransmissibility:badT', 'T must be a nonempty m x nF numeric matrix.');
end
Z = [real(T); imag(T)];                    % 2m x nF
bad = isnan(real(T)) | isnan(imag(T));
Z([bad; bad]) = NaN;
f = Z(:);                                  % frequency-major
end