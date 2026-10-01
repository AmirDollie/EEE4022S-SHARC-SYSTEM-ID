function [T, f, J, dT] = transmissibilityFromField(H, dH, layout, ref)
%TRANSMISSIBILITYFROMFIELD  Transmissibilities, stacked features and their Jacobian for any layout.
%   [T, f, J, dT] = transmissibilityFromField(H, dH, layout, ref)
%
%   H       nP x nF complex FRF at all candidate points (emmField.H at the bins)
%   dH      nP x nF x nTheta derivatives (emmField.dH), or [] (then J, dT empty)
%   layout  indices of the sensors used (into the rows of H), any order
%   ref     the reference, as an index INTO layout
%
%   T       m x nF, T_j = H_j / H_r for the non-reference sensors in ASCENDING layout position
%           (the extractTransmissibility / stackTransmissibility convention)
%   f       stacked features (stackTransmissibility: per frequency [Re T; Im T])
%   J       numel(f) x nTheta, d f / d theta, from dT_j = (dH_j H_r - H_j dH_r) / H_r^2
%   dT      m x nF x nTheta complex
%
% Lives in Inverse Mapping/Parameter Inversion/Design and Robustness/Engine/.

layout = layout(:).';
nS = numel(layout);
others = setdiff(1:nS, ref);
Hl = H(layout, :);
Hr = Hl(ref, :);
T = Hl(others, :) ./ Hr;
f = stackTransmissibility(T);
J = []; dT = [];
if isempty(dH), return; end
nT = size(dH, 3);
m = numel(others); nF = size(H, 2);
dT = complex(zeros(m, nF, nT));
J = zeros(2 * m * nF, nT);
for i = 1:nT
    dHl = dH(layout, :, i);
    dT(:, :, i) = (dHl(others, :) .* Hr - Hl(others, :) .* dHl(ref, :)) ./ Hr.^2;
    J(:, i) = stackTransmissibility(dT(:, :, i));
end
end