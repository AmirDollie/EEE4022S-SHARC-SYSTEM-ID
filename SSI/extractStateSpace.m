function [A, C, singularValues] = extractStateSpace(P, l, n)
%EXTRACTSTATESPACE Recovers the state matrix A and output matrix C from
%the Hankel projection, at a chosen model order n
%
%  usage: [A, C, SINGULARVALUES] = EXTRACTSTATESPACE(P, L, N)
%
%   P - (l*i) x j projection matrix from hankelProjection.m
%   L - number of output channels (all channels, not just references ->
%       must match the l used when building P via buildHankelMatrix)
%   N - model order to truncate to (number of states to keep)
%
%   Returns:
%   A - n x n state matrix
%   C - l x n output matrix (first l rows of the observability matrix)
%   SINGULARVALUES - the full singular value vector of P, for inspecting
%       where a natural truncation point (if any) actually is, and for
%       driving a stabilization diagram (run this at several N and
%       compare) rather than trusting a single N blindly.
%
%   Method: SVD of P gives P = U*S*V'. Truncating to the top N singular
%   values/vectors gives a balanced realisation of the extended
%   observability matrix, Oi = U1*S1^(1/2). A is then recovered from the
%   shift-invariance of Oi (the fact that dropping the first block row of
%   Oi and dropping the last give two matrices related by exactly one
%   factor of A), via a least-squares (pseudo-inverse) solve, NOT by
%   forming a second, shifted projection matrix separately, since Oi's
%   own shift structure already contains everything needed.

    [U, S, ~] = svd(P);
    singularValues = diag(S);

    i = size(P,1) / l;
    if mod(i,1) ~= 0
        error('extractStateSpace:badDimensions', ...
              'size(P,1)=%d is not divisible by l=%d; check l matches buildHankelMatrix.', size(P,1), l);
    end

    U1 = U(:, 1:n);
    S1 = S(1:n, 1:n);

    Oi = U1 * sqrt(S1);        % (l*i) x n

    Oi_top    = Oi(1:end-l, :);   % rows 1 .. l*(i-1)   (all but last block)
    Oi_bottom = Oi(l+1:end, :);   % rows l+1 .. l*i     (all but first block)

    A = pinv(Oi_top) * Oi_bottom;
    C = Oi(1:l, :);
end