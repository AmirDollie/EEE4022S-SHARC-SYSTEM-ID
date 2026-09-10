function MACval = computeMAC(Phi_a, Phi_b)
%COMPUTEMAC Modal Assurance Criterion between two (possibly complex) mode
%shape vectors or matrices (Allemang & Brown, 1982).
%
%   usage MACVAL = COMPUTEMAC(PHI_A, PHI_B)
%
%   PHI_A - l x na matrix of mode shape columns
%   PHI_B - l x nb matrix of mode shape columns
%
%   Returns an na x nb matrix, MACVAL(p,q) = MAC between column p of
%   PHI_A and column q of PHI_B:
%
%       MAC(phi_a,phi_b) = |phi_a^H * phi_b|^2 / [(phi_a^H*phi_a)*(phi_b^H*phi_b)]
%
%   in [0,1]. 1 indicates the two shapes are consistent up to an
%   arbitrary complex scale factor (exactly what output-only
%   identification leaves undetermined -> see the note in
%   modalParameters.m); 0 indicates orthogonal, unrelated shapes.
%
%   NOTE: the numerator is squared. A version of this found online
%   (referenced during this project) omitted the square, silently
%   computing sqrt(MAC) instead, which inflates every value toward 1
%   and makes any literature-stated threshold (e.g. Peeters & De Roeck's
%   99% mode shape consistency criterion) compare against the wrong
%   scale. Verified against this specific error in the tester below,
%   not just checked "looks plausible". Shame man I will not mention their
%   name but eish.

%   Once again, a huuuuge text bubble above for the simplest, intuitive bit
%   of poorly written code. I am seriously debating getting Claude to write
%   these mini functions at this point lol.

    [~, na] = size(Phi_a);
    [~, nb] = size(Phi_b);

    MACval = zeros(na, nb);
    for p = 1:na
        for q = 1:nb
            num = abs(Phi_a(:,p)' * Phi_b(:,q))^2;
            den = (Phi_a(:,p)' * Phi_a(:,p)) * (Phi_b(:,q)' * Phi_b(:,q));
            MACval(p,q) = num / den;
        end
    end
end