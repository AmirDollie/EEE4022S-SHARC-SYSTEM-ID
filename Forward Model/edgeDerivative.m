%% Closed form derivative phi'_m(-gamma)
%Use in defelction sum!
function dphi = edgeDerivative(xi, gamma)

% Note derivation;
% d/dz[ cos(xi(z+1))/cos(xi(1-gamma))]
    % = -xi*sin(k(z+1))/cos(k(1-gamma))

    %and at z = -gamma
    % = -xo*tan(xi(1-gamma))


%NOTE:
%The thesis starts by printing phi'_m(0) in a single plaec, but the
%pipeline as reconstructed only makes sense if we use phi'_m(-gamma), which
%is what this function uses. 

%If an issue, simple fix and only needs ammending here!

    dphi = -xi .* tan(xi .* (1 - gamma));
end
