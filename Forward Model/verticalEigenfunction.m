%% Insanely simple function, encapsulates the vertical bases for disc covered as well as open

function phi = verticalEigenfunction(z,xi,gamma)

% This is the normalised vertical eigenfunction
% Shared for open water phi0_m(z)
% and for disc: phi_m(z), phi_{-j}(z)

%Note if gamma = 0, then its for the open water case
% see phi(z) = cos(z+1)/cos(1-gamma) :)

% Note that this file doesn't distinguish which domain!
% Call it correctly...

%  Z may be a vector, xi is a scalar
% phi may be complex

phi = cos(xi .* (z+1)) ./ cos(xi .* (1-gamma));
end