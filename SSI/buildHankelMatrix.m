function [Yp_ref, Yf, H] = buildHankelMatrix(y, i, refIdx)
%BUILDHANKELMATRIX Stacks a multi-channel output time series into the
%block Hankel matrix used by SSI-DATA 
%
% example usage:   [YP_REF, YF, H] = BUILDHANKELMATRIX(Y, I, REFIDX)

% WHAT EVERYTHING MEANS:
% 1. Yp_ref:
%   this is the 'past'. It captures the last 'i' samples of output history.
%   Each column of this is this window. It only contains data from the
%   REFERENCE sensors, ie the sensors we trust most... 

% 2. Yf:
%   this is the 'future'. It captures the next 'i' samples from EVERY
%   sensor, not just references. Each column represents 'what happens right
%   after' the corresponding column in Yp_ref!

% In other words:
% Yp_ref is "recent history, from the sensors we trust most, used to 
% PREDICT"; Yf is "what happens next, at every sensor, that we're trying
%  to explain." The projection step that follows in other code blocks 
% is precisely the operation that measures 
% how much of Yf can be accounted for by Yp_ref, 
% and whatever's left over after removing that predictable part is the raw 
% material the rest of SSI turns into your identified system.




    [l, n] = size(y);

    if nargin < 3 || isempty(refIdx)
        refIdx = 1:l;
    end
    if any(refIdx < 1) || any(refIdx > l)
        error('buildHankelMatrix:badRefIdx', 'refIdx must index rows 1..%d of y.', l);
    end

    r = length(refIdx);
    j = n - 2*i + 1;
    if j < 1
        error('buildHankelMatrix:tooFewSamples', ...
              'Need n >= 2*i; got n=%d, i=%d (j=%d).', n, i, j);
    end

    yRef = y(refIdx, :);

    Yp_ref = zeros(r*i, j);
    for k = 0:(i-1)
        Yp_ref(k*r+1 : (k+1)*r, :) = yRef(:, k+1 : k+j);
    end

    Yf = zeros(l*i, j);
    for k = i:(2*i-1)
        Yf((k-i)*l+1 : (k-i+1)*l, :) = y(:, k+1 : k+j);
    end

    scale = 1/sqrt(j);
    Yp_ref = Yp_ref * scale;
    Yf = Yf * scale;

    H = [Yp_ref; Yf];
end