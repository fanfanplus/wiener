function [ A, C, INDEX, ERRORS ] = spiker( B, LA )
% SPIKER Is a more efficient spiking-filter subroutine
%      It make use of Simpson sideways recursion as performed by SIDE

%   Inputs:    B            input sequence
%                  LA          length of filter
%   Outputs:  A           filter for optimum spike position
%                  C            actual output for optimum spike position
%                  INDEX    optimum spike position as MATLAB subscript 
%                  ERRORS  average squared error from 1 to LC

% Author: Felix Zhang
% Last modified: 2018-4-9

% References:
% [1] M. T. Silvia, and E. A. Robinson (1979) "Deconvolution of Geophysical Time
% Series in the Exploration for Oil and Natural Gas".

LB = length(B);
LC = LA+LB-1;
[R, lags] = xcorr(B, B, LA-1);
R = R(lags>=0);
for I = 1: LC
    SPACE = impuls(LC, I);
    [G, lags] = xcorr(SPACE, B, LA-1);
    G = G(lags>=0);
    if (I-1) <= 0
        [A, PEOC] = eureka(R, G);
        Q = dot(A, G);
        C = conv(A, B);
        ERRORS(I) = 1.0 - Q;
    else
        A = side(G, A, PEOC, R);
        Q = dot(A, G);
        C = conv(A, B);
        ERRORS(I) = 1.0 - Q;
    end
end

[EMIN, INDEX] = min(ERRORS);
SPACE = impuls(LC, INDEX);
[A, C] = shape(B, SPACE, LA);
end