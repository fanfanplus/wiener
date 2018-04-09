function [ D ] = impuls( LD, K )
%IMPULS Is Kronecker impulse function

% Author: Felix Zhang
% Last modified: 2018-4-9

% References:
% [1] M. T. Silvia, and E. A. Robinson (1979) "Deconvolution of Geophysical Time
% Series in the Exploration for Oil and Natural Gas".

for I = 1: LD
    D(I) = 0;
    if I == K
        D(I) = 1;
    end
end
end

