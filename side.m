function [ O ] = side( G, F, A, R )
%SIDE Is Simpson sideways recursion
%   SIDE is used in conjunction with EUREKA to recur sideways.
%   Inputs:    H     cross-correlation
%                  F     filter coefficients
%                  A    prediction-error operator coefficients
%                  R     autocorrelation from lag zero to lag LR-1
%   Outputs:  O    new filter coefficients

% Author: Felix Zhang
% Last modified: 2018-4-9

% References:
% [1] M. T. Silvia, and E. A. Robinson (1979) "Deconvolution of Geophysical Time
% Series in the Exploration for Oil and Natural Gas".

LF = length(F);
V = R(1);
H = G(1);
S = 0.0;
T = 0.0;
if LF == 1
    FLF = F(LF);
    W = (H-S+FLF*T)/V;
    if LF == 1
        O(1) = W;
        return;
    end
    for I = 2: LF
        J = LF-I+2;
        O(J) = F(J-1)+W*A(J)-FLF*A(I);
    end
    O(1) = W;
    return;
end

for I = 2: LF
    J = LF+2-I;
    S = S+F(I-1)*R(I);
    T = T+A(J)*R(I);
    V = V+A(I)*R(I);
end
FLF = F(LF);
W = (H-S+FLF*T)/V;
if LF == 1
    O(1) = W;
    return;
end

for I = 2: LF
    J = LF-I+2;
    O(J) = F(J-1)+W*A(J)-FLF*A(I);
end
O(1) = W;
end