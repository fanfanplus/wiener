function [ F, A ] = eureka( R, G )
%EUREKA Solves the least square normal equation
%   F filter coefficients and A prediction-error operator

% Author: Felix Zhang
% Last modified: 2018-4-9

% References:
% [1] M. T. Silvia, and E. A. Robinson (1979) "Deconvolution of Geophysical Time
% Series in the Exploration for Oil and Natural Gas".

LR = length(R);
V = R(1);
A(1) = 1.0;
F(1) = G(1)/V;
if LR == 1
    return;
end
for L = 2: LR
    D = 0.0;
    Q = 0.0;
    L3 = L-1;
    for J = 1 : L3
         K = L-J+1;
         D = D + A(J)*R(K);
         Q = Q + F(J)*R(K);
    end
    C = D/V;
    if L == 2
        A(L) = -C;
        V = V-C*D;
        S = (Q-G(L))/V;
        for J = 1: L3
            K = L-J+1;
            F(J) = F(J)-S*A(K);
        end
        F(L) = -S;
        continue;
    end
    L1 = floor((L-2)/2);
    L2 = L1+1;
    if L2 < 2
        LT3 = L2+1;
        A(LT3) = A(LT3) - C*A(LT3);
        A(L) = -C;
        V = V-C*D;
        S = (Q-G(L))/V;
        for J = 1: L3
            K = L-J+1;
            F(J) = F(J)-S*A(K);
        end
        F(L) = -S;
        continue;
    end
    for J = 2 : L2
        HOLD = A(J);
        K = L-J+1;
        A(J) = A(J)-C*A(K);
        A(K) = A(K)-C*HOLD;
    end
    if (2*L1) == (L-2)
        A(L) = -C;
        V = V-C*D;
        S = (Q-G(L))/V;
        for J = 1: L3
            K = L-J+1;
            F(J) = F(J)-S*A(K);
        end
        F(L) = -S;
        continue;
    end
    LT3 = L2+1;
    A(LT3) = A(LT3) - C*A(LT3);
    A(L) = -C;
    V = V-C*D;
    S = (Q-G(L))/V;
    for J = 1: L3
        K = L-J+1;
        F(J) = F(J)-S*A(K);
    end
    F(L) = -S;
end
end

