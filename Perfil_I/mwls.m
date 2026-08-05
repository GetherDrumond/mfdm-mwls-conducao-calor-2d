function [M] = mwls(x,y,X,Y,S,m,p,g)
% MWLS  Moving Weighted Least Squares — shape function derivatives
%
%   M = mwls(x,y,X,Y,S,m,p,g)
%
%   Inputs:
%     x,y  - evaluation point
%     X,Y  - all nodal coordinates
%     S    - star (sorted neighbor indices + distances), from star()
%     m    - number of star nodes
%     p    - polynomial order (1 = linear, 2 = quadratic). Default: 2
%     g    - regularization parameter for weight function. Default: 0.1
%
%   Output:
%     M    - (s x m) matrix of shape function derivatives
%            For p=2: rows = [T, Tx, Ty, Txx, Txy, Tyy]
%            For p=1: rows = [T, Tx, Ty]

if nargin < 8
    g = 0.1;       % MELHORIA: g=0.1 em vez de 0 (regulariza pesos)
end
if nargin < 7
    p = 2;
end

s = 6;
if p == 1
    s = 3;
end

h = X(S(1:m,1))' - x;
k = Y(S(1:m,1))' - y;
d = sqrt(h.^2 + k.^2);

z = [ones(1,m); h; k; 0.5*h.^2; h.*k; 0.5*k.^2]';
P = z(:, 1:s);

% Peso com regularizacao: evita singularidade quando d -> 0
W = diag(1./((d.^2 + (g^4)./(d.^2 + g^2)).^(p+1)));

% Resolver com \ e regularizacao de Tikhonov
PtWP = P' * W * P;
reg  = 1e-12 * trace(PtWP) * eye(s);
M    = (PtWP + reg) \ (P' * W);
