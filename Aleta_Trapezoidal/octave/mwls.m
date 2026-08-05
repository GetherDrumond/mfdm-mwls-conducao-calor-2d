function [M] = mwls(x, y, X, Y, S, m, p, g)
%MWLS Formulas de diferencas por Minimos Quadrados Moveis Ponderados.
%
%   M = mwls(x, y, X, Y, S, m)        p = 2, g = 0.1
%   M = mwls(x, y, X, Y, S, m, p, g)
%
%   Retorna M (s x m), s = 6 para p = 2. Cada linha aproxima um operador
%   a partir dos m valores nodais da estrela:
%       M(1,:) -> T          M(2,:) -> dT/dx      M(3,:) -> dT/dy
%       M(4,:) -> d2T/dx2    M(5,:) -> d2T/dxdy   M(6,:) -> d2T/dy2
%
%   CORRECAO (Bug 5) -- normalizacao das coordenadas locais.
%   A base [1, h, k, h^2/2, h*k, k^2/2] tem termos de ordens de grandeza
%   diferentes, e cond(P'*W*P) ~ dm^-6. Sem normalizar, o erro CRESCE sob
%   refino da nuvem (divergencia silenciosa: so aparece em estudo de malha).
%   Aqui adimensionaliza-se por dm = max(d) e as derivadas sao reescaladas
%   ao final por dm^-1 e dm^-2.

  if nargin < 8, g = 0.1; end
  if nargin < 7, p = 2;   end

  s = 6;
  if p == 1, s = 3; end

  idx = S(1:m,1);
  h = X(idx)' - x;
  k = Y(idx)' - y;
  d = sqrt(h.^2 + k.^2);

  % ---- adimensionalizacao ----
  dm = max(d);
  if dm < 1e-14, dm = 1.0; end
  hn = h/dm;  kn = k/dm;  dn = d/dm;

  geff = max(g, 1e-6);          % evita peso infinito quando d = 0

  z = [ones(1,m); hn; kn; 0.5*hn.^2; hn.*kn; 0.5*kn.^2]';
  P = z(:, 1:s);

  W = diag(1 ./ ((dn.^2 + geff^4 ./ (dn.^2 + geff^2 + 1e-15)).^(p+1) + 1e-15));

  A = P' * W * P;
  A = A + 1e-12 * trace(A) * eye(s);      % regularizacao de Tikhonov

  Mn = A \ (P' * W);                      % operador \ e nao ^-1

  % ---- volta as unidades fisicas ----
  esc_all = [1; 1/dm; 1/dm; 1/dm^2; 1/dm^2; 1/dm^2];
  M = diag(esc_all(1:s)) * Mn;
end
