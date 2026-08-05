function [qx, qy, qm] = fluxos_aleta(N, X, Y, T, kcond, m, g, furos)
%FLUXOS_ALETA Recupera os fluxos nodais pela Lei de Fourier via MWLS.
%
%   [qx, qy, qm] = fluxos_aleta(N, X, Y, T, kcond, m, g, furos)
%       qx = -k*dT/dx    qy = -k*dT/dy    qm = sqrt(qx^2 + qy^2)
%
%   Usa as linhas 2 e 3 da matriz M (derivadas de primeira ordem) com a
%   mesma estrela por visibilidade empregada na montagem do sistema.

  if nargin < 6 || isempty(m), m = 16;  end
  if nargin < 7 || isempty(g), g = 0.1; end
  if nargin < 8,               furos = []; end

  qx = zeros(N,1);  qy = zeros(N,1);

  for i = 1:N
    S = star_visivel(N, X(i), Y(i), X, Y, m, furos);
    M = mwls(X(i), Y(i), X, Y, S, m, 2, g);
    idx = S(1:m,1);
    qx(i) = -kcond * (M(2,:) * T(idx));
    qy(i) = -kcond * (M(3,:) * T(idx));
  end

  qm = sqrt(qx.^2 + qy.^2);
end
