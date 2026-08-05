function [pts, bc] = gerar_pontos_U(n_base, th, W, H)
% GERAR_PONTOS_U  Gera nuvem de pontos para o Perfil U
%
%   [pts, bc] = gerar_pontos_U(n_base, th, W, H)
%
%   Parametros:
%     n_base - densidade da malha (numero de divisoes na menor dimensao)
%     th     - espessura das paredes e base do U [m]
%     W      - largura total do U [m]
%     H      - altura total do U [m]
%
%   Saidas:
%     pts - [N x 2] coordenadas (x, y)
%     bc  - [N x 1] flags de contorno (1 = Dirichlet, 0 = interior)
%
%   O Perfil U consiste em:
%     - Base:          0 <= x <= W,       0 <= y <= th
%     - Parede esq.:   0 <= x <= th,      0 <= y <= H
%     - Parede dir.:   W-th <= x <= W,    0 <= y <= H
%     - Cavidade:      th < x < W-th,     th < y  (vazio)

  % Espacamento da grade
  ds = th / max(n_base, 2);

  % Grade regular cobrindo o retangulo [0,W] x [0,H]
  xv = 0:ds:W;
  yv = 0:ds:H;
  % Garantir que os limites exatos estejam incluidos
  if xv(end) < W; xv(end+1) = W; end
  if yv(end) < H; yv(end+1) = H; end

  [Xg, Yg] = meshgrid(xv, yv);
  Xg = Xg(:);
  Yg = Yg(:);

  % Mascara: manter apenas pontos dentro do dominio U
  % O dominio U = retangulo total MENOS a cavidade
  % Cavidade: th < x < W-th  E  y > th
  tol = ds * 1e-6;
  na_cavidade = (Xg > th + tol) & (Xg < W - th - tol) & (Yg > th + tol);
  no_dominio = ~na_cavidade;

  Xg = Xg(no_dominio);
  Yg = Yg(no_dominio);
  Np = length(Xg);

  % Remover pontos duplicados (por seguranca)
  pts_raw = [Xg, Yg];
  pts_raw = unique(round(pts_raw / (ds*1e-4)) * (ds*1e-4), 'rows');
  Np = size(pts_raw, 1);

  % Classificar contorno vs interior
  bc = zeros(Np, 1);
  for i = 1:Np
    x = pts_raw(i, 1);
    y = pts_raw(i, 2);

    % Borda externa
    if abs(x) < tol || abs(x - W) < tol
      bc(i) = 1;  % paredes laterais externas
    elseif abs(y) < tol
      bc(i) = 1;  % base inferior
    elseif abs(y - H) < tol && (x <= th + tol || x >= W - th - tol)
      bc(i) = 1;  % topo das paredes
    end

    % Borda interna (cavidade)
    if abs(x - th) < tol && y >= th - tol
      bc(i) = 1;  % parede interna esquerda
    elseif abs(x - (W - th)) < tol && y >= th - tol
      bc(i) = 1;  % parede interna direita
    elseif abs(y - th) < tol && x >= th - tol && x <= W - th + tol
      bc(i) = 1;  % base da cavidade (topo da base do U)
    end
  end

  pts = pts_raw;

  fprintf('  gerar_pontos_U: %d nos (contorno: %d, interior: %d)\n', ...
          Np, sum(bc==1), sum(bc==0));
  fprintf('  Espacamento ds = %.4f m\n', ds);
end
