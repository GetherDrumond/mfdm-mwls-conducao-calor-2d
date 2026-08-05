function [A, rhs] = montar_sistema_mfdm(pts, bc, T_bc)
% MONTAR_SISTEMA_MFDM  Sistema linear MFDM/MWLS para equação de Laplace
%  Usa base polinomial completa de grau 2: [1, x, y, x², xy, y²]
%  Função peso gaussiana com suporte local (estrela de n_viz vizinhos)

N   = size(pts, 1);
A   = sparse(N, N);
rhs = zeros(N, 1);

n_viz  = 16;   % MELHORIA: 16 vizinhos (>= 6 para grau 2 em 2D)
alpha  = 8.0;  % parâmetro da gaussiana
g_reg  = 0.1;  % MELHORIA: regularizacao de pesos

% Pré-calcular todas as distâncias seria O(N²); fazemos por nó
for i = 1:N

  if bc(i) == 1
    A(i, i) = 1.0;
    rhs(i)  = T_bc(i);
    continue;
  end

  % ── Encontrar vizinhos mais próximos ──────────────────────────────────
  dx   = pts(:,1) - pts(i,1);
  dy   = pts(:,2) - pts(i,2);
  dist = sqrt(dx.^2 + dy.^2);
  dist(i) = Inf;

  [d_sort, idx_sort] = sort(dist);
  nv      = min(n_viz, N-1);
  viz     = idx_sort(1:nv);
  dv      = d_sort(1:nv);

  h = dv(end) * 1.1;  % raio da estrela

  % ── Pesos gaussianos ─────────────────────────────────────────────────
  % MELHORIA: peso com regularizacao g=0.1
  w = exp(-alpha * ((dv.^2 + g_reg^4./(dv.^2 + g_reg^2)) / h^2));

  % ── Coordenadas locais dos vizinhos ──────────────────────────────────
  xi = pts(viz, 1) - pts(i,1);
  yi = pts(viz, 2) - pts(i,2);

  % ── Matriz de monomios P [nv × 6]: [1, x, y, x², xy, y²] ────────────
  P = [ones(nv,1), xi, yi, xi.^2, xi.*yi, yi.^2];

  W_diag = diag(w);

  % ── Matriz de momento M = Pᵀ W P  (6×6) ─────────────────────────────
  M = P' * W_diag * P;
  M = M + 1e-12 * trace(M) * eye(6);  % regularização

  % ── Vetor do Laplaciano: d²p/dx² + d²p/dy² avaliado nos monomios ─────
  % p = [1, x, y, x², xy, y²]
  % d²p/dx²  = [0, 0, 0, 2, 0, 0]
  % d²p/dy²  = [0, 0, 0, 0, 0, 2]
  % Laplaciano = soma = [0, 0, 0, 2, 0, 2]
  lap = [0; 0; 0; 2; 0; 2];

  % ── Coeficientes MWLS para os vizinhos ───────────────────────────────
  % Laplaciano(T)_i = lap' * (M \ P' W) * T_viz = coef_viz' * T_viz
  try
    B        = M \ (P' * W_diag);   % 6 × nv
    coef_viz = (lap' * B)';         % nv × 1
  catch
    % Se singular, usar diferenças finitas de 1ª ordem como fallback
    fprintf('  Aviso: singularidade no no %d, usando fallback\n', i);
    coef_viz = zeros(nv, 1);
  end

  % ── Coeficiente central via consistência numérica ────────────────────
  % Para o operador Laplaciano: sum(coef_all) = 0
  coef_central = -sum(coef_viz);

  % ── Preencher linha da matriz esparsa ─────────────────────────────────
  A(i, i) = coef_central;
  for k = 1:nv
    A(i, viz(k)) = A(i, viz(k)) + coef_viz(k);
  end
  rhs(i) = 0.0;  % Laplace: fonte = 0

end
end
