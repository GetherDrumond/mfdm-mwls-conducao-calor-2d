function T = localMFDM_heat_aleta(N, X, Y, BCtipo, BCval, kcond, qvol, m, g, furos, nrm)
%LOCALMFDM_HEAT_ALETA Conducao de calor 2D estacionaria por MFDM/MWLS,
%   com suporte a condicoes de Dirichlet E de Neumann.
%
%   T = localMFDM_heat_aleta(N,X,Y,BCtipo,BCval,kcond,qvol,m,g,furos,nrm)
%
%   BCtipo(i): 0 = no interno
%              1 = Dirichlet  (T = BCval(i))
%              2 = Neumann    (dT/dn = 0, normal em nrm(i,:))
%   furos    : (nf x 3) [xc yc r] para o criterio de visibilidade; [] se nao houver
%   nrm      : (N x 2) normais externas (use normais_contorno.m)
%
%   Equacao resolvida nos nos internos:  k*(d2T/dx2 + d2T/dy2) + q = 0

  if nargin < 8  || isempty(m),     m     = 16;  end
  if nargin < 9  || isempty(g),     g     = 0.1; end
  if nargin < 10,                   furos = [];  end
  if nargin < 11,                   nrm   = zeros(N,2); end

  % montagem esparsa por triplas (I,J,V) -- muito mais rapido que A(i,j)=...
  nnzmax = N*m + N;
  I = zeros(nnzmax,1);  J = zeros(nnzmax,1);  V = zeros(nnzmax,1);
  cnt = 0;
  B = zeros(N,1);

  for i = 1:N

    if BCtipo(i) == 1                      % ---- Dirichlet ----
      cnt = cnt + 1;
      I(cnt) = i;  J(cnt) = i;  V(cnt) = 1;
      B(i) = BCval(i);
      continue;
    end

    S = star_visivel(N, X(i), Y(i), X, Y, m, furos);
    M = mwls(X(i), Y(i), X, Y, S, m, 2, g);
    idx = S(1:m,1);

    if BCtipo(i) == 2                      % ---- Neumann: n . grad(T) = 0 ----
      coef = nrm(i,1)*M(2,:) + nrm(i,2)*M(3,:);
      B(i) = 0;
    else                                   % ---- interno: laplaciano ----
      coef = M(4,:) + M(6,:);
      B(i) = -qvol/kcond;
    end

    for jj = 1:m
      cnt = cnt + 1;
      I(cnt) = i;  J(cnt) = idx(jj);  V(cnt) = coef(jj);
    end
  end

  A = sparse(I(1:cnt), J(1:cnt), V(1:cnt), N, N);
  T = A \ B;
end
