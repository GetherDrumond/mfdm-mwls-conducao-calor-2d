function nrm = normais_contorno(X, Y, Elin, Ttri, N)
%NORMAIS_CONTORNO Normal externa unitaria nos nos de contorno.
%
%   nrm = normais_contorno(X, Y, Elin, Ttri, N)
%   Elin: (nl x 3) com [codigo  n1  n2] -- arestas de contorno do .msh
%   Ttri: (nt x 3) conectividade triangular
%   nrm : (N x 2) normal externa por no (zeros nos nos internos)
%
%   CORRECAO (Bug 6) -- orientacao em fronteiras internas.
%   Orientar a normal pelo centroide do dominio so vale para contorno
%   convexo. Em um FURO a normal externa do dominio aponta para DENTRO do
%   furo, e o teste do centroide a inverte (produziu K = 50 no lugar de 2).
%   Criterio correto: cada aresta de contorno pertence a exatamente um
%   triangulo; a normal externa aponta no sentido OPOSTO ao terceiro vertice.

  nt = size(Ttri, 1);

  % ---- tabela aresta -> triangulo ----
  E = zeros(3*nt, 3);                       % [nmin  nmax  itri]
  for t = 1:nt
    v = Ttri(t,:);
    pr = [v(1) v(2); v(2) v(3); v(3) v(1)];
    for e = 1:3
      E(3*(t-1)+e, :) = [min(pr(e,:)), max(pr(e,:)), t];
    end
  end
  chave = E(:,1)*(N+1) + E(:,2);            % chave unica por aresta

  nrm = zeros(N, 2);

  for i = 1:size(Elin, 1)
    a = Elin(i,2);  b = Elin(i,3);
    kk = min(a,b)*(N+1) + max(a,b);
    pos = find(chave == kk, 1);
    if isempty(pos), continue; end

    t  = E(pos,3);
    v  = Ttri(t,:);
    op = v(v ~= a & v ~= b);                % vertice oposto a aresta
    if isempty(op), continue; end
    op = op(1);

    tx = X(b) - X(a);  ty = Y(b) - Y(a);
    L  = sqrt(tx^2 + ty^2);
    if L < 1e-14, continue; end

    nx =  ty/L;  ny = -tx/L;                % normal candidata
    xm = (X(a) + X(b))/2;  ym = (Y(a) + Y(b))/2;

    if (X(op) - xm)*nx + (Y(op) - ym)*ny > 0   % apontou para o interior
      nx = -nx;  ny = -ny;
    end

    nrm(a,:) = nrm(a,:) + [nx ny];
    nrm(b,:) = nrm(b,:) + [nx ny];          % media nos vertices
  end

  L = sqrt(sum(nrm.^2, 2));
  L(L < 1e-14) = 1;
  nrm = nrm ./ [L L];
end
