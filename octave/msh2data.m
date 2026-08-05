function [N, X, Y, BC, nt, Ttri, Elin] = msh2data(mshfile, outfile)
%MSH2DATA Converte malha Gmsh (.msh v2.2 ASCII) para o formato MFDM.
%
%   [N,X,Y,BC,nt,Ttri,Elin] = msh2data('aleta.msh')
%   Elin (nl x 3) = [codigo n1 n2] das arestas de contorno (p/ normais)
%   [...] = msh2data('aleta.msh', 'aletaDATA.txt')  tambem grava o arquivo
%
%   Formato de saida (mesmo do railroadrailDATA.txt):
%       N
%       x1  y1  BC1
%       ...
%       nt
%       n1  n2  n3
%       ...
%
%   Codigos de BC (vindos das Physical Curves do .geo):
%       0 = no interno
%       1 = topo   (T = 500 C)
%       2 = base   (T = 0 C)
%       3 = lados inclinados (isolado / Neumann)
%       4 = borda dos furos  (T = 0 C)
%
%   Regra de prioridade nos cantos: o menor codigo prevalece, ou seja,
%   um no na juncao topo/lado e classificado como topo (Dirichlet ganha
%   de Neumann, evitando sistema mal-posto).

  if nargin < 2, outfile = ''; end

  fid = fopen(mshfile, 'rt');
  if fid < 0
    error('msh2data: nao consegui abrir %s', mshfile);
  end

  X = []; Y = []; BC = []; Ttri = []; Elin = []; N = 0; nt = 0;

  while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);

    % ---------------- NOS ----------------
    if strcmp(line, '$Nodes')
      N = str2double(fgetl(fid));
      X = zeros(N,1); Y = zeros(N,1); BC = zeros(N,1);
      for i = 1:N
        v = sscanf(fgetl(fid), '%f');
        idx = v(1);            % tag do no (Gmsh numera a partir de 1)
        X(idx) = v(2);
        Y(idx) = v(3);
      end

    % ---------------- ELEMENTOS ----------------
    elseif strcmp(line, '$Elements')
      ne = str2double(fgetl(fid));
      Ttri = zeros(ne, 3);  Elin = zeros(ne, 3);
      kt = 0;  kl = 0;
      for i = 1:ne
        v = sscanf(fgetl(fid), '%d');
        etype = v(2);          % 1 = linha (2 nos), 2 = triangulo (3 nos)
        ntags = v(3);
        phys  = v(4);          % 1a tag = Physical group = nosso codigo BC
        conn  = v(4+ntags : end)';

        if etype == 1
          % elemento de linha -> aresta de contorno + marca os nos
          kl = kl + 1;
          Elin(kl,:) = [phys, conn(1), conn(2)];
          for k = 1:numel(conn)
            n = conn(k);
            if BC(n) == 0
              BC(n) = phys;
            else
              BC(n) = min(BC(n), phys);   % prioridade: menor codigo
            end
          end
        elseif etype == 2
          kt = kt + 1;
          Ttri(kt,:) = conn(1:3);
        end
      end
      Ttri = Ttri(1:kt, :);
      Elin = Elin(1:kl, :);
      nt = kt;
    end
  end
  fclose(fid);

  % ---------------- VERIFICACOES ----------------
  if N == 0,  error('msh2data: nenhum no lido de %s', mshfile); end
  if nt == 0, error('msh2data: nenhum triangulo lido de %s', mshfile); end
  if any(Ttri(:) < 1) || any(Ttri(:) > N)
    error('msh2data: conectividade fora de [1..%d]', N);
  end

  nb = sum(BC ~= 0);
  fprintf('Malha lida: %d nos (%d contorno, %d internos), %d triangulos\n', ...
          N, nb, N-nb, nt);
  for c = 1:4
    fprintf('   codigo %d: %d nos\n', c, sum(BC == c));
  end

  % ---------------- GRAVACAO ----------------
  if ~isempty(outfile)
    fo = fopen(outfile, 'wt');
    fprintf(fo, '%d\n', N);
    for i = 1:N
      fprintf(fo, '%.10f %.10f %d\n', X(i), Y(i), BC(i));
    end
    fprintf(fo, '%d\n', nt);
    for i = 1:nt
      fprintf(fo, '%d %d %d\n', Ttri(i,1), Ttri(i,2), Ttri(i,3));
    end
    fclose(fo);
    fprintf('Gravado: %s\n', outfile);
  end
end
