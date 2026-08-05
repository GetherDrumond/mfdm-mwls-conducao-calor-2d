function [pts, bc] = carregar_malha_gmsh(arquivo)
% CARREGAR_MALHA_GMSH  Le o arquivo perfilUDATA.txt gerado pelo Gmsh
%
%  [pts, bc] = carregar_malha_gmsh('perfilUDATA.txt')
%
%  Formato do arquivo:
%    Linha 1:  N  (numero de nos)
%    Linhas 2..N+1:  x  y  bc_flag
%      bc_flag = 1 -> no de contorno (Dirichlet)
%      bc_flag = 0 -> no interior
%
%  Saidas:
%    pts  - [N x 2] coordenadas (x, y)
%    bc   - [N x 1] flags de contorno

fid = fopen(arquivo, 'r');
if fid < 0
  error('Arquivo nao encontrado: %s', arquivo);
end

N = fscanf(fid, '%d', 1);
dados = fscanf(fid, '%f %f %d', [3, N])';
fclose(fid);

pts = dados(:, 1:2);
bc  = dados(:, 3);

fprintf('  Malha lida: %d nos  (contorno: %d, interior: %d)\n', ...
        N, sum(bc==1), sum(bc==0));
end
