%% perfilU_main_single_figure.m  (v2 corrigida)
%  Conducao de calor 2D - Perfil U - MFDM/MWLS
%  TODOS OS GRAFICOS EM UMA UNICA FIGURA (layout 2x2)
%  Malha carregada do perfilUDATA.txt (gerado pelo Gmsh)
%
%  Correcoes:
%    - Visualizacao com trisurf + view(2) (sem griddata na cavidade)
%    - Variavel Wdom para evitar conflito
%    - Protecao de limites de caxis contra zero

clear; clc; close all;
fprintf('=== Perfil U - MFDM/MWLS (Unica Figura 2x2) ===\n\n');

Wdom   = 4.0;
H      = 3.0;
th     = 0.6;
k_cond = 1.0;

%% 1. Carregar malha do perfilUDATA.txt
fprintf('1) Carregando perfilUDATA.txt...\n');
[pts, bc] = carregar_malha_gmsh('perfilUDATA.txt');
N = size(pts, 1);

xmin = min(pts(:,1));
xmax = max(pts(:,1));
T_bc = 500.0 * (pts(:,1) - xmin) / (xmax - xmin);

fprintf('2) Montando sistema MFDM (%d nos)...\n', N);
[A, rhs] = montar_sistema_mfdm(pts, bc, T_bc);

fprintf('3) Resolvendo sistema linear...\n');
T = A \ rhs;
fprintf('   Tmin = %.2f C    Tmax = %.2f C\n', min(T), max(T));

fprintf('4) Calculando fluxos de calor...\n');
qx = zeros(N,1);
qy = zeros(N,1);
n_flux = 10;
for i = 1:N
  dx   = pts(:,1) - pts(i,1);
  dy   = pts(:,2) - pts(i,2);
  dist = sqrt(dx.^2 + dy.^2);
  dist(i) = Inf;
  [dv, idx] = sort(dist);
  viz = idx(1:min(n_flux, N-1));
  dv  = dv(1:min(n_flux, N-1));
  xi  = pts(viz,1) - pts(i,1);
  yi  = pts(viz,2) - pts(i,2);
  dT  = T(viz) - T(i);
  h   = max(dv) * 1.1;
  w   = exp(-6*(dv/h).^2);
  Wm  = diag(w);
  Pm  = [xi, yi];
  try
    g = (Pm'*Wm*Pm + 1e-12*eye(2)) \ (Pm'*Wm*dT);
    qx(i) = -k_cond * g(1);
    qy(i) = -k_cond * g(2);
  catch
    qx(i) = 0; qy(i) = 0;
  end
end
q_mag = sqrt(qx.^2 + qy.^2);
fprintf('   |q|max = %.2f W/m2\n', max(q_mag));

fprintf('5) Triangulando e filtrando vazio...\n');
tri_all = delaunay(pts(:,1), pts(:,2));
cx = (pts(tri_all(:,1),1) + pts(tri_all(:,2),1) + pts(tri_all(:,3),1)) / 3;
cy = (pts(tri_all(:,1),2) + pts(tri_all(:,2),2) + pts(tri_all(:,3),2)) / 3;
% FIX: filtro correto da cavidade
no_vazio = (cx > th + 1e-9) & (cx < Wdom - th - 1e-9) & (cy > th + 1e-9);
tri = tri_all(~no_vazio, :);
fprintf('   %d triangulos apos filtro\n', size(tri,1));

%% Contorno do perfil U para plotagem
xe = [0, Wdom, Wdom, Wdom-th, Wdom-th, th, th, 0, 0];
ye = [0, 0,    H,    H,       th,      th, H,  H, 0];

%% Colormap RdBu simetrico
nc = 256; half = nc/2;
r1 = linspace(0.017, 1.0, half)'; g1 = linspace(0.114, 1.0, half)'; b1 = linspace(0.388, 1.0, half)';
r2 = linspace(1.0, 0.403, half)'; g2 = linspace(1.0, 0.000, half)'; b2 = linspace(1.0, 0.122, half)';
cmap_rdbu = [r1, g1, b1; r2, g2, b2];

%% Limites simetricos para fluxos (com protecao)
lim_qx  = max(abs(qx))  * 1.02;
lim_qy  = max(abs(qy))  * 1.02;
lim_qmg = max(q_mag)    * 1.02;

if lim_qx  < 1e-10; lim_qx  = 1; end
if lim_qy  < 1e-10; lim_qy  = 1; end
if lim_qmg < 1e-10; lim_qmg = 1; end

campos  = {T,       qx,        qy,        q_mag};
vmins   = {min(T), -lim_qx,  -lim_qy,   0};
vmaxs   = {max(T),  lim_qx,   lim_qy,   lim_qmg};
cmaps   = {jet(256), cmap_rdbu, cmap_rdbu, hot(256)};
titulos = {'(a) TEMPERATURA T', '(b) FLUXO q_x', '(c) FLUXO q_y', '(d) MAGNITUDE |q|'};
cunits  = {'T [graus C]', 'q_x [W/m2]', 'q_y [W/m2]', '|q| [W/m2]'};

fprintf('6) Gerando figura unica 2x2...\n');
figure('Position', [50 50 1280 920], 'Name', 'Perfil U - Analise Termica');

for p = 1:4
    subplot(2, 2, p);

    % FIX: trisurf com Z=0, CData=campo, view(2)
    hs = trisurf(tri, pts(:,1), pts(:,2), zeros(N,1), campos{p}, ...
                 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);

    colormap(gca, cmaps{p});
    caxis([vmins{p}, vmaxs{p}]);
    cb = colorbar;
    tks = linspace(vmins{p}, vmaxs{p}, 6);
    set(cb, 'Ticks', tks, ...
        'TickLabels', arrayfun(@(v) sprintf('%.1f', v), tks, 'UniformOutput', false));
    ylabel(cb, cunits{p}, 'FontSize', 10);

    hold on;
    plot(xe, ye, 'k-', 'LineWidth', 2.0);
    hold off;

    title(titulos{p}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x [m]', 'FontSize', 10);
    ylabel('y [m]', 'FontSize', 10);
    axis equal tight;
    xlim([-0.2 Wdom+0.2]); ylim([-0.2 H+0.2]);
    set(gca, 'FontSize', 9);
    grid on; box on;
end

sgtitle('Perfil U - Analise de Transferencia de Calor  |  MFDM/MWLS  |  BC: T = 500*x/W', ...
        'FontSize', 13, 'FontWeight', 'bold');

fprintf('\n=== FIGURA UNICA GERADA ===\n');
fprintf('  (1,1) Temperatura T [graus C]\n');
fprintf('  (1,2) Fluxo de calor qx [W/m2]\n');
fprintf('  (2,1) Fluxo de calor qy [W/m2]\n');
fprintf('  (2,2) Magnitude |q| [W/m2]\n');
