%% perfilU_main.m  (v6 — melhorias de precisao e visualizacao)
%  Conducao de calor 2D - Perfil U - MFDM/MWLS
%  Versao compativel com Octave - 4 subplots em uma unica figura
%
%  Melhorias v6 (orientador):
%    - n_flux = 16 (mais vizinhos para calculo de fluxo)
%    - g = 0.1 (regularizacao de pesos)
%    - Tratamento de flutuacao de ponto flutuante nos graficos
%    - Escalas fixas baseadas nos valores analiticos esperados

clear; clc; close all;
fprintf('=== Perfil U - MFDM/MWLS v6 (Octave - melhorias) ===\n\n');

Wdom   = 4.0;   % largura total
H      = 3.0;
th     = 0.6;
n_base = 14;
k_cond = 1.0;

fprintf('1) Gerando nos...\n');
[pts, bc] = gerar_pontos_U(n_base, th, Wdom, H);
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
n_flux = 16;       % MELHORIA: 16 em vez de 10 (mais vizinhos = mais preciso)
g_reg  = 0.1;      % MELHORIA: regularizacao de pesos
for i = 1:N
  dx   = pts(:,1) - pts(i,1);
  dy   = pts(:,2) - pts(i,2);
  dist = sqrt(dx.^2 + dy.^2);
  dist(i) = Inf;
  [dv, idx] = sort(dist);
  nv  = min(n_flux, N-1);
  viz = idx(1:nv);
  dv  = dv(1:nv);
  xi  = pts(viz,1) - pts(i,1);
  yi  = pts(viz,2) - pts(i,2);
  dT  = T(viz) - T(i);
  h   = max(dv) * 1.1;
  % MELHORIA: peso com regularizacao g=0.1
  w   = exp(-6*((dv.^2 + g_reg^4./(dv.^2 + g_reg^2))/(h^2)).^1);
  Wm  = diag(w);
  Pm  = [xi, yi];
  try
    g_coef = (Pm'*Wm*Pm + 1e-12*eye(2)) \ (Pm'*Wm*dT);
    qx(i) = -k_cond * g_coef(1);
    qy(i) = -k_cond * g_coef(2);
  catch
    qx(i) = 0; qy(i) = 0;
  end
end
q_mag = sqrt(qx.^2 + qy.^2);
fprintf('   |q|max = %.2f W/m2\n', max(q_mag));

% ============================================================
% MELHORIA: Tratamento de flutuacao de ponto flutuante
% ============================================================
qx_analitico = -k_cond * 500 / (xmax - xmin);
fprintf('\n   Valor analitico esperado: qx = %.4f W/m2\n', qx_analitico);
fprintf('   Range numerico obtido:    qx = [%.4f, %.4f]\n', min(qx), max(qx));
fprintf('   Range numerico obtido:    qy = [%.6e, %.6e]\n', min(qy), max(qy));

% Threshold: se qy < 1% de |qx_analitico|, tratar como zero
thr_qy = max(abs(qy));
if thr_qy < abs(qx_analitico) * 0.01
    fprintf('   qy tratado como zero (flutuacao < 1%% de |qx_analitico|)\n');
    qy_plot = zeros(size(qy));
    lim_qy = abs(qx_analitico) * 0.05;
else
    qy_plot = qy;
    lim_qy = max(abs(qy)) * 1.05;
end

fprintf('\n5) Triangulando e filtrando vazio...\n');
tri_all = delaunay(pts(:,1), pts(:,2));
cx = (pts(tri_all(:,1),1) + pts(tri_all(:,2),1) + pts(tri_all(:,3),1)) / 3;
cy = (pts(tri_all(:,1),2) + pts(tri_all(:,2),2) + pts(tri_all(:,3),2)) / 3;
no_vazio = (cx > th + 1e-9) & (cx < Wdom - th - 1e-9) & (cy > th + 1e-9);
tri = tri_all(~no_vazio, :);
fprintf('   %d triangulos apos filtro\n', size(tri,1));

%% Contorno do perfil U
xe = [0, Wdom, Wdom, Wdom-th, Wdom-th, th, th, 0, 0];
ye = [0, 0,    H,    H,       th,      th, H,  H, 0];

%% Colormap RdBu simetrico
nc = 256; half = nc/2;
r1 = linspace(0.017, 1.0, half)'; g1 = linspace(0.114, 1.0, half)'; b1 = linspace(0.388, 1.0, half)';
r2 = linspace(1.0, 0.403, half)'; g2 = linspace(1.0, 0.000, half)'; b2 = linspace(1.0, 0.122, half)';
cmap_rdbu = [r1, g1, b1; r2, g2, b2];

%% Limites dos campos
lim_qx = max(abs(qx)) * 1.05;
lim_qm = max(q_mag)    * 1.05;
if lim_qx < 1e-10; lim_qx = 1; end
if lim_qm < 1e-10; lim_qm = 1; end

campos  = {T,       qx,        qy_plot,   q_mag};
vmins   = {min(T), -lim_qx,   -lim_qy,   min(q_mag)*0.95};
vmaxs   = {max(T),  lim_qx,    lim_qy,   lim_qm};
cmaps   = {jet(256), cmap_rdbu, cmap_rdbu, hot(256)};
titulos = {'TEMPERATURE T', 'HEAT FLUX q_x', ...
           'HEAT FLUX q_y', 'HEAT FLUX MAGNITUDE |q|'};

fprintf('6) Gerando figura unica com 4 subplots...\n');

%% FIGURA UNICA COM 4 SUBPLOTS (2x2)
figure('Position', [100 100 1200 900], 'Name', 'Perfil U - Analise Termica v6');

for kfig = 1:4
    subplot(2, 2, kfig);
    hs = trisurf(tri, pts(:,1), pts(:,2), zeros(N,1), campos{kfig}, ...
                 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    colormap(gca, cmaps{kfig});

    vlo = vmins{kfig};
    vhi = vmaxs{kfig};
    if isnan(vlo) || isnan(vhi) || vlo >= vhi
        vlo = -1; vhi = 1;
    end
    caxis([vlo, vhi]);
    colorbar;

    hold on;
    plot(xe, ye, 'k-', 'LineWidth', 2);
    hold off;

    title(titulos{kfig}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x', 'FontSize', 10);
    ylabel('y', 'FontSize', 10);
    axis equal;
    xlim([-0.2 Wdom+0.2]);
    ylim([-0.2 H+0.2]);
    grid on; box on;
    set(gca, 'FontSize', 9);
end

sgtitle(sprintf('Perfil U - MFDM/MWLS  [qx analitico = %.2f W/m2]', qx_analitico), ...
        'FontSize', 13, 'FontWeight', 'bold');

%% FIGURA 2: Malha
figure('Position', [100 100 800 600], 'Name', 'Perfil U - Malha');
triplot(tri, pts(:,1), pts(:,2), 'b');
hold on;
bc_idx  = find(bc == 1);
int_idx = find(bc == 0);
plot(pts(bc_idx,1),  pts(bc_idx,2),  'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'DisplayName', 'Contorno');
plot(pts(int_idx,1), pts(int_idx,2), 'b.', 'MarkerSize', 3, 'DisplayName', 'Interior');
plot(xe, ye, 'k-', 'LineWidth', 2, 'DisplayName', 'Perfil U');
legend('Location', 'northeast', 'FontSize', 9);
hold off;
axis equal; grid on; box on;
title(sprintf('Malha Delaunay - Perfil U (%d nos)', N), ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('x [m]'); ylabel('y [m]');

fprintf('\nPronto! Figuras geradas.\n');
