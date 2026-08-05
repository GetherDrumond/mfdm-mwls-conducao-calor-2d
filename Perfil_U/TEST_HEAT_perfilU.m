%% TEST_HEAT_perfilU.m  (v3 — melhorias de precisao e visualizacao)
%  Conducao de calor 2D estacionaria - Perfil U (canal)
%  MFDM/MWLS — mesmo metodo do Perfil I (trilho)
%  Le malha do perfilUDATA.txt (gerado pelo Gmsh)
%
%  Melhorias v3 (orientador):
%    - m=16, sparse, g=0.1 (via localMFDM_heat e mwls atualizados)
%    - Tratamento de flutuacao de ponto flutuante nos graficos
%    - Escalas fixas baseadas nos valores analiticos esperados

clear; close all; clc;
fprintf('=== Perfil U - MFDM/MWLS v3 (melhorias) ===\n\n');

%% Parametros fisicos
kcond = 1.0;
qvol  = 0.0;

%% Parametros geometricos
Wdom = 4.0;
H    = 3.0;
th   = 0.6;

%% Ler perfilUDATA.txt
fname = 'perfilUDATA.txt';
fid = fopen(fname, 'rt');
if fid < 0
    error('Arquivo nao encontrado: %s\nColoque o arquivo na pasta atual.', fname);
end

N = fscanf(fid, '%d', [1 1]);
fprintf('Nos: %d\n', N);
Z = fscanf(fid, '%f', [3*N 1]);
X     = Z(1:3:end);
Y     = Z(2:3:end);
BCraw = Z(3:3:end);

nt = fscanf(fid, '%d', [1 1]);
fprintf('Triangulos: %d\n', nt);
Zt   = fscanf(fid, '%d', [3*nt 1]);
Ttri = [Zt(1:3:end)'; Zt(2:3:end)'; Zt(3:3:end)']';
Ttri = double(Ttri);
fclose(fid);

if any(Ttri(:) < 1) || any(Ttri(:) > N)
    error('Indices de triangulos fora do intervalo [1..%d].', N);
end

%% Condicao de contorno: T = 500*(x-xmin)/(xmax-xmin)
BCtype = zeros(N, 1);
BCtype(BCraw ~= 0) = 1;

xmin = min(X);
xmax = max(X);
BCval = zeros(N, 1);
BCval(BCraw ~= 0) = 500 * (X(BCraw ~= 0) - xmin) / (xmax - xmin);

fprintf('BC: T = 500*(x - %.2f)/(%.2f - %.2f)\n', xmin, xmax, xmin);
fprintf('Nos de contorno: %d   Nos interiores: %d\n\n', sum(BCtype), sum(BCtype==0));

%% Resolver conducao de calor (MFDM)
fprintf('Resolvendo sistema MFDM...\n');
Tnod = localMFDM_heat(N, X, Y, BCtype, BCval, kcond, qvol);
fprintf('   Tmin = %.2f C    Tmax = %.2f C\n\n', min(Tnod), max(Tnod));

%% Pos-processamento: fluxos via MWLS
fprintf('Calculando fluxos de calor (MWLS)...\n');
[FluxN, FluxT, Res] = postprocessingHeatMFDM(N, X, Y, Ttri, nt, Tnod, kcond, qvol);
fprintf('   |q|max = %.4f W/m2\n\n', max(FluxN(:,3)));

zT  = Tnod(:);
zqx = FluxN(:,1);
zqy = FluxN(:,2);
zq  = FluxN(:,3);

% ============================================================
% MELHORIA: Tratamento de flutuacao de ponto flutuante
% ============================================================
qx_analitico = -kcond * 500 / (xmax - xmin);
fprintf('Valor analitico esperado: qx = %.4f W/m2\n', qx_analitico);
fprintf('Range numerico obtido:    qx = [%.4f, %.4f]\n', min(zqx), max(zqx));
fprintf('Range numerico obtido:    qy = [%.6e, %.6e]\n', min(zqy), max(zqy));

% Threshold: se qy < 1% de |qx_analitico|, tratar como zero
thr_qy = max(abs(zqy));
if thr_qy < abs(qx_analitico) * 0.01
    fprintf('qy tratado como zero (flutuacao < 1%% de |qx_analitico|)\n');
    zqy_plot = zeros(size(zqy));
    lim_qy = abs(qx_analitico) * 0.05;
else
    zqy_plot = zqy;
    lim_qy = max(abs(zqy)) * 1.05;
end

lim_qx = max(abs(zqx)) * 1.05;
lim_qm = max(zq)        * 1.05;
if lim_qx < 1e-10; lim_qx = 1; end
if lim_qm < 1e-10; lim_qm = 1; end

%% Colormap RdBu simetrico
nc = 256; half = nc/2;
r1 = linspace(0.017, 1.0, half)'; g1 = linspace(0.114, 1.0, half)'; b1 = linspace(0.388, 1.0, half)';
r2 = linspace(1.0, 0.403, half)'; g2 = linspace(1.0, 0.000, half)'; b2 = linspace(1.0, 0.122, half)';
cmap_rdbu = [r1, g1, b1; r2, g2, b2];

%% Contorno do Perfil U
xe = [0, Wdom, Wdom, Wdom-th, Wdom-th, th, th, 0, 0];
ye = [0, 0,    H,    H,       th,      th, H,  H, 0];

%% Campos a plotar
campos  = {zT,      zqx,       zqy_plot,  zq};
vmins   = {min(zT), -lim_qx,  -lim_qy,   min(zq)*0.95};
vmaxs   = {max(zT),  lim_qx,   lim_qy,   lim_qm};
cmaps   = {jet(256), cmap_rdbu, cmap_rdbu, hot(256)};
titulos = {'TEMPERATURE T', 'HEAT FLUX q_x', ...
           'HEAT FLUX q_y', 'HEAT FLUX MAGNITUDE |q|'};

%% FIGURA 1: 4 campos em 2x2
fprintf('Gerando figura...\n');
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [50 50 1280 920], 'Name', 'Perfil U - Analise Termica v3');

for kfig = 1:4
    subplot(2, 2, kfig);

    hs = trisurf(Ttri, X, Y, zeros(N,1), campos{kfig}, ...
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
    plot(xe, ye, 'k-', 'LineWidth', 2.0);
    hold off;

    title(titulos{kfig}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x', 'FontSize', 10);
    ylabel('y', 'FontSize', 10);
    axis equal tight;
    xlim([-0.2 Wdom+0.2]); ylim([-0.2 H+0.2]);
    set(gca, 'FontSize', 9, 'Box', 'on');
    grid on;
end

sgtitle(sprintf('Perfil U - MFDM/MWLS  [qx analitico = %.2f W/m2]', qx_analitico), ...
        'FontSize', 13, 'FontWeight', 'bold');

%% FIGURA 2: Malha triangular
figure(2); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 800 600], 'Name', 'Perfil U - Malha');
triplot(Ttri, X, Y, 'b');
hold on;
bc_idx  = find(BCtype == 1);
int_idx = find(BCtype == 0);
plot(X(bc_idx),  Y(bc_idx),  'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'DisplayName', 'Contorno');
plot(X(int_idx), Y(int_idx), 'b.', 'MarkerSize', 3, 'DisplayName', 'Interior');
plot(xe, ye, 'k-', 'LineWidth', 2, 'DisplayName', 'Perfil U');
legend('Location', 'northeast', 'FontSize', 9);
hold off;
axis equal; grid on; box on;
title(sprintf('Malha - Perfil U (%d nos, %d triangulos)', N, nt), ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('x [m]'); ylabel('y [m]');

fprintf('\n=== CONCLUIDO ===\n');
fprintf('Nos: %d  |  Triangulos: %d\n', N, nt);
fprintf('Tmin=%.2f C  Tmax=%.2f C\n', min(zT), max(zT));
fprintf('|q|max=%.4f W/m2  |q|medio=%.4f W/m2\n', max(zq), mean(zq));
