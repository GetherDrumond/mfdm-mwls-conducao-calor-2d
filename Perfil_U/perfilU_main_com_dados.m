%% perfilU_main_com_dados.m  (v2 corrigida)
%  Conducao de calor 2D - Perfil U - MFDM/MWLS
%  Usando dados do arquivo perfilUDATA.txt
%  TODOS OS GRAFICOS EM UMA UNICA FIGURA (layout 2x2)
%
%  Correcoes:
%    - Visualizacao com trisurf + view(2) (sem griddata na cavidade)
%    - Variavel Wdom para evitar conflito
%    - Protecao de limites de caxis contra zero

clear; clc; close all;
fprintf('=== Perfil U - MFDM/MWLS com dados do arquivo ===\n\n');

Wdom = 4.0; H = 3.0; th = 0.6;
kcond = 1.0; qvol = 0.0;

%% 1) Carregar perfilUDATA.txt  (formato: N, nos x y bc, nt, triangulos)
fprintf('1) Carregando perfilUDATA.txt...\n');
fid = fopen('perfilUDATA.txt', 'r');
if fid == -1, error('Arquivo perfilUDATA.txt nao encontrado!'); end

N  = fscanf(fid, '%d', 1);
Z  = fscanf(fid, '%f', [3*N 1]);
X     = Z(1:3:end);
Y     = Z(2:3:end);
BCraw = Z(3:3:end);

nt = fscanf(fid, '%d', 1);
Zt   = fscanf(fid, '%d', [3*nt 1]);
Ttri = [Zt(1:3:end)'; Zt(2:3:end)'; Zt(3:3:end)']';
Ttri = double(Ttri);
fclose(fid);

fprintf('   %d nos  (contorno: %d, interior: %d)\n', N, sum(BCraw~=0), sum(BCraw==0));
fprintf('   %d triangulos\n', nt);

%% 2) Condicao de contorno: T = 500*(x-xmin)/(xmax-xmin)
BCtype = zeros(N,1);
BCtype(BCraw ~= 0) = 1;
xmin = min(X); xmax = max(X);
BCval = zeros(N,1);
BCval(BCraw ~= 0) = 500 * (X(BCraw ~= 0) - xmin) / (xmax - xmin);

%% 3) Resolver MFDM
fprintf('2) Resolvendo MFDM...\n');
Tnod = localMFDM_heat(N, X, Y, BCtype, BCval, kcond, qvol);
fprintf('   Tmin=%.2f C  Tmax=%.2f C\n', min(Tnod), max(Tnod));

%% 4) Calcular fluxos via postprocessingHeatMFDM (mesmo metodo do Perfil I)
fprintf('3) Calculando fluxos (MWLS)...\n');
[FluxN, ~, ~] = postprocessingHeatMFDM(N, X, Y, Ttri, nt, Tnod, kcond, qvol);
zT  = Tnod(:);
zqx = FluxN(:,1);
zqy = FluxN(:,2);
zq  = FluxN(:,3);
fprintf('   |q|max = %.4f W/m2\n', max(zq));

%% 5) Contorno do perfil U
xe = [0, Wdom, Wdom, Wdom-th, Wdom-th, th, th, 0, 0];
ye = [0, 0,    H,    H,       th,      th, H,  H, 0];

%% 6) Colormap RdBu manual
nc = 256; half = nc/2;
r1 = linspace(0.017, 1.0, half)'; g1 = linspace(0.114, 1.0, half)'; b1 = linspace(0.388, 1.0, half)';
r2 = linspace(1.0, 0.403, half)'; g2 = linspace(1.0, 0.000, half)'; b2 = linspace(1.0, 0.122, half)';
cmap_rdbu = [r1, g1, b1; r2, g2, b2];

%% 7) FIGURA UNICA 2x2 — usando trisurf + view(2)
fprintf('4) Gerando figura 2x2...\n');
lim_qx = max(abs(zqx)) * 1.05;
lim_qy = max(abs(zqy)) * 1.05;
lim_qm = max(zq) * 1.05;

if lim_qx < 1e-10; lim_qx = 1; end
if lim_qy < 1e-10; lim_qy = 1; end
if lim_qm < 1e-10; lim_qm = 1; end

campos  = {zT,      zqx,        zqy,        zq};
vmins   = {min(zT), -lim_qx,   -lim_qy,    0};
vmaxs   = {max(zT),  lim_qx,    lim_qy,    lim_qm};
cmaps   = {jet(256), cmap_rdbu,  cmap_rdbu,  hot(256)};
titulos = {'(a) TEMPERATURA T', '(b) FLUXO DE CALOR q_x', ...
           '(c) FLUXO DE CALOR q_y', '(d) MAGNITUDE DO FLUXO |q|'};
cunits  = {'T [graus C]', 'q_x [W/m2]', 'q_y [W/m2]', '|q| [W/m2]'};

figure('Position', [50 50 1400 1000], 'Name', 'Perfil U - Analise Termica', 'Color', 'white');

for p = 1:4
    subplot(2, 2, p);

    % FIX: trisurf com Z=0, CData=campo, view(2)
    hs = trisurf(Ttri, X, Y, zeros(N,1), campos{p}, ...
                 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);

    colormap(gca, cmaps{p});
    caxis([vmins{p}, vmaxs{p}]);
    cb = colorbar;
    ylabel(cb, cunits{p}, 'FontSize', 11, 'FontWeight', 'bold');

    hold on;
    plot(xe, ye, 'k-', 'LineWidth', 2.5);
    if p == 1
        ci = find(BCtype == 1);
        plot(X(ci), Y(ci), 'k.', 'MarkerSize', 6);
    end
    hold off;

    title(titulos{p}, 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('x [m]', 'FontSize', 11);
    ylabel('y [m]', 'FontSize', 11);
    axis equal tight;
    xlim([-0.2 Wdom+0.2]); ylim([-0.2 H+0.2]);
    set(gca, 'FontSize', 10, 'Box', 'on', 'LineWidth', 1.2);
    grid on;
end

sgtitle('Perfil U - Analise de Transferencia de Calor (MFDM/MWLS)', ...
        'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n=== CONCLUIDO ===\n');
fprintf('Nos: %d  |  Triangulos: %d\n', N, nt);
fprintf('Tmin=%.2f C  Tmax=%.2f C\n', min(zT), max(zT));
fprintf('|q|max=%.4f W/m2  |q|medio=%.4f W/m2\n', max(zq), mean(zq));
