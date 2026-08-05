% TEST_HEAT.m — Perfil I (Railroad Rail)
% MFDM/MWLS — Conducao de calor 2D estacionaria
% UFES — Transferencia de Calor Computacional
%
% Melhorias aplicadas:
%   - m=16 (mais vizinhos), sparse, g=0.1 (regularizacao de pesos)
%   - Tratamento de flutuacao de ponto flutuante nos graficos
%   - Escalas fixas baseadas nos valores analiticos esperados

clear; close all; clc;

%% ---- Parametros fisicos ----
kcond = 1.0;   % condutividade termica k [W/(m.C)]
qvol  = 0.0;   % geracao volumetrica de calor (Laplace se = 0)

%% ---- Ler dominio irregular ----
fname = 'railroadrailDATA.txt';
fid = fopen(fname, 'rt');
if fid < 0
    error('Could not open %s. Put the file in the current folder.', fname);
end

N = fscanf(fid, '%d', [1 1]);
if isempty(N) || N <= 0
    fclose(fid);
    error('Invalid N read from %s.', fname);
end

Z = fscanf(fid, '%f', [3*N 1]);
if numel(Z) < 3*N
    fclose(fid);
    error('Not enough node data in %s.', fname);
end

X     = Z(1:3:end);
Y     = Z(2:3:end);
BCraw = Z(3:3:end);

nt = fscanf(fid, '%d', [1 1]);
if isempty(nt) || nt <= 0
    fclose(fid);
    error('Invalid nt read from %s.', fname);
end

Zt = fscanf(fid, '%d', [3*nt 1]);
if numel(Zt) < 3*nt
    fclose(fid);
    error('Not enough triangle data in %s.', fname);
end

Ttri = [Zt(1:3:end)'; Zt(2:3:end)'; Zt(3:3:end)']';
Ttri = double(Ttri);
fclose(fid);

if any(Ttri(:) < 1) || any(Ttri(:) > N)
    error('Ttri has indices outside [1..N]. Check %s.', fname);
end

% ---- BC: gradiente linear em x ----
BCtype = zeros(N, 1);
BCtype(BCraw ~= 0) = 1;

xmin = min(X);
xmax = max(X);
BCval = zeros(N, 1);
BCval(BCraw ~= 0) = 500 * (X(BCraw ~= 0) - xmin) / (xmax - xmin);

%% ---- Resolver conducao de calor ----
disp('PRIMARY NODAL SOLUTION (STEADY HEAT CONDUCTION)');
Tnod = localMFDM_heat(N, X, Y, BCtype, BCval, kcond, qvol);

%% ---- Pos-processamento: fluxos ----
disp('RECOVERED NODAL HEAT FLUXES');
[FluxN, FluxT, Res] = postprocessingHeatMFDM(N, X, Y, Ttri, nt, Tnod, kcond, qvol);

zT  = Tnod(:);
zqx = FluxN(:,1);
zqy = FluxN(:,2);
zq  = FluxN(:,3);

% ============================================================
% MELHORIA: Tratamento de flutuacao de ponto flutuante
% ============================================================
% Valor analitico esperado: qx = -k * 500/(xmax-xmin), qy = 0
% Flutuacoes numericas da ordem de 10^-12 devem ser tratadas como zero.
%
% Criterio: se a variacao relativa de um campo for < 1% do valor
% analitico esperado, considerar o campo como uniforme (sem manchas).

qx_analitico = -kcond * 500 / (xmax - xmin);
fprintf('Valor analitico esperado: qx = %.4f W/m2\n', qx_analitico);
fprintf('Range numerico obtido:    qx = [%.4f, %.4f]\n', min(zqx), max(zqx));
fprintf('Range numerico obtido:    qy = [%.6e, %.6e]\n', min(zqy), max(zqy));

% Threshold para considerar flutuacao como zero
thr_qy = max(abs(zqy));
if thr_qy < abs(qx_analitico) * 0.01
    % qy e essencialmente zero — flutuacao de ponto flutuante
    fprintf('qy tratado como zero (flutuacao < 1%% de |qx_analitico|)\n');
    zqy_plot = zeros(size(zqy));
    lim_qy = abs(qx_analitico) * 0.05;  % escala simetrica pequena para referencia
else
    zqy_plot = zqy;
    lim_qy = max(abs(zqy)) * 1.05;
end

% Para qx: usar escala baseada no range real
lim_qx = max(abs(zqx)) * 1.05;

% Para |q|: usar escala baseada no range real
lim_qm = max(zq) * 1.05;
if lim_qm < 1e-10; lim_qm = 1; end

Tmin_val = min(zT);
Tmax_val = max(zT);
if isnan(Tmin_val) || isnan(Tmax_val) || Tmin_val == Tmax_val
    Tmin_val = 0; Tmax_val = 500;
end

%% ---- Colormap RdBu simetrico ----
nc = 256; half = nc/2;
r1 = linspace(0.017, 1.0, half)'; g1 = linspace(0.114, 1.0, half)'; b1 = linspace(0.388, 1.0, half)';
r2 = linspace(1.0, 0.403, half)'; g2 = linspace(1.0, 0.000, half)'; b2 = linspace(1.0, 0.122, half)';
cmap_rdbu = [r1, g1, b1; r2, g2, b2];

%% ---- Visualizacao 2D com trisurf + view(2) ----
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [50 50 1280 920]);

% Campos a plotar (usando zqy_plot para qy tratado)
fields  = {zT, zqx, zqy_plot, zq};
titles  = {'TEMPERATURE T', 'HEAT FLUX q_x', ...
           'HEAT FLUX q_y', 'HEAT FLUX MAGNITUDE |q|'};
vmins   = {Tmin_val, -lim_qx, -lim_qy, min(zq)*0.95};
vmaxs   = {Tmax_val,  lim_qx,  lim_qy, lim_qm};
cmaps   = {jet(256), cmap_rdbu, cmap_rdbu, hot(256)};

for kfig = 1:4
    subplot(2, 2, kfig);
    hs = trisurf(Ttri, X, Y, zeros(N, 1), fields{kfig});
    set(hs, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);

    colormap(gca, cmaps{kfig});

    % Protecao: garantir que vmin < vmax
    vlo = vmins{kfig};
    vhi = vmaxs{kfig};
    if isnan(vlo) || isnan(vhi) || vlo >= vhi
        vlo = -1; vhi = 1;
    end
    caxis([vlo, vhi]);

    cb = colorbar;
    title(titles{kfig}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x', 'FontSize', 10);
    ylabel('y', 'FontSize', 10);
    axis equal tight;
    grid on; box on;
    set(gca, 'FontSize', 9);
end

sgtitle(sprintf('Perfil I (Railroad Rail) - MFDM/MWLS  [qx analitico = %.2f W/m2]', qx_analitico), ...
        'FontSize', 13, 'FontWeight', 'bold');

%% ---- Malha 2D ----
figure(2); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 800 600]);
triplot(Ttri, X, Y, 'b');
hold on;
bc_idx  = find(BCtype == 1);
int_idx = find(BCtype == 0);
plot(X(bc_idx),  Y(bc_idx),  'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'DisplayName', 'Contorno');
plot(X(int_idx), Y(int_idx), 'b.', 'MarkerSize', 3, 'DisplayName', 'Interior');
legend('Location', 'northeast', 'FontSize', 9);
hold off;
axis equal; grid on; box on;
title(sprintf('Malha - Perfil I (%d nos, %d triangulos)', N, nt), ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('x'); ylabel('y');

disp('Concluido! Verifique as figuras 1 e 2.');
