%% plotar_4campos.m  (v2 corrigida)
%  Plota T, qx, qy, |q| em uma unica janela (2x2)
%  Usa trisurf + view(2) para visualizacao 2D correta
%  (sem griddata que interpolava na cavidade vazia)
%
%  Executar APOS rodar perfilU_main (variaveis T, qx, qy, q_mag, tri, pts, N ja existem)

if ~exist('T','var') || ~exist('tri','var')
  error('Execute perfilU_main primeiro para gerar T, qx, qy, q_mag, tri, pts.');
end

%% Parametros
Wdom = 4.0; H = 3.0; th = 0.6;
xe = [0,Wdom,Wdom,Wdom-th,Wdom-th,th,th,0,0];
ye = [0,0,   H,   H,      th,     th,H, H,0];

%% Colormap RdBu
nc = 256; half = nc/2;
r1 = linspace(0.017,1.0,half)'; g1 = linspace(0.114,1.0,half)'; b1 = linspace(0.388,1.0,half)';
r2 = linspace(1.0,0.403,half)'; g2 = linspace(1.0,0.000,half)'; b2 = linspace(1.0,0.122,half)';
cmap_rdbu = [r1,g1,b1; r2,g2,b2];

%% Limites de cada campo (com protecao)
lim_qx   = max(abs(qx))  * 1.02;
lim_qy   = max(abs(qy))  * 1.02;
lim_qmag = max(q_mag)    * 1.02;

if lim_qx   < 1e-10; lim_qx   = 1; end
if lim_qy   < 1e-10; lim_qy   = 1; end
if lim_qmag < 1e-10; lim_qmag = 1; end

campos  = {T,         qx,       qy,       q_mag};
vmins   = {min(T),   -lim_qx,  -lim_qy,  0};
vmaxs   = {max(T),    lim_qx,   lim_qy,  lim_qmag};
cmaps   = {jet(256),  cmap_rdbu, cmap_rdbu, hot(256)};
titulos = {'Temperatura  T  [graus C]', ...
           'Fluxo de Calor  qx  [W/m2]', ...
           'Fluxo de Calor  qy  [W/m2]', ...
           'Magnitude do Fluxo  |q|  [W/m2]'};
cunits  = {'T [graus C]', 'q_x [W/m2]', 'q_y [W/m2]', '|q| [W/m2]'};

%% Figura principal
figure('Position', [30 30 1300 950], 'Name', 'Perfil U - Campos Termicos');

N = size(pts, 1);

for p = 1:4
  subplot(2, 2, p);

  % FIX: trisurf com Z=0, CData=campo, view(2)
  hs = trisurf(tri, pts(:,1), pts(:,2), zeros(N,1), campos{p}, ...
               'EdgeColor', 'none', 'FaceColor', 'interp');
  view(2);

  hold on;
  plot(xe, ye, 'k-', 'LineWidth', 2.0);
  hold off;

  axis equal tight;
  xlim([-0.2, Wdom+0.2]); ylim([-0.2, H+0.2]);

  title(titulos{p}, 'FontSize', 11, 'FontWeight', 'bold');
  xlabel('x [m]', 'FontSize', 10);
  ylabel('y [m]', 'FontSize', 10);
  set(gca, 'FontSize', 9);

  colormap(gca, cmaps{p});
  caxis([vmins{p}, vmaxs{p}]);
  cb = colorbar;
  tks = linspace(vmins{p}, vmaxs{p}, 6);
  set(cb, 'Ticks', tks, ...
      'TickLabels', arrayfun(@(v)sprintf('%.1f',v), tks, 'UniformOutput', false));
  ylabel(cb, cunits{p}, 'FontSize', 9);

  grid on; box on;
end

sgtitle('MFDM - Conducao de Calor 2D  |  Perfil U  |  BC: T linear em x (0-500 graus C)', ...
        'FontSize', 12, 'FontWeight', 'bold');
