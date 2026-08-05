%% perfilU_gmsh_main.m  (v2 corrigida)
%  Versao que usa a malha gerada pelo Gmsh (perfilUDATA.txt)
%  em vez de gerar os pontos internamente.
%
%  Correcoes:
%    - Visualizacao com trisurf + view(2) (sem griddata na cavidade)
%    - Variavel Wdom para evitar conflito
%    - Protecao de limites de caxis contra zero
%    - 4 figuras separadas com colormaps corretos
%
%  Uso: coloque perfilUDATA.txt na mesma pasta e execute este script.

clear; clc; close all;
fprintf('=== Perfil U - MFDM com malha Gmsh ===\n\n');

Wdom = 4.0;
H    = 3.0;
th   = 0.6;
k_cond = 1.0;

%% 1. Carregar malha do Gmsh
fprintf('1) Carregando malha do Gmsh...\n');
[pts, bc] = carregar_malha_gmsh('perfilUDATA.txt');
N = size(pts, 1);

%% 2. Condicao de contorno: T = 500*(x-xmin)/(xmax-xmin)
xmin = min(pts(:,1));
xmax = max(pts(:,1));
T_bc = 500.0 * (pts(:,1) - xmin) / (xmax - xmin);

%% 3. Montar e resolver sistema MFDM
fprintf('2) Montando sistema MFDM (%d nos)...\n', N);
[A, rhs] = montar_sistema_mfdm(pts, bc, T_bc);

fprintf('3) Resolvendo...\n');
T = A \ rhs;
fprintf('   Tmin = %.2f C    Tmax = %.2f C\n', min(T), max(T));

%% 4. Calcular fluxos
fprintf('4) Calculando fluxos...\n');
qx = zeros(N,1); qy = zeros(N,1);
n_flux = 10;
for i = 1:N
  dx = pts(:,1)-pts(i,1); dy = pts(:,2)-pts(i,2);
  dist = sqrt(dx.^2+dy.^2); dist(i) = Inf;
  [dv,idx] = sort(dist);
  viz = idx(1:min(n_flux,N-1)); dv = dv(1:min(n_flux,N-1));
  xi = pts(viz,1)-pts(i,1); yi = pts(viz,2)-pts(i,2);
  dT = T(viz)-T(i); h = max(dv)*1.1;
  w = exp(-6*(dv/h).^2); Wm = diag(w); Pm = [xi,yi];
  try
    g = (Pm'*Wm*Pm+1e-12*eye(2))\(Pm'*Wm*dT);
    qx(i) = -k_cond*g(1); qy(i) = -k_cond*g(2);
  catch; end
end
q_mag = sqrt(qx.^2+qy.^2);

%% 5. Triangulacao + filtro do vazio
tri_all = delaunay(pts(:,1), pts(:,2));
cx = (pts(tri_all(:,1),1)+pts(tri_all(:,2),1)+pts(tri_all(:,3),1))/3;
cy = (pts(tri_all(:,1),2)+pts(tri_all(:,2),2)+pts(tri_all(:,3),2))/3;
% FIX: filtro correto da cavidade
no_vazio = (cx>th+1e-9)&(cx<Wdom-th-1e-9)&(cy>th+1e-9);
tri = tri_all(~no_vazio,:);

%% Contorno para plotagem
xe = [0,Wdom,Wdom,Wdom-th,Wdom-th,th,th,0,0];
ye = [0,0,   H,   H,      th,     th,H, H,0];

%% Colormap RdBu
nc=256; half=nc/2;
r1=linspace(0.017,1.0,half)'; g1=linspace(0.114,1.0,half)'; b1=linspace(0.388,1.0,half)';
r2=linspace(1.0,0.403,half)'; g2=linspace(1.0,0.000,half)'; b2=linspace(1.0,0.122,half)';
cmap_rdbu=[r1,g1,b1;r2,g2,b2];

%% Limites (com protecao)
lim_qx = max(abs(qx))*1.02;
lim_qy = max(abs(qy))*1.02;
lim_qm = max(q_mag)*1.02;
if lim_qx < 1e-10; lim_qx = 1; end
if lim_qy < 1e-10; lim_qy = 1; end
if lim_qm < 1e-10; lim_qm = 1; end

%% FIGURA 1 - Temperatura
figure('Position',[10 400 680 500],'Name','Temperatura T');
trisurf(tri,pts(:,1),pts(:,2),zeros(N,1),T,'EdgeColor','none','FaceColor','interp');
view(2); colormap('jet'); caxis([min(T) max(T)]);
cb=colorbar; tks=linspace(min(T),max(T),6);
set(cb,'Ticks',tks,'TickLabels',arrayfun(@(v)sprintf('%.0f',v),tks,'UniformOutput',false));
ylabel(cb,'T [graus C]','FontSize',11);
hold on; plot(xe,ye,'k-','LineWidth',1.8); hold off;
title('Temperatura  T  [graus C]','FontSize',13,'FontWeight','bold');
xlabel('x [m]','FontSize',11); ylabel('y [m]','FontSize',11);
axis equal; xlim([-0.2 Wdom+0.2]); ylim([-0.2 H+0.2]);
grid on; box on;

%% FIGURA 2 - Fluxo qx
figure('Position',[700 400 680 500],'Name','Fluxo qx');
trisurf(tri,pts(:,1),pts(:,2),zeros(N,1),qx,'EdgeColor','none','FaceColor','interp');
view(2); colormap(gca, cmap_rdbu);
caxis([-lim_qx lim_qx]);
cb=colorbar; tks=linspace(-lim_qx,lim_qx,7);
set(cb,'Ticks',tks,'TickLabels',arrayfun(@(v)sprintf('%.1f',v),tks,'UniformOutput',false));
ylabel(cb,'qx [W/m2]','FontSize',11);
hold on; plot(xe,ye,'k-','LineWidth',1.8); hold off;
title('Fluxo de Calor  qx  [W/m2]','FontSize',13,'FontWeight','bold');
xlabel('x [m]','FontSize',11); ylabel('y [m]','FontSize',11);
axis equal; xlim([-0.2 Wdom+0.2]); ylim([-0.2 H+0.2]);
grid on; box on;

%% FIGURA 3 - Fluxo qy
figure('Position',[10 30 680 500],'Name','Fluxo qy');
trisurf(tri,pts(:,1),pts(:,2),zeros(N,1),qy,'EdgeColor','none','FaceColor','interp');
view(2); colormap(gca, cmap_rdbu);
caxis([-lim_qy lim_qy]);
cb=colorbar; tks=linspace(-lim_qy,lim_qy,7);
set(cb,'Ticks',tks,'TickLabels',arrayfun(@(v)sprintf('%.1f',v),tks,'UniformOutput',false));
ylabel(cb,'qy [W/m2]','FontSize',11);
hold on; plot(xe,ye,'k-','LineWidth',1.8); hold off;
title('Fluxo de Calor  qy  [W/m2]','FontSize',13,'FontWeight','bold');
xlabel('x [m]','FontSize',11); ylabel('y [m]','FontSize',11);
axis equal; xlim([-0.2 Wdom+0.2]); ylim([-0.2 H+0.2]);
grid on; box on;

%% FIGURA 4 - Magnitude |q|
figure('Position',[700 30 680 500],'Name','Magnitude |q|');
trisurf(tri,pts(:,1),pts(:,2),zeros(N,1),q_mag,'EdgeColor','none','FaceColor','interp');
view(2); colormap('hot'); caxis([0 lim_qm]);
cb=colorbar; tks=linspace(0,max(q_mag),6);
set(cb,'Ticks',tks,'TickLabels',arrayfun(@(v)sprintf('%.1f',v),tks,'UniformOutput',false));
ylabel(cb,'|q| [W/m2]','FontSize',11);
hold on; plot(xe,ye,'k-','LineWidth',1.8); hold off;
title('Magnitude do Fluxo de Calor  |q|  [W/m2]','FontSize',13,'FontWeight','bold');
xlabel('x [m]','FontSize',11); ylabel('y [m]','FontSize',11);
axis equal; xlim([-0.2 Wdom+0.2]); ylim([-0.2 H+0.2]);
grid on; box on;

fprintf('\nPronto! 4 figuras geradas.\n');
