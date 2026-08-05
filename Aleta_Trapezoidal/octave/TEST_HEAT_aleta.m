%TEST_HEAT_ALETA  Aleta trapezoidal com furos circulares -- MFDM/MWLS
%
%   Rode este script no Octave com todos os .m e o .msh na mesma pasta:
%       >> TEST_HEAT_aleta
%
%   O script executa DOIS casos:
%     CASO A -- VERIFICACAO (teste de aceitacao)
%        Dirichlet linear em TODO o contorno, inclusive nos furos.
%        Como funcoes lineares sao harmonicas, a solucao exata e conhecida
%        e os furos nao perturbam o campo. O erro deve cair a precisao de
%        maquina (~1e-10). Se NAO cair, ha erro de porte -- pare e investigue
%        antes de olhar o Caso B.
%     CASO B -- FISICO
%        Contorno externo Dirichlet linear + furos ISOLADOS (Neumann).
%        Produz a concentracao de fluxo na borda dos furos.
%
%   Referencia (obtida com a implementacao em Python, malha aleta.msh):
%       Caso A : erro L2 = 1,47e-10  |  qx medio = -83,3333 W/m2
%       Caso B : K = |q|max/|q|nom = 1,79   pico a 90 graus do eixo do fluxo

clear; close all; clc;

%% ---------------- parametros ----------------
kcond = 1.0;      % condutividade termica [W/(m.C)]
qvol  = 0.0;      % geracao volumetrica [W/m3]
m     = 16;       % vizinhos na estrela MWLS
g     = 0.1;      % parametro de suavizacao dos pesos

r_furo = 0.45;
furos  = [-1.5, 1.75, r_furo;
           1.5, 1.75, r_furo];

%% ---------------- malha ----------------
arq = 'aleta.msh';
if exist(arq, 'file') ~= 2
  error('Coloque %s na mesma pasta deste script.', arq);
end
[N, X, Y, BC, nt, Ttri, Elin] = msh2data(arq);

% checagem de sanidade da geometria
fprintf('\nx em [%.3f, %.3f]   y em [%.3f, %.3f]\n', ...
        min(X), max(X), min(Y), max(Y));

nrm = normais_contorno(X, Y, Elin, Ttri, N);

xmin = min(X);  xmax = max(X);
Tref = @(x) 500*(x - xmin)/(xmax - xmin);      % campo linear de referencia

%% =====================================================================
%  CASO A -- VERIFICACAO
%  =====================================================================
disp(' ');
disp('===== CASO A: VERIFICACAO (Dirichlet linear em todo o contorno) =====');

BCtipo = zeros(N,1);
BCtipo(BC ~= 0) = 1;                 % todo o contorno e Dirichlet
BCval  = zeros(N,1);
BCval(BC ~= 0) = Tref(X(BC ~= 0));

tic;
TA = localMFDM_heat_aleta(N, X, Y, BCtipo, BCval, kcond, qvol, m, g, furos, nrm);
tA = toc;

[qxA, qyA, qmA] = fluxos_aleta(N, X, Y, TA, kcond, m, g, furos);

Texato = Tref(X);
errL2  = norm(TA - Texato)/norm(Texato);
qx_an  = -500/(xmax - xmin);

fprintf('Tmin = %.4f   Tmax = %.4f       (exato: 0 / 500)\n', min(TA), max(TA));
fprintf('erro L2 relativo = %.3e         (esperado ~1e-10)\n', errL2);
fprintf('qx medio = %.5f  | analitico = %.5f W/m2\n', mean(qxA), qx_an);
fprintf('erro relativo em qx = %.4f%%\n', abs(mean(qxA)-qx_an)/abs(qx_an)*100);
fprintf('tempo = %.2f s\n', tA);

if errL2 < 1e-6
  disp('>>> TESTE DE ACEITACAO: APROVADO');
else
  disp('>>> TESTE DE ACEITACAO: REPROVADO -- ha erro na implementacao.');
  disp('    Verifique mwls.m (normalizacao) e msh2data.m (codigos de BC).');
end

%% =====================================================================
%  CASO B -- FISICO (furos isolados)
%  =====================================================================
disp(' ');
disp('===== CASO B: FUROS ISOLADOS (Neumann na fronteira interna) =====');

BCtipo = zeros(N,1);
BCtipo(BC == 1 | BC == 2 | BC == 3) = 1;      % topo, base, lados -> Dirichlet
BCtipo(BC == 4) = 2;                          % furos -> Neumann
BCval  = zeros(N,1);
sel = (BCtipo == 1);
BCval(sel) = Tref(X(sel));

tic;
TB = localMFDM_heat_aleta(N, X, Y, BCtipo, BCval, kcond, qvol, m, g, furos, nrm);
tB = toc;

[qxB, qyB, qmB] = fluxos_aleta(N, X, Y, TB, kcond, m, g, furos);

q_nom = 500/(xmax - xmin);
ifuro = find(BC == 4);
[qmax, ip] = max(qmB(ifuro));
K = qmax/q_nom;

ino = ifuro(ip);
xc  = furos(1 + (X(ino) > 0), 1);
ang = atan2(Y(ino) - 1.75, X(ino) - xc) * 180/pi;

fprintf('|q| nominal    = %.4f W/m2\n', q_nom);
fprintf('|q| max borda  = %.4f W/m2\n', qmax);
fprintf('K = |q|max/|q|nom = %.4f      (teoria furo isolado: K = 2)\n', K);
fprintf('pico a %.1f graus do eixo do fluxo  (esperado +-90)\n', ang);
fprintf('tempo = %.2f s\n', tB);

%% ---------------- graficos ----------------
figure(1); clf; set(gcf, 'Color', 'w');
campos = {TB, qxB, qyB, qmB};
titulos = {'(a) Temperatura T [C]', '(b) Fluxo q_x [W/m^2]', ...
           '(c) Fluxo q_y [W/m^2]', '(d) Magnitude |q| [W/m^2]'};
for kf = 1:4
  subplot(2,2,kf);
  h = trisurf(Ttri, X, Y, zeros(N,1), campos{kf});
  set(h, 'EdgeColor', 'none');
  shading interp;  view(2);
  axis equal tight;  colorbar;
  title(titulos{kf}, 'FontSize', 10, 'FontWeight', 'bold');
  xlabel('x [m]');  ylabel('y [m]');
end
colormap jet;

figure(2); clf; set(gcf, 'Color', 'w');
triplot(Ttri, X, Y, 'Color', [.4 .5 .8]);  hold on;
plot(X(BC ~= 0), Y(BC ~= 0), 'r.', 'MarkerSize', 8);
axis equal;  grid on;
title(sprintf('Malha -- aleta trapezoidal (%d nos, %d triangulos)', N, nt));
xlabel('x [m]');  ylabel('y [m]');

disp(' ');
disp('Concluido. Figuras 1 e 2 geradas.');
