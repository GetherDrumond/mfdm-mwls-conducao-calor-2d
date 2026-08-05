function [S] = star_visivel(N, x, y, X, Y, m, furos)
%STAR_VISIVEL Estrela MFD com criterio de VISIBILIDADE.
%
%   S = star_visivel(N, x, y, X, Y, m, furos)
%   furos: matriz (nf x 3) com [xc  yc  r]. Passe [] para dominio sem furos,
%          caso em que a funcao se reduz ao star.m original.
%
%   CORRECAO (Bug 7) -- dominios multiplamente conexos.
%   O criterio de distancia pura do star.m permite que a estrela ATRAVESSE
%   um furo e capture nos do lado oposto, que nao tem relacao fisica com o
%   no central. Em furos pequenos o efeito e sistematico: para r = 0,15 m
%   obteve-se K = 5,58 (distancia pura) contra K = 1,81 (visibilidade).
%   Aqui um vizinho so e aceito se o segmento que o liga ao no central nao
%   cruzar o interior de nenhum furo.

  D = sqrt((x - X(:)').^2 + (y - Y(:)').^2);
  Sall = [1:N; D]';
  Sall = sortrows(Sall, 2);

  if isempty(furos)
    S = Sall(1:m, :);
    return;
  end

  ncand = min(N, 4*m);            % candidatos examinados
  S  = zeros(m, 2);
  na = 0;

  for c = 1:ncand
    j  = Sall(c,1);
    qx = X(j);  qy = Y(j);
    dx = qx - x;  dy = qy - y;
    L2 = dx^2 + dy^2;

    ok = true;
    if L2 >= 1e-18
      for f = 1:size(furos,1)
        xc = furos(f,1);  yc = furos(f,2);  r = furos(f,3);
        t  = ((xc - x)*dx + (yc - y)*dy) / L2;    % projecao no segmento
        t  = max(0, min(1, t));
        dc = sqrt((x + t*dx - xc)^2 + (y + t*dy - yc)^2);
        if dc < r*(1 - 1e-6)                      % cruza o interior do furo
          ok = false;
          break;
        end
      end
    end

    if ok
      na = na + 1;
      S(na,:) = Sall(c,:);
      if na == m, break; end
    end
  end

  % ---- fallback: se a visibilidade nao completou a estrela, preenche
  %      com os mais proximos ainda nao usados (evita estrela deficiente)
  if na < m
    warning('star_visivel: apenas %d vizinhos visiveis em (%.3f, %.3f)', na, x, y);
    for c = 1:size(Sall,1)
      if na == m, break; end
      if ~any(S(1:na,1) == Sall(c,1))
        na = na + 1;
        S(na,:) = Sall(c,:);
      end
    end
  end
end
