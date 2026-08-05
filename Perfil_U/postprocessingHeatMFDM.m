function [FluxN,FluxT,Res] = postprocessingHeatMFDM(N,X,Y,Ttri,nt,T,k,q)
%POSTPROCESSINGHEATMFDM - Postprocessing for MFDM steady heat conduction
%
%   [FluxN,FluxT,Res] = postprocessingHeatMFDM(N,X,Y,Ttri,nt,T,k,q)
%
%   Outputs
%     FluxN (N x 3): nodal heat flux [qx, qy, |q|]
%     FluxT (nt x 3): triangle-integrated flux
%     Res (nt x 1): residual indicator

    m = 16;            % MELHORIA: m=16 (consistente com localMFDM_heat)
    if m > N
        m = N;
    end

    if ~isscalar(q)
        error('For residual computation, q is expected to be a scalar.');
    end

    % ---- Nodal flux recovery ----
    FluxN = zeros(N,3);
    for i = 1:N
        S = star(N, X(i), Y(i), X, Y, m);
        M = mwls(X(i), Y(i), X, Y, S, m, 1);  % p=1, g=0.1 (default)
        Tvals = T(S(1:m,1));

        if any(isnan(M(:))) || any(isnan(Tvals))
            FluxN(i,:) = [0, 0, 0];
            continue;
        end

        Tx = M(2,:) * Tvals;
        Ty = M(3,:) * Tvals;
        qx = -k * Tx;
        qy = -k * Ty;
        FluxN(i,:) = [qx, qy, sqrt(qx^2 + qy^2)];
    end

    % ---- Triangle-integrated flux + residual indicator ----
    nG = 3;
    Gw = [1/3 1/3 1/3];
    FluxT = zeros(nt,3);
    Res   = zeros(nt,1);

    for ktri = 1:nt
        w  = Ttri(ktri,:)';
        Xt = X(w);
        Yt = Y(w);

        Gpx = [0.5*(Xt(1)+Xt(2)) 0.5*(Xt(2)+Xt(3)) 0.5*(Xt(1)+Xt(3))];
        Gpy = [0.5*(Yt(1)+Yt(2)) 0.5*(Yt(2)+Yt(3)) 0.5*(Yt(1)+Yt(3))];

        wek = [Xt(1)-Xt(2), Yt(1)-Yt(2);
               Xt(2)-Xt(3), Yt(2)-Yt(3);
               Xt(3)-Xt(1), Yt(3)-Yt(1)];
        dl = [sqrt(wek(1,1)^2 + wek(1,2)^2),
              sqrt(wek(2,1)^2 + wek(2,2)^2),
              sqrt(wek(3,1)^2 + wek(3,2)^2)];
        p = 0.5*sum(dl);
        J = sqrt(p*(p-dl(1))*(p-dl(2))*(p-dl(3)));

        int_qx = 0;
        int_qy = 0;
        res    = 0;

        for pg = 1:nG
            xg = Gpx(pg);
            yg = Gpy(pg);
            S  = star(N, xg, yg, X, Y, m);
            M  = mwls(xg, yg, X, Y, S, m, 2);  % p=2, g=0.1 (default)
            Tvals = T(S(1:m,1));

            if any(isnan(M(:))) || any(isnan(Tvals))
                continue;
            end

            Tder = M * Tvals;

            Tx  = Tder(2);
            Ty  = Tder(3);
            Txx = Tder(4);
            Tyy = Tder(6);

            int_qx = int_qx + J*Gw(pg) * (-k*Tx);
            int_qy = int_qy + J*Gw(pg) * (-k*Ty);

            res = res + abs((Txx + Tyy) + q/k);
        end

        FluxT(ktri,1) = int_qx;
        FluxT(ktri,2) = int_qy;
        FluxT(ktri,3) = sqrt(int_qx^2 + int_qy^2);
        Res(ktri)     = res;
    end
end
