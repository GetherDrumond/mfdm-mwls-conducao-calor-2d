function T = localMFDM_heat(N,X,Y,BCtype,BCval,k,q)
%LOCALMFDM_HEAT - MFDM solver for steady 2D heat conduction (constant k)
%
%   Governing equation (steady state, isotropic, constant conductivity):
%       nabla^2 T = - q/k
%
%   T = localMFDM_heat(N,X,Y,BCtype,BCval,k,q) returns nodal temperatures T.
%
%   Inputs
%     N        - number of nodes
%     X,Y      - nodal coordinates (N x 1)
%     BCtype   - boundary code per node:
%                  0 = internal node (PDE collocation)
%                  1 = Dirichlet boundary node (prescribed temperature)
%     BCval    - prescribed temperature for BCtype=1 (N x 1)
%     k        - thermal conductivity (scalar, constant)
%     q        - volumetric heat generation (scalar or N x 1)

    A = sparse(N,N);   % MELHORIA: sparse em vez de zeros (eficiencia + precisao)
    B = zeros(N,1);

    m = 16;            % MELHORIA: m=16 em vez de 9 (mais vizinhos = mais preciso)
    if N < m
        m = N;
    end

    if isscalar(q)
        qv = q*ones(N,1);
    else
        qv = q(:);
        if numel(qv) ~= N
            error('q must be scalar or a vector of length N');
        end
    end

    for i = 1:N
        if BCtype(i) == 0
            % Internal node: Laplacian(T) = -q/k
            S = star(N, X(i), Y(i), X, Y, m);
            M = mwls(X(i), Y(i), X, Y, S, m);

            row = M(4,:) + M(6,:);    % T_xx + T_yy
            if any(isnan(row)) || any(isinf(row))
                % Fallback: tentar com mais vizinhos
                m2 = min(m + 8, N);
                S2 = star(N, X(i), Y(i), X, Y, m2);
                M2 = mwls(X(i), Y(i), X, Y, S2, m2);
                row2 = M2(4,:) + M2(6,:);
                if ~any(isnan(row2)) && ~any(isinf(row2))
                    A(i, S2(1:m2,1)) = row2;
                else
                    A(i,i) = 1;
                end
            else
                A(i, S(1:m,1)) = row;
            end
            B(i) = -qv(i)/k;
        elseif BCtype(i) == 1
            % Dirichlet boundary node: T = prescribed
            A(i,i) = 1;
            B(i)   = BCval(i);
        else
            error('Unsupported BCtype=%d at node %d.', BCtype(i), i);
        end
    end

    T = A \ B;

    % Tratar NaN residuais por interpolacao de vizinhos
    nan_idx = find(isnan(T));
    if ~isempty(nan_idx)
        warning('localMFDM_heat: %d nos com NaN detectados, interpolando...', length(nan_idx));
        for ii = 1:length(nan_idx)
            ni = nan_idx(ii);
            S = star(N, X(ni), Y(ni), X, Y, m);
            viz = S(2:m, 1);
            valid = viz(~isnan(T(viz)));
            if ~isempty(valid)
                T(ni) = mean(T(valid));
            else
                T(ni) = 0;
            end
        end
    end
end
