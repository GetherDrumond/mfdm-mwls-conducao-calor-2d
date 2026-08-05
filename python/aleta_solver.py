"""
Aleta trapezoidal com furos - MFDM/MWLS
Reproduz a metodologia do relatorio: m=16, p=2, g=0.1, Tikhonov 1e-12.
Itens: (i) simulacao (ii) convergencia (iii) raio dos furos (iv) MFDM vs FEM
"""
import numpy as np, time, json
from scipy.spatial import cKDTree
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

# =====================================================================
# LEITURA DA MALHA
# =====================================================================
def read_msh(fn):
    with open(fn) as f:
        lines = f.read().split('\n')
    i = 0; Ttri = []; Tlin = []
    while i < len(lines):
        L = lines[i].strip()
        if L == '$Nodes':
            i += 1; N = int(lines[i])
            X = np.zeros(N+1); Y = np.zeros(N+1); BC = np.zeros(N+1, int)
            for _ in range(N):
                i += 1; v = lines[i].split(); k = int(v[0])
                X[k] = float(v[1]); Y[k] = float(v[2])
        elif L == '$Elements':
            i += 1; ne = int(lines[i])
            for _ in range(ne):
                i += 1; v = [int(float(t)) for t in lines[i].split()]
                et, ntg, ph = v[1], v[2], v[3]; conn = v[3+ntg:]
                if et == 1:
                    Tlin.append((ph, conn[0], conn[1]))
                    for n in conn:
                        BC[n] = ph if BC[n] == 0 else min(BC[n], ph)
                elif et == 2:
                    Ttri.append(conn[:3])
        i += 1
    return N, X[1:], Y[1:], BC[1:], np.array(Ttri)-1, Tlin


# =====================================================================
# MWLS  (identico ao mwls.m corrigido: \ + Tikhonov + g_eff)
# =====================================================================
def mwls_M(xi, yi, Xs, Ys, p=2, g=0.1):
    """
    Coordenadas locais NORMALIZADAS por dm = max(d).
    Sem isso cond(P'WP) ~ dm^-6 e a solucao DIVERGE sob refino de malha.
    Base:  [1, h, k, h^2/2, h*k, k^2/2]  com h = dm*h_chapeu
        => P = P_chapeu * diag(1, dm, dm, dm^2, dm^2, dm^2)
        => M = diag(1, 1/dm, 1/dm, 1/dm^2, 1/dm^2, 1/dm^2) * M_chapeu
    """
    h = Xs - xi; k = Ys - yi
    d = np.sqrt(h**2 + k**2)
    dm = d.max()
    if dm < 1e-14:
        dm = 1.0
    hn, kn, dn = h/dm, k/dm, d/dm                 # adimensionalizacao
    geff = max(g, 1e-6)
    P = np.column_stack([np.ones_like(hn), hn, kn,
                         0.5*hn**2, hn*kn, 0.5*kn**2])
    s = 6 if p == 2 else 3
    P = P[:, :s]
    w = 1.0/((dn**2 + geff**4/(dn**2 + geff**2 + 1e-15))**(p+1) + 1e-15)
    W = np.diag(w)
    A = P.T @ W @ P
    A += 1e-12*np.trace(A)*np.eye(s)              # regularizacao de Tikhonov
    Mn = np.linalg.solve(A, P.T @ W)
    esc = np.array([1.0, 1/dm, 1/dm, 1/dm**2, 1/dm**2, 1/dm**2])[:s]
    return Mn * esc[:, None]                      # volta as unidades fisicas


# =====================================================================
# SOLVER MFDM
#   bc_map: {codigo: ('D', func) | ('N', None)}
# =====================================================================
def visible_stars(X, Y, m, holes, ncand=48):
    """
    Estrela MFD com criterio de VISIBILIDADE.
    O criterio de distancia pura (star.m original) permite que a estrela
    atravesse um furo e capture nos do lado oposto, que nao tem relacao
    fisica com o no central. Em furos pequenos isso e sistematico e
    corrompe o operador. Aqui um vizinho so e aceito se o segmento que
    o liga ao no central nao cruzar o interior de nenhum furo.
    holes: lista de (xc, yc, r)
    """
    pts = np.column_stack([X, Y])
    tree = cKDTree(pts)
    N = len(X)
    _, CAND = tree.query(pts, k=min(ncand, N))
    IDX = np.zeros((N, m), dtype=int)

    for i in range(N):
        p = pts[i]
        aceitos = []
        for j in CAND[i]:
            q = pts[j]
            ok = True
            for (xc, yc, r) in holes:
                c = np.array([xc, yc])
                d = q - p
                L2 = d @ d
                if L2 < 1e-18:
                    continue
                t = np.clip((c - p) @ d / L2, 0.0, 1.0)   # ponto mais proximo
                if np.linalg.norm(p + t*d - c) < r*(1 - 1e-6):
                    ok = False                            # cruza o furo
                    break
            if ok:
                aceitos.append(j)
            if len(aceitos) == m:
                break
        while len(aceitos) < m:                           # fallback de seguranca
            for j in CAND[i]:
                if j not in aceitos:
                    aceitos.append(j)
                    break
            else:
                break
        IDX[i, :len(aceitos)] = aceitos
    return tree, IDX


def solve_mfdm(X, Y, BC, bc_map, m=16, g=0.1, normals=None, holes=None):
    N = len(X)
    pts = np.column_stack([X, Y])
    if holes:
        tree, IDX = visible_stars(X, Y, m, holes)
    else:
        tree = cKDTree(pts)
        _, IDX = tree.query(pts, k=min(m, N))
    A = lil_matrix((N, N)); b = np.zeros(N)

    for i in range(N):
        code = BC[i]
        kind = bc_map.get(code, ('I', None))[0] if code != 0 else 'I'

        if kind == 'D':
            A[i, i] = 1.0
            b[i] = bc_map[code][1](X[i], Y[i])
            continue

        idx = IDX[i]
        M = mwls_M(X[i], Y[i], X[idx], Y[idx], p=2, g=g)

        if kind == 'N':                            # n . grad(T) = 0
            nx, ny = normals[i]
            A[i, idx] += nx*M[1, :] + ny*M[2, :]
            b[i] = 0.0
        else:                                      # interno: laplaciano
            A[i, idx] += M[3, :] + M[5, :]
            b[i] = 0.0

    return spsolve(csr_matrix(A), b), tree, IDX


def fluxes_mfdm(X, Y, T, IDX, k=1.0, g=0.1):
    N = len(X); qx = np.zeros(N); qy = np.zeros(N)
    for i in range(N):
        idx = IDX[i]
        M = mwls_M(X[i], Y[i], X[idx], Y[idx], p=2, g=g)
        qx[i] = -k*(M[1, :] @ T[idx])
        qy[i] = -k*(M[2, :] @ T[idx])
    return qx, qy, np.hypot(qx, qy)


# =====================================================================
# FEM P1  (item iv)
# =====================================================================
def solve_fem(X, Y, Ttri, BC, bc_map):
    N = len(X)
    A = lil_matrix((N, N)); b = np.zeros(N)
    for tri in Ttri:
        i, j, k = tri
        xi, yi = X[i], Y[i]; xj, yj = X[j], Y[j]; xk, yk = X[k], Y[k]
        bb = np.array([yj-yk, yk-yi, yi-yj])
        cc = np.array([xk-xj, xi-xk, xj-xi])
        Ae = 0.5*abs((xj-xi)*(yk-yi) - (xk-xi)*(yj-yi))
        if Ae < 1e-14:
            continue
        Ke = (np.outer(bb, bb) + np.outer(cc, cc))/(4*Ae)
        for a in range(3):
            for c in range(3):
                A[tri[a], tri[c]] += Ke[a, c]
    for i in range(N):
        code = BC[i]
        if code != 0 and bc_map.get(code, ('I',))[0] == 'D':
            A[i, :] = 0; A[i, i] = 1.0
            b[i] = bc_map[code][1](X[i], Y[i])
    return spsolve(csr_matrix(A), b)


# =====================================================================
# NORMAIS NOS NOS DE CONTORNO (media das normais das arestas)
# =====================================================================
def boundary_normals(X, Y, Tlin, N, Ttri):
    """
    Normal externa nos nos de contorno.
    A orientacao NAO pode ser pelo centroide do dominio: em um furo a normal
    externa aponta para DENTRO do furo, e o teste do centroide a inverte.
    Criterio correto: cada aresta de contorno pertence a exatamente um
    triangulo; a normal externa aponta no sentido oposto ao 3o vertice.
    """
    edge2tri = {}
    for t in Ttri:
        for a, b_ in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edge2tri[(min(a, b_), max(a, b_))] = t

    nrm = np.zeros((N, 2))
    for ph, a, b_ in Tlin:
        a -= 1; b_ -= 1
        tx, ty = X[b_]-X[a], Y[b_]-Y[a]
        L = np.hypot(tx, ty)
        if L < 1e-14:
            continue
        nx, ny = ty/L, -tx/L
        t = edge2tri.get((min(a, b_), max(a, b_)))
        if t is None:
            continue
        op = [v for v in t if v not in (a, b_)]    # vertice oposto
        if not op:
            continue
        xm, ym = (X[a]+X[b_])/2, (Y[a]+Y[b_])/2
        if (X[op[0]]-xm)*nx + (Y[op[0]]-ym)*ny > 0:   # aponta p/ o interior
            nx, ny = -nx, -ny
        nrm[a] += (nx, ny); nrm[b_] += (nx, ny)

    L = np.hypot(nrm[:, 0], nrm[:, 1])
    L[L < 1e-14] = 1.0
    return nrm/L[:, None]
