"""
rodar_tudo.py -- reproduz os quatro itens da secao da aleta trapezoidal.

    python3 rodar_tudo.py

Requisitos: numpy, scipy   (matplotlib apenas para --figuras)
As malhas devem estar em ../malhas/
"""
import numpy as np, time, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aleta_solver import (read_msh, solve_mfdm, fluxes_mfdm,
                          solve_fem, boundary_normals)

M = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'malhas')
msh = lambda n: os.path.join(M, n)

R_PADRAO = 0.45
HOLES = lambda r: [(-1.5, 1.75, r), (1.5, 1.75, r)]


# =====================================================================
def item_i():
    print('\n' + '='*70)
    print('ITEM (i) -- VERIFICACAO NA ALETA COM FUROS')
    print('='*70)
    N, X, Y, BC, Ttri, Tlin = read_msh(msh('aleta.msh'))
    xmin, xmax = X.min(), X.max()
    Tex = lambda x, y: 500*(x - xmin)/(xmax - xmin)
    bc = {c: ('D', Tex) for c in (1, 2, 3, 4)}

    t0 = time.time()
    T, _, IDX = solve_mfdm(X, Y, BC, bc, m=16, g=0.1, holes=HOLES(R_PADRAO))
    dt = time.time() - t0
    qx, qy, qm = fluxes_mfdm(X, Y, T, IDX)

    e = np.linalg.norm(T - Tex(X, Y))/np.linalg.norm(Tex(X, Y))
    qan = -500/(xmax - xmin)
    print(f'N = {N}   nt = {len(Ttri)}   contorno = {int((BC!=0).sum())}')
    print(f'Tmin = {T.min():.4f}   Tmax = {T.max():.4f}   (exato 0 / 500)')
    print(f'erro L2 relativo = {e:.3e}')
    print(f'qx medio = {qx.mean():.5f}  | analitico = {qan:.5f} W/m2')
    print(f'erro relativo qx = {abs(qx.mean()-qan)/abs(qan)*100:.4f}%')
    print(f'tempo = {dt:.2f} s')


# =====================================================================
def itens_ii_iv():
    print('\n' + '='*70)
    print('ITEM (ii) CONVERGENCIA  +  ITEM (iv) MFDM vs FEM')
    print('MMS harmonica nao-polinomial:  T = e^x cos(y)')
    print('='*70)
    Tex = lambda x, y: np.exp(x)*np.cos(y)
    bc = {c: ('D', Tex) for c in (1, 2, 3, 4)}
    holes = HOLES(R_PADRAO)

    print(f"{'N':>6} {'h':>7} | {'MFDM e_L2':>10} {'ord':>5} {'t[s]':>6}"
          f" | {'FEM e_L2':>10} {'ord':>5} {'t[s]':>6}")
    print('-'*74)
    res = []; pm = pf = None
    for tag in ['M2.20', 'M1.70', 'M1.30', 'M1.00', 'M0.75', 'M0.55']:
        N, X, Y, BC, Ttri, Tlin = read_msh(msh(f'aleta_{tag}.msh'))
        P = np.column_stack([X, Y])
        a, b, c = P[Ttri[:, 0]], P[Ttri[:, 1]], P[Ttri[:, 2]]
        ar = 0.5*np.abs((b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1])
                        - (c[:, 0]-a[:, 0])*(b[:, 1]-a[:, 1]))
        h = np.sqrt(2*ar.mean())
        Te = Tex(X, Y); nTe = np.linalg.norm(Te)

        t0 = time.time()
        T, _, _ = solve_mfdm(X, Y, BC, bc, m=16, g=0.1, holes=holes)
        tm = time.time()-t0
        em = np.linalg.norm(T-Te)/nTe

        t0 = time.time(); Tf = solve_fem(X, Y, Ttri, BC, bc); tf = time.time()-t0
        ef = np.linalg.norm(Tf-Te)/nTe

        om = np.log(pm[1]/em)/np.log(pm[0]/h) if pm else np.nan
        of = np.log(pf[1]/ef)/np.log(pf[0]/h) if pf else np.nan
        print(f"{N:6d} {h:7.4f} | {em:10.3e} {om:5.2f} {tm:6.2f}"
              f" | {ef:10.3e} {of:5.2f} {tf:6.2f}")
        res.append((N, len(Ttri), h, em, tm, ef, tf)); pm = (h, em); pf = (h, ef)

    r = np.array(res)
    np.save('conv.npy', r)
    print(f"\nTaxa media MFDM : p = {np.polyfit(np.log(r[:,2]), np.log(r[:,3]),1)[0]:.3f}")
    print(f"Taxa media FEM  : p = {np.polyfit(np.log(r[:,2]), np.log(r[:,5]),1)[0]:.3f}")
    print(f"Razao de erro FEM/MFDM (malha mais fina): {r[-1,5]/r[-1,3]:.1f}x")
    print(f"Razao de tempo MFDM/FEM (malha mais fina): {r[-1,4]/r[-1,6]:.1f}x")


# =====================================================================
def item_iii():
    print('\n' + '='*70)
    print('ITEM (iii) -- INFLUENCIA DO RAIO NA CONCENTRACAO DE FLUXO')
    print('='*70)
    print(f"{'r':>6} {'N':>6} {'|q|max':>9} {'K':>8}    (teoria r->0: K = 2)")
    print('-'*48)
    out = []
    for r in [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75]:
        N, X, Y, BC, Ttri, Tlin = read_msh(msh(f'fino_r{r:.2f}.msh'))
        nrm = boundary_normals(X, Y, Tlin, N, Ttri)
        xmin, xmax = X.min(), X.max()
        f = lambda x, y: 500*(x - xmin)/(xmax - xmin)
        bc = {1: ('D', f), 2: ('D', f), 3: ('D', f), 4: ('N', None)}
        T, _, IDX = solve_mfdm(X, Y, BC, bc, m=16, g=0.1,
                               normals=nrm, holes=HOLES(r))
        qx, qy, qm = fluxes_mfdm(X, Y, T, IDX)
        furo = np.where(BC == 4)[0]
        K = qm[furo].max()/(500/(xmax-xmin))
        print(f"{r:6.2f} {N:6d} {qm[furo].max():9.3f} {K:8.4f}")
        out.append((r, N, qm[furo].max(), K))
    np.save('raios.npy', np.array(out))


# =====================================================================
def convergencia_de_K():
    """Mostra por que grandezas de PICO exigem malha mais fina que medias."""
    print('\n' + '='*70)
    print('CONVERGENCIA DO PROPRIO K  (r = 0,45 fixo)')
    print('='*70)
    print(f"{'malha':>8} {'N':>6} {'K':>8}")
    print('-'*26)
    for tag in ['M2.20', 'M1.70', 'M1.30', 'M1.00', 'M0.75', 'M0.55']:
        N, X, Y, BC, Ttri, Tlin = read_msh(msh(f'aleta_{tag}.msh'))
        nrm = boundary_normals(X, Y, Tlin, N, Ttri)
        xmin, xmax = X.min(), X.max()
        f = lambda x, y: 500*(x-xmin)/(xmax-xmin)
        bc = {1: ('D', f), 2: ('D', f), 3: ('D', f), 4: ('N', None)}
        T, _, IDX = solve_mfdm(X, Y, BC, bc, m=16, g=0.1,
                               normals=nrm, holes=HOLES(R_PADRAO))
        qx, qy, qm = fluxes_mfdm(X, Y, T, IDX)
        furo = np.where(BC == 4)[0]
        print(f"{tag:>8} {N:6d} {qm[furo].max()/(500/(xmax-xmin)):8.4f}")
    print('K converge por baixo: a malha grossa SUBESTIMA o pico em ~10%.')


if __name__ == '__main__':
    item_i()
    itens_ii_iv()
    item_iii()
    convergencia_de_K()
    print('\nConcluido.')
