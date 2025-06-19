# ============================================================
#   PROGRAMA 6 - Métodos Numéricos de Otimização para a Função Objetivo f(x,y)=4x+2y+x²−2x⁴+2xy−3y²
#   Métodos implementados:
#   1) Busca Aleatória (grade ou uniforme randômica)
#   2) Aclive Máximo (gradiente normalizado + linha de busca exata)
#   3) Gradientes Conjugados (Fletcher–Reeves) com
# ============================================================


# ---------------------------------------------------------------------
# 0 | Imports
# ---------------------------------------------------------------------
from __future__ import annotations
import time
from typing import Literal, NamedTuple, Callable, List
import matplotlib.pyplot as plt 

import numpy as np
import pandas as pd
import math

# ---------------------------------------------------------------------
# 1 | Função-alvo, gradiente e hessiana
# ---------------------------------------------------------------------
def f(points: np.ndarray) -> np.ndarray:
    """Função objetivo:  f(x,y)=4x+2y+x²−2x⁴+2xy−3y²"""
    x = points[..., 0]
    y = points[..., 1]
    return 4*x + 2*y + x**2 - 2*x**4 + 2*x*y - 3*y**2


def grad(p: np.ndarray) -> np.ndarray:        # usado nos métodos 2 e 3
    x, y = p
    return np.array([4 + 2*x - 8*x**3 + 2*y,
                     2 + 2*x - 6*y])


def hessian(p: np.ndarray) -> np.ndarray:     # para classificação do ponto crítico
    x, _ = p
    return np.array([[2 - 24*x**2, 2],
                     [2,           -6]])

# ---------------------------------------------------------------------
# 2 | Tipos auxiliares
# ---------------------------------------------------------------------
Distribution = Literal["regular", "random"]


class SearchResult(NamedTuple):
    method:       str
    n_points:     int
    distribution: Distribution
    point:        np.ndarray
    value:        float
    cost:         float        # wall-clock em segundos
    rel_error:    float | None

class SearchResultSA(NamedTuple):
    method: str
    n_iter: int
    point: np.ndarray
    value: float
    cost: float
    traj: list[np.ndarray]

class SearchResultCG(NamedTuple):
    method: str
    n_iter: int
    point: np.ndarray
    value: float
    cost: float
    traj: List[np.ndarray]

# ---------------------------------------------------------------------
# 3 | Método 1 – Busca Aleatória
# ---------------------------------------------------------------------
def _generate_points(n: int, dist: Distribution,
                     low: float = -3.0, high: float = 3.0,
                     rng: np.random.Generator | None = None) -> np.ndarray:
    """Cria (n,2) pontos no domínio especificado."""
    if rng is None:
        rng = np.random.default_rng()

    if dist == "regular":
        k = int(np.ceil(np.sqrt(n)))
        grid = np.linspace(low, high, k)
        xx, yy = np.meshgrid(grid, grid, indexing="xy")
        pts = np.column_stack([xx.ravel(), yy.ravel()])[:n]
    else:                                    # dist == "random"
        pts = rng.uniform(low, high, size=(n, 2))
    return pts


def random_search(n: int, dist: Distribution,
                  cost_fn: Callable[[], float] = time.perf_counter,
                  rng: np.random.Generator | None = None) -> SearchResult:
    """Executa a busca aleatória (grade ou uniforme randômica)."""
    t0 = cost_fn()
    pts  = _generate_points(n, dist, rng=rng)
    vals = np.apply_along_axis(f, 1, pts)
    idx  = np.argmax(vals)                              # buscamos o máximo
    elapsed = cost_fn() - t0
    return SearchResult("RandomSearch", n, dist,
                        pts[idx], vals[idx], elapsed, rel_error=None)

# ---------------------------------------------------------------------
# 4 | Método 2 – Aclive Máximo (linha de busca exata)
#      – gradiente normalizado
#      – Newton–Raphson 1-D (até 50 it.) com hessiana numérica de alta
#        precisão (eps = 1e-10)
#      – fallback: Golden Section (máx 200 it., tol = 1e-12)
# ---------------------------------------------------------------------

def _golden_section_max(phi, a: float, b: float,
                        tol: float = 1e-12, max_iter: int = 200) -> float:
    """Máximo de phi em [a,b] (phi unimodal) via Golden Section."""
    gr = (math.sqrt(5) - 1) / 2           # razão áurea conjugada ≈ 0.618…
    c  = b - gr * (b - a)
    d  = a + gr * (b - a)
    fc, fd = phi(c), phi(d)

    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if fc < fd:                      # mantém [c,b]
            a, c = c, d
            fc  = fd
            d   = a + gr * (b - a)
            fd  = phi(d)
        else:                            # mantém [a,d]
            b, d = d, c
            fd  = fc
            c   = b - gr * (b - a)
            fc  = phi(c)
    return 0.5 * (a + b)


def _find_h_star(xk: np.ndarray, gk: np.ndarray,
                 newton_max: int = 50) -> float:
    """
    Resolve g'(h)=0 via Newton com segunda derivada numérica;
    se falhar, aplica Golden Section no intervalo [0, h_try*4].
    """
    def dphi(h: float) -> float:
        return grad(xk + h * gk).dot(gk)

    h = 1.0e-2                                # chute inicial
    for _ in range(newton_max):
        d1 = dphi(h)
        if abs(d1) < 1.0e-14:                 # atingiu precisão
            return h
        eps = 1.0e-10
        d2 = (dphi(h + eps) - d1) / eps       # segunda derivada
        if abs(d2) < 1.0e-20:                 # Hessiana ~0 → instável
            break
        h_new = h - d1 / d2
        if h_new <= 0 or h_new > 10:          # passo inviável
            break
        h = h_new

    # --- fallback -----------------------------------------------------
    phi = lambda hh: f(xk + hh * gk)
    return _golden_section_max(phi, 0.0, h * 4.0)


def steepest_ascent(
        x0: np.ndarray,
        tol: float = 1e-8,
        max_iter: int = 1000
) -> SearchResultSA:
    """
    Método do Aclive Máximo com linha de busca exata.
    Direção = gradiente normalizado; passo ótimo = h*.
    """
    t0 = time.perf_counter()
    x   = np.asarray(x0, dtype=float)
    traj = [x.copy()]

    for k in range(1, max_iter + 1):
        gk   = grad(x)
        gnorm = np.linalg.norm(gk)
        if gnorm < tol:
            break

        gk_unit = gk / gnorm                   # normaliza direção
        h_star  = _find_h_star(x, gk_unit)
        x       = x + h_star * gk_unit
        traj.append(x.copy())

    elapsed = time.perf_counter() - t0
    return SearchResultSA("SteepestAscent", k, x, f(x), elapsed, traj)  # type: ignore[arg-type]

# -----------------------------------------------------------
# Método dos Gradientes Conjugados (Fletcher–Reeves)
# -----------------------------------------------------------
def conjugate_gradients(x0: np.ndarray,
                        tol: float = 1e-6,
                        max_iter: int = 200,
                        restart_m: int = 10) -> SearchResultCG:
    """
    Versão para maximização (d_k orientado pelo gradiente).
    Linha de busca exata via _find_h_star.
    """
    t0 = time.perf_counter()
    x  = np.asarray(x0, dtype=float)
    g  = grad(x)
    d  = g.copy()                  # direção inicial = gradiente
    traj = [x.copy()]

    for k in range(1, max_iter + 1):
        if np.linalg.norm(g) < tol:
            break
       
        h_star = _find_h_star(x, d)
        x_new  = x + h_star * d
        g_new  = grad(x_new)

        beta = (g_new @ g_new) / (g @ g)      # Fletcher–Reeves
        if beta < 0 or k % restart_m == 0:    # reinicia a direção
            beta = 0.0
        d    = g_new + beta * d               # nova direção

        x, g = x_new, g_new
        traj.append(x.copy())

    elapsed = time.perf_counter() - t0
    return SearchResultCG("ConjugateGradients", k, x, f(x), elapsed, traj)

# ---------------------------------------------------------------------
# 6 | Pipeline – executa sucessivamente os três métodos
# ---------------------------------------------------------------------
NP_VECTOR = [50, 200, 350, 500, 700, 1000]


def run_random_block() -> pd.DataFrame:
    rows: list[dict] = []
    for dist in ("regular", "random"):
        prev_val: float | None = None
        for n in NP_VECTOR:
            res = random_search(n, dist)
            if prev_val is not None:
                rel = abs(res.value - prev_val) / abs(prev_val)
                res = res._replace(rel_error=rel)
            rows.append(res._asdict())
            prev_val = res.value
    return pd.DataFrame(rows)


def banner():
    print("\n" + "="*66)
    print("||{:^62}||".format(""))
    print("||{:^62}||".format("PROGRAMA 6"))
    print("||{:^62}||".format("MÉTODOS NUMÉRICOS DE OTIMIZAÇÃO"))
    print("||{:^62}||".format("f(x,y) = 4x + 2y + x² − 2x⁴ + 2xy − 3y²"))
    print("||{:^62}||".format("1) Busca Aleatória"))
    print("||{:^62}||".format("2) Aclive Máximo (linha de busca exata)"))
    print("||{:^62}||".format("3) Gradientes Conjugados (Fletcher–Reeves)"))
    print("||{:^62}||".format("ENG. LUCAS WANICK  —  MESTRADO EM ENG. MECÂNICA • UnB"))
    print("||{:^62}||".format(""))
    print("="*66 + "\n")


def executar_metodo_1():
    df = run_random_block()
    print("\n==== MÉTODO 1 – BUSCA ALEATÓRIA ====")
    pd.set_option("display.float_format", "{:.6e}".format)
    pd.set_option("display.width", 120)
    print(df.to_string(index=False))


def executar_metodo_2():
    res = steepest_ascent(np.array([0.0, 0.0]))
    print("\n==== MÉTODO 2 – ACLIVE MÁXIMO ====")
    print(f"Método:        {res.method}")
    print(f"Iterações:     {res.n_iter}")
    print(f"Ponto ótimo:   ({res.point[0]:.6f}, {res.point[1]:.6f})")
    print(f"Valor máximo:  {res.value:.6f}")
    print(f"Custo (s):     {res.cost:.3e}")
    print(f"Tamanho trilha:{len(res.traj)}")

    # plot
    xx = np.linspace(-3, 3, 300)
    yy = np.linspace(-3, 3, 300)
    X, Y = np.meshgrid(xx, yy)
    Z = f(np.stack([X, Y], axis=-1))
    xs, ys = zip(*res.traj)
    plt.contour(X, Y, Z, levels=25)
    plt.plot(xs, ys, "-o", ms=3)
    plt.xlabel("x"); plt.ylabel("y")
    plt.title("Caminho – Aclive Máximo")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.show()


def executar_metodo_3():
    res = conjugate_gradients(np.array([0.0, 0.0]))
    print("\n==== MÉTODO 3 – GRADIENTES CONJUGADOS ====")
    print(f"Método:        {res.method}")
    print(f"Iterações:     {res.n_iter}")
    print(f"Ponto ótimo:   ({res.point[0]:.6f}, {res.point[1]:.6f})")
    print(f"Valor máximo:  {res.value:.6f}")
    print(f"Custo (s):     {res.cost:.3e}")
    print(f"Tamanho trilha:{len(res.traj)}")

    # plot
    xx = np.linspace(-3, 3, 300)
    yy = np.linspace(-3, 3, 300)
    X, Y = np.meshgrid(xx, yy)
    Z = f(np.stack([X, Y], axis=-1))
    xs, ys = zip(*res.traj)
    plt.contour(X, Y, Z, levels=25)
    plt.plot(xs, ys, "-o", ms=3)
    plt.xlabel("x"); plt.ylabel("y")
    plt.title("Caminho – Gradientes Conjugados")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.show()


def main():
    banner()
    input("Pressione <Enter> para iniciar o Método 1...")

    executar_metodo_1()
    resp = input("\nExecutar Método 2 (Aclive Máximo)? [s/N] ").strip().lower()
    if resp != 's':
        print("\nPrograma encerrado pelo usuário.\n")
        return

    executar_metodo_2()
    resp = input("\nExecutar Método 3 (Gradientes Conjugados)? [s/N] ").strip().lower()
    if resp != 's':
        print("\nPrograma encerrado pelo usuário.\n")
        return

    executar_metodo_3()
    print("\nTodos os métodos foram executados. Programa finalizado.\n")

# ---------------------------------------------------------------------
# 7 | Ponto de entrada
# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()
