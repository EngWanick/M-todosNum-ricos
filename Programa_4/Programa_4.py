# Programa 4:
# Feito em Python 
# Método de Bairstow para cálculo de raízes de funções polinomiais de grau η e plotagem de fractais de iterações para cada par de chutes iniciais (r,s).
# Aluno : Lucas Wanick 242104541
# Repositório do GitHub -> https://github.com/EngWanick/M-todosNum-ricos/tree/main/Programa_4
# -------------------------------------------------------------------------------------------

# ====== SECTION 1 =======
#1.1 Importando bibliotecas
import os  # Configura o ambiente para limitar o número de threads paralelos

# Evita uso excessivo de threads por bibliotecas numéricas em sistemas multicore:
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # OpenBLAS: implementação otimizada da BLAS
os.environ["OMP_NUM_THREADS"] = "1"       # OpenMP: interface para paralelismo baseado em threads
os.environ["MKL_NUM_THREADS"] = "1"       # MKL (Intel): biblioteca de álgebra linear altamente otimizada

import numpy as np
import pandas as pd                       # O Pandas é uma biblioteca Python para manipulação e análise de dados
from joblib import Parallel, delayed      # O Joblib é uma biblioteca Python para paralelização de tarefas
from tqdm import tqdm                     # O TQDM é uma biblioteca Python para exibir barras de progresso em loops
import sys                                # Para tratamento de erros e exceções

# 1.2 Configurações de exibição
def banner_inicial():
    print("\n" + "*" * 70)
    print("*{:^68}*".format(""))
    print("*{:^68}*".format("PROGRAMA 4 - MÉTODO DE BAIRSTOW"))
    print("*{:^68}*".format("LOCALIZAÇÃO DE RAÍZES REAIS E COMPLEXAS"))
    print("*{:^68}*".format("DE POLINÔMIOS DE GRAU \u03B7"))
    print("*{:^68}*".format(""))
    print("*{:^68}*".format("ENG.WANICK - MÉTODOS NUMÉRICOS AVANÇADOS"))
    print("*{:^68}*".format("MESTRADO EM CIÊNCIAS MECÂNICAS"))
    print("*{:^68}*".format("UnB"))
    print("*{:^68}*".format(""))
    print("*" * 70)
    input("Pressione ENTER para iniciar o programa...\n")

banner_inicial()

# ====== SECTION 2 =======
#2.1 Definindo funções
# 2.1.1 quadroots - Calcula as raízes de um polinômio quadrático
def quadroots(r, s):
    disc = r**2 + 4 * s
    if disc >= 0:
        root1 = (r + np.sqrt(disc)) / 2
        root2 = (r - np.sqrt(disc)) / 2
        return [root1, root2]
    else:
        real = r / 2
        imag = np.sqrt(abs(disc)) / 2
        return [complex(real, imag), complex(real, -imag)]

# 2.1.2 evaluate_polynomial - Avalia o polinômio em função da raíz encontrada para avaliar a convergência
def evaluate_polynomial(coeffs, x):
    result = 0
    for coef in coeffs:
        result = result * x + coef
    return result

# 2.1.3 bairstow
def bairstow(a, r, s, es=1e-6, maxit=1000):
    n = len(a) - 1
    roots = []
    a = a.copy()                             # Evita modificar a lista original
    iter_total = 0

    while n >= 3:                            # Loop principal   
        iter = 0
        kicks = 0
        ea1 = ea2 = 100.0
        while iter < maxit and (ea1 > es or ea2 > es):
            b = [0.0] * (n + 1)
            c = [0.0] * (n + 1)

            b[n] = a[n]
            b[n - 1] = a[n - 1] + r * b[n]

            for i in range(n - 2, -1, -1):
                b[i] = a[i] + r * b[i + 1] + s * b[i + 2]

            c[n] = b[n]
            c[n - 1] = b[n - 1] + r * c[n]

            for i in range(n - 2, -1, -1):
                c[i] = b[i] + r * c[i + 1] + s * c[i + 2]

            det = c[2] * c[2] - c[3] * c[1]
            if abs(det) < 1e-12:
                kicks += 1                    # Se o determinante for muito pequeno, incrementa-se o valor dos chutes iniciais
                if kicks > 50:
                    print(f"[INFO] Bairstow falhou em convergir para o polinômio: {print_polynomial(a)}")
                    return [], iter_total
                r += 0.05
                s += 0.05            
                continue

            dr = (-b[1] * c[2] + b[0] * c[3]) / det
            ds = (-b[0] * c[2] + b[1] * c[1]) / det

            r += dr
            s += ds
            iter += 1

            if r != 0:
                ea1 = abs(dr / r) * 100
            if s != 0:
                ea2 = abs(ds / s) * 100
        iter_total += iter

        quad = quadroots(r, s)                 # Calcula as raízes do polinômio quadrático de bairstow em termos de r e s
        
        res_quad = sum(abs(evaluate_polynomial(a, root)) for root in quad) # Avalia o polinômio nos raízes encontrada em quadroots() p/ verificar a convergência
        if res_quad < 1e-3:                                                # Se o "resíduo" for menor que uma tolerância, considera-se que a raiz foi encontrada
            roots.extend(quad)
            a = b[2:]
            n -= 2
        else:                                                              # Se não, resolve para um caso linear de um polinômio (x-r) - Horner
            root_guess = r                                                 # Adotou-se este método pois Bairstow estava retornando raízes diferentes das esperadas
            roots.append(root_guess)
            b = [0.0] * (n + 1)
            b[n] = a[n]
            for i in range(n - 1, -1, -1):
                b[i] = a[i] + root_guess * b[i + 1]
            a = b[1:]
            n -= 1

    if n == 2:                                 # Quando n chega em 2, para um polinômio de grau par, resolve diretamente com os novos a's.
        if abs(a[2]) > 1e-12:
            r = -a[1] / a[2]
            s = -a[0] / a[2]
            roots.extend(quadroots(r, s))
        elif abs(a[1]) > 1e-12:
            roots.append(-a[0] / a[1])
        else:
            pass
    elif n == 1:                               # Para polinômios ímpares, ao final, restará um polinômio linear
        roots.append(-a[0] / a[1])
    elif n > 0:
        try:
            fallback_roots = np.roots(a)
            roots.extend(fallback_roots)
        except Exception as e:
            print(f"Fallback com no.roots falhou: {e}")
    return roots, iter_total

# 2.1.4 grid_search_bairstow - Cria e percorre a malha de chutes (r, s)
def grid_search_bairstow(coeffs, r_range, s_range, step=0.0, es=1e-6, maxit=1000):
    if all(abs(c) < 1e-12 for c in coeffs):
        print(f"[ERRO] Muitos coeficientes são ~0. Ajuste o polinômio e tente novamente.")
        sys.exit(1)  # Todos coeficientes são ~0, não há mais raízes a extrair
    
    r_values = np.arange(r_range[0], r_range[1] + step, step)
    s_values = np.arange(s_range[0], s_range[1] + step, step)
    pontos = [(r, s) for r in r_values for s in s_values]

    def avaliar_ponto(r, s):
        roots, iter_count = bairstow(coeffs, r, s, es, maxit)
        residual = sum(abs(evaluate_polynomial(coeffs, root)) for root in roots) \
        if roots else np.inf
        return (r, s, residual, iter_count, roots)

        
    print(f"[INFO] Polinômio em análise: {print_polynomial(coeffs)} de grau {len(coeffs)-1}")      
    print(f"[INFO] Processando {len(pontos)} pontos...")
    resultados = Parallel(n_jobs=-1)(delayed(avaliar_ponto)(r, s) for r, s in tqdm(pontos))

    resultados = [res for res in resultados if res is not None]

    if not resultados:
        try:
            fallback_roots = np.roots(coeffs)
            if len(fallback_roots) > 0:
                print("\n AVISO: fallback activado com np.roots - método de Bairstow falhou em convergir para este polinômio.")
                return 0, 0, fallback_roots, 0.0, [(r, s, float('inf'), -1) for r, s in pontos]
        except:
            return None, None, [], float('inf'), []
        
    best = min(resultados, key=lambda x: x[2])
    best_r, best_s, min_residual, _, best_roots = best
    residual_map = [(r, s, residual, iters) for r, s, residual, iters, _ in resultados]
    
    return best_r, best_s, best_roots, min_residual, residual_map

# 2.1.5 print_polynomial - apenas para printar o polinômio no terminal
def print_polynomial(coeffs):
    n = len(coeffs) - 1
    terms = []

    for i, c in enumerate(coeffs):
        exp = n - i
        if c == 0:
            continue
        term = (
            f"{c}" if exp == 0 else
            (f"x" if c == 1 and exp == 1 else
             f"-x" if c == -1 and exp == 1 else
             f"x^{exp}" if c == 1 else
             f"-x^{exp}" if c == -1 else
             f"{c}x" if exp == 1 else
             f"{c}x^{exp}")
        )
        terms.append(term)

    return " + ".join(terms).replace("+ -", "- ").replace(" 1x", " x").replace(" -1x", " -x").replace(" + -", " - ").replace(" 1x^", " x^").replace(" -1x^", " -x^")

# 2.1.6 exportar_csv_iteracoes
def exportar_csv_iteracoes(residual_map, filename="bairstow_iteracoes.csv", ativar=True):
    if not ativar:
        return
    df = pd.DataFrame(residual_map, columns=["r", "s", "residual", "iterations"])
    df.to_csv(filename, index=False)
    print(f"[INFO] CSV de iterações exportado: {filename}")

# 2.1.7 export_txt
def export_txt(residual_map, filename="bairstow_iteracoes.txt", activate=True):
    """
    Gera um arquivo de texto com três colunas:
      r   s   iterations
    para uso direto no gnuplot:
      splot "bairstow_iteracoes.txt" using 1:2:3 with dots
    """
    if not activate:
        return
    residual_map.sort(key=lambda x: (x[0], x[1]))  # Ordena por r e s
    with open(filename, "w") as f:
        # cabeçalho (opcional)
        f.write("# r s iterations\n")
        last_r = None
        for r, s, _, iters in residual_map:
            if last_r is None or r != last_r:
                f.write(f"{r:.6f} {s:.6f} {int(iters)}\n")
                last_r = r
    print(f"[INFO] TXT de iterações exportado como: {filename}")

# 2.1.8 Fechar o programa
def safe_input(prompt=""):
    resposta = input(prompt).strip()
    if resposta.lower() == "exit":
        print("Encerrando o programa conforme solicitado.")
        sys.exit(0)
    return resposta

# ====== SECTION 3 =======
#3.1 Plotagem do Fractal

# 3.1.1 plot_iteracoes_fractal
def plot_iteracoes_fractal(residual_map, r_range, s_range, resolution=(1000, 1000)):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from datetime import datetime
    
    print("\n[INFO] Gerando gráfico de alta resolução (1000x1000) com o número de iterações...")

    # Inicializa a matriz de iterações
    iter_map = np.full((resolution[1], resolution[0]), np.nan)

    # Preenche a matriz com os dados de iteração
    for r, s, _, iters in tqdm(residual_map, desc="Montando matriz de iterações", unit="pts"):
        i = int((s - s_range[0]) / (s_range[1] - s_range[0]) * (resolution[1] - 1))
        j = int((r - r_range[0]) / (r_range[1] - r_range[0]) * (resolution[0] - 1))
        if 0 <= i < resolution[1] and 0 <= j < resolution[0]:
            iter_map[i, j] = iters

    # Plotagem
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(
        iter_map,
        extent=(r_range[0], r_range[1], s_range[0], s_range[1]),
        origin='lower',
        cmap='plasma',          # Paleta profissional
        interpolation='bilinear',         # Suavização
        aspect='equal',
        alpha=0.9                          # Transparência leve
    )

    # Título e rótulos
    ax.set_title(f"Fractal de Bairstow - Polinômio {print_polynomial(coeffs_test)}", fontsize=14)
    ax.axhline(0, color='black', linewidth=1, alpha=0.6)  # eixo s
    ax.axvline(0, color='black', linewidth=1, alpha=0.6)  # eixo r
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))



    # Barra de cor
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Número de Iterações", fontsize=12)

    plt.tight_layout()

    # Nome do arquivo com timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"bairstow_fractal_iteracoes_{timestamp}.png"

    # Salvamento com fundo branco e ajustes visuais
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[INFO] Fractal salvo como: {filename}")

    plt.show()

# ====== SECTION 4 =======

#4.1 Testes de polinômios
#[6, 11, -33, -33, 11, 6] [1, 2, 4, 8] [1, 0, 0, 0, 0, 1, 0] [1, -25, 230, -950, 1689, -945] [1, -10, 35, -50, 24] [1, 0, -1]
#[1, 2, -5, -6] [1,-5 ,6] [1, 1, -6] [1, 2, -7, -8, 12] [1, -5, 5, 5, -6, 0]

#4.2 Execução
if __name__ == "__main__":
    # Solicita e valida coeficientes do polinômio com validação
    print("Entre com os coeficientes do polinômio separados por vírgula, do maior para o menor grau(ex: 1,2,3):")
    while True:
        coeffs_input = safe_input()
        try:
            coeffs_test = [float(c.strip()) for c in coeffs_input.split(",")]
            if len(coeffs_test) < 2:
                raise ValueError("O polinômio deve ter pelo menos dois coeficientes.")
            break
        except Exception as e:
            print(f"Entrada inválida para os coeficientes. Formato esperado: a,b,c,...")
            # Solicita e valida intervalo de busca
    print("Entre com o intervalo de busca dos chutes iniciais de r e s (mín, máx): ")
    while True:
        interval_input = safe_input()
        try:
            lo, hi = [float(x) for x in interval_input.split(",")]
            if lo >= hi:
                raise ValueError("O valor mínimo deve ser menor que o máximo.")
            break
        except Exception as e:
            print(f"Entrada inválida para intervalo. Formato esperado: min, max (ex: -3, 3)")
    r_range = (lo, hi)
    s_range = (lo, hi)

    # Resolução fixa da malha (1000x1000)
    step = (hi - lo) / 1000.0

    # execução da malha de busca de raízes e coleta de iterações
    best_r, best_s, best_roots, min_residual, residual_map = grid_search_bairstow(
        coeffs_test, r_range, s_range, step=step
    )

    # Print de Resultados
    print("== Resultados da Busca de Raízes ==")
    if best_roots is None or len(best_roots) == 0:
        print("Nenhuma raiz válida foi encontrada para esse polinômio.")
    else:
        df = pd.DataFrame({
            "Raiz": [float(r.real) if isinstance(r, complex) and r.imag == 0 else r for r in best_roots]
        })
        print(df)
        print(f"Melhor r: {best_r:.4f}, Melhor s: {best_s:.4f}, Residual mínimo: {min_residual:.4f}")

    # Plotagem do Fractal
    gerar = safe_input("Deseja gerar o gráfico de iterações? [Y / N]: ").strip().lower()
    if gerar == "y" or gerar == "Y":
        plot_iteracoes_fractal(residual_map, r_range, s_range)
    else:
        print("Plotagem do fractal cancelada.")

    # Exportação de dados
    exportar = safe_input("Deseja exportar os dados de iterações? [Y / N]: ").strip().lower()
    if exportar == "y" or exportar == "Y":
        print("Escolha o formato de exportação:")
        print("1. CSV")
        print("2. TXT")
        escolha = safe_input().strip()
        if escolha == "1":
            exportar_csv_iteracoes(residual_map, ativar=True)
        elif escolha == "2":
            export_txt(residual_map, filename="bairstow_iteracoes.txt", activate=True)
        else:
            print("Escolha inválida. Exportação cancelada.")
            print("FIM DO PROGRAMA.")
    else:
        print("Exportação cancelada.")
        print("FIM DO PROGRAMA.")


