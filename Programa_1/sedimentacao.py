# Simulação com o método de Runge-Kutta de 4ª ordem (RK4)
# Aplicado à equação de movimento de uma partícula sedimentando em meio viscoso
# Considerando os casos linear (Stokes) e com força quadrática (inercial)
# Implementa controle de overflow para garantir estabilidade numérica

# =============================
# Cálculo da exponencial via série de Taylor
# =============================
# Implementação da série de Taylor para e^x
def exp_taylor(x, n=100, tol=1e-15): # Foram feitos testes com 20 e 60 termos, mas a série estava divergindo com x = t/St pelas limitações do float64. A solução foi limitar em 10^15.
    sum = 1.0
    term = 1.0
    for i in range(1, n):
        term *= x / i
        sum += term
        if abs(term) < tol:
            break
    return sum
def exp_negative(x):
    # Para establizar numericamente, foi necessário utilizar a identidade e^(-x) = 1 / e^x.
    return 1 / exp_taylor(x)

# =============================
# Solução analítica para o regime de Stokes (Re → 0)
# =============================
# Solução exata da EDO linear: y(t) = 1 - exp(-t / St)
def analytical(t, St):
    return 1 - exp_negative(t / St) # Para o resultado, utilizamos a definição imposta.

# =============================
# Solução numérica com RK4 para o modelo linear (Stokes)
# =============================
def rk4_stokes(y, t, h, St):
    def f(y):
        return (-y + 1) / St
    k1 = h * f(y)
    k2 = h * f(y + 0.5 * k1)
    k3 = h * f(y + 0.5 * k2)
    k4 = h * f(y + k3)
    return y + (k1 + 2*k2 + 2*k3 + k4) / 6.0

# =============================
# Solução numérica com RK4 para o modelo com força quadrática
# =============================
def rk4_quadratic(y, t, h, St, Re_s):
    def f(y):
        # Proteção contra explosão numérica
        if abs(y) > 1e6:
            return 0.0
        return (-y - (Re_s / 2.0) * y**2 + 1.0) / St

    try:
        k1 = h * f(y)
        k2 = h * f(y + 0.5 * k1)
        k3 = h * f(y + 0.5 * k2)
        k4 = h * f(y + k3)
        result = y + (k1 + 2*k2 + 2*k3 + k4) / 6.0

        # Verificação adicional contra overflow
        if abs(result) > 1e6:
            return float('nan')
        return result
    except OverflowError:
        return float('nan')

# =============================
# Parametrização da simulação
# =============================
St_values = [0.1, 0.5, 1, 2, 5]                  # Números de Stokes representando regimes de arrasto distintos
Re_s_values = [0.01, 0.1, 1, 10, 100]            # Números de Reynolds específicos adimensionalizados pela vel. de Stokes
h_values = [0.1, 0.01, 0.001]                    # Passos de tempo para análise de estabilidade e precisão
t_max = 5.0                                      # Tempo final da simulação

# =============================
# Geração do arquivo CSV com os resultados
# =============================
archive = "sedimentacao_validacao.csv"
with open(archive, "w") as f:
    f.write("t,h,St,Re_s,y_analitico,y_rk4_stokes,y_rk4_quadratico,erro_stokes,erro_quadratico\n")

    for h in h_values:
        for St in St_values:
            for Re_s in Re_s_values:
                n_steps = int(t_max / h)
                y_a = 0.0
                y_rk4 = 0.0
                y_quad = 0.0
                t = 0.0

                for _ in range(n_steps + 1):
                    # Solução analítica (Stokes)
                    y_a = analytical(t, St)

                    # Solução numérica (RK4 linear)
                    y_rk4 = rk4_stokes(y_rk4, t, h, St)

                    # Solução numérica (RK4 quadrático)
                    y_quad = rk4_quadratic(y_quad, t, h, St, Re_s)

                    # Erros absolutos em relação à solução exata
                    erro_s = abs(y_rk4 - y_a) if not str(y_rk4) == 'nan' else float('nan')
                    erro_q = abs(y_quad - y_a) if not str(y_quad) == 'nan' else float('nan')

                    # Escrita dos dados no CSV
                    f.write(f"{t:.4f},{h},{St},{Re_s},{y_a:.8f},{y_rk4:.8f},{y_quad:.8f},{erro_s:.8f},{erro_q:.8f}\n")
                    t += h
