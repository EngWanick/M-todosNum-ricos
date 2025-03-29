# ============================
# SIMULAÇÃO DA EDO DE SEDIMENTAÇÃO
# Simulação numérica do problema de sedimentação fluido partícula 
# em sua solução analítica via série de Taylor truncada e
# fator integrante e método de Euler para comparação do resultado.
# Plotagem da curva de erro absoluto para ajuste de parâmetros
# ============================

# --- 1. Parâmetros de simulação ---
Stokes = [0.1, 0.5, 1, 2, 5]  # Valores de Stokes
h = 0.01                      # Passo de integração
t0 = 0
tf = 5
n_steps = int((tf - t0) / h) + 1

# --- 2. Função da exponencial via série de Taylor ---
def exp_taylor(x, n=100, tol=1e-15): # Foram feitos testes com 20 e 60 termos, mas a série estava divergindo com x = t/St pelas limitações do float64. A solução foi limitar em 10^15.
     sum = 0.0
     term = 1.0
     for i in range(n):
         sum += term
         term *= x / (i + 1)
         if abs(term) < tol:
            break
     return sum

def exp_negative(x): # Para establizar numericamente, foi necessário utilizar a identidade e^(-x) = 1 / e^x.
   return 1 / exp_taylor(x)


# Função da solução analítica da EDO de sedimentação - y(t) = 1 - e^(-t/St)
def analytical(t, St):
    return 1 - exp_negative(t / St)  # Para o resultado, utilizamos a definição imposta.

# --- 3. Função para aplicar o método de Euler ---
def euler(y, h, St):
    return y + h * (-y + 1) / St # Implementa a iteração y_{n+1} = y_n + h * f(t, y)
                                 # com f(t, y) = (-y + 1)/St derivado da equação de movimento

# --- 4. Função principal da simulação ---
def sim(St, h, t0, tf):
    t = t0
    y = 0
    results = []
    for _ in range(n_steps):
        y_analytical = analytical(t, St)
        error = abs(y - y_analytical)
        results.append([St, round(t, 4), round(y, 10), round(y_analytical, 10), round(error, 10)])
        y = euler(y, h, St)
        t += h
    return results

# --- 5. Executa simulação para todos os valores de St ---
all_results = []
for St in Stokes:
    all_results.extend(sim(St, h, t0, tf))

# --- 6. Exporta para CSV ---
import csv
filename = "stokes_curves_final.csv"

with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["St", "t", "y_euler", "y_analytical", "error_abs"])
    writer.writerows(all_results)

print(f"Arquivo '{filename}' gerado com sucesso.")
