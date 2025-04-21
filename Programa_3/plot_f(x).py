import numpy as np
import matplotlib.pyplot as plt

# Função polinomial corrigida
def f(x):
    return x**5 - 25*x**4 + 230*x**3 - 950*x**2 + 1689*x - 945

# Geração de pontos
x_vals = np.linspace(0, 10, 400)
y_vals = f(x_vals)

# Plotagem com ticks personalizados
plt.figure(figsize=(10, 6))
plt.plot(x_vals, y_vals, label=r'$P(x) = x^5 - 25x^4 + 230x^3 - 950x^2 + 1689x - 945$', linewidth=2)
plt.axhline(0, color='gray', linestyle='--')
plt.axvline(0, color='gray', linestyle='--')
plt.xticks(np.arange(0, 11, 1))  # Mostra ticks de 1 em 1 até 10
plt.title('Gráfico do Polinômio do Quinto Grau')
plt.xlabel('x')
plt.ylabel('P(x)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
