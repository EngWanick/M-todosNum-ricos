import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define a função f(x)
def f(x):
    return x**5 - 25*x**4 + 230*x**3 - 950*x**2 + 1689*x - 945

# Função para calcular os coeficientes da parábola em cada iteração
def muller_iterations(f, x0, x1, x2, tol=1e-5, max_iter=10):
    iterations = []
    for _ in range(max_iter):
        h0 = x1 - x0
        h1 = x2 - x1
        delta0 = (f(x1) - f(x0)) / h0
        delta1 = (f(x2) - f(x1)) / h1
        a = (delta1 - delta0) / (h1 + h0)
        b = a * h1 + delta1
        c = f(x2)

        delta = b**2 - 4*a*c
        if delta < 0:
            break

        sqrt_delta = np.sqrt(delta)
        denom = b + sqrt_delta
        dx = -2 * c / denom
        x3 = x2 + dx
        iterations.append((x0, x1, x2, a, b, c))

        if abs(x3 - x2) < tol:
            break

        x0, x1, x2 = x1, x2, x3
    return iterations

# Prepara dados
iterations = muller_iterations(f, 0.5, 1.5, 2.0)
x_vals = np.linspace(0, 10, 500)
y_vals = f(x_vals)

# Setup da figura
fig, ax = plt.subplots()
ax.plot(x_vals, y_vals, 'k-', label='f(x)')
parabola_line, = ax.plot([], [], 'gray', label='Parábola (Müller)')
vlines = [ax.axvline(0, color='blue', linestyle='--') for _ in range(3)]
ax.axhline(0, color='black', linewidth=0.5)
ax.set_xlim(0, 10)
ax.set_ylim(-1100, 1100)
ax.legend()

# Função de animação
def update(frame):
    x0, x1, x2, a, b, c = iterations[frame]
    x_parab = np.linspace(min(x0, x2)-1, max(x0, x2)+1, 200)
    y_parab = a*(x_parab - x2)**2 + b*(x_parab - x2) + c
    parabola_line.set_data(x_parab, y_parab)

    vlines[0].set_xdata([x0, x0])
    vlines[1].set_xdata([x1, x1])
    vlines[2].set_xdata([x2, x2])
    return parabola_line, *vlines

ani = FuncAnimation(fig, update, frames=len(iterations), interval=1000, blit=True)
plt.close()
ani_path = "D:/Engenharia Mecânica/Mestrado/2º Sem/MetNum/Tarefa_4/Programa_3/img/animacao_muller.mp4"
ani.save(ani_path)

ani_path
