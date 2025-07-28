import numpy as np
import matplotlib.pyplot as plt

def gerar_malha_tipo_corrigida(n, L=1.50, H=2.50, R=0.25):
    dx = L / (n - 1)
    dy = H / (n - 1)
    malha_tipo = np.full((n, n), 'I', dtype=str)

    for i in range(n):       
        for j in range(n):   
            x = j * dx
            y = H - i * dy  

            # Superior esquerdo
            if x < R and y > (H - R):
                dist_sup_esq = (x - R)**2 + (y - (H - R))**2
                if dist_sup_esq > R**2:
                    malha_tipo[i, j] = 'A'
                    continue

            # Inferior direito
            if x > (L - R) and y < R:
                dist_inf_dir = (x - (L - R))**2 + (y - R)**2
                if dist_inf_dir > R**2:
                    malha_tipo[i, j] = 'A'
                    continue

            # Bordas Dirichlet
            if i == 0 and malha_tipo[i, j] != 'A':  
                malha_tipo[i, j] = 'D'
                continue
            if i == n - 1 and malha_tipo[i, j] != 'A':  
                malha_tipo[i, j] = 'D'
                continue

            # Bordas Neumann
            if j == 0 and malha_tipo[i, j] != 'A':
                malha_tipo[i, j] = 'N'
                continue
            if j == n - 1 and malha_tipo[i, j] != 'A':
                malha_tipo[i, j] = 'N'
                continue

    return malha_tipo

# Gerar e plotar
n = 81
malha_tipo = gerar_malha_tipo_corrigida(n)

color_map = {'I': 0, 'N': 1, 'D': 2, 'A': 3}
numeric_malha = np.vectorize(color_map.get)(malha_tipo)

plt.figure(figsize=(8, 6))
cax = plt.imshow(numeric_malha, cmap='tab10', origin='upper')
cbar = plt.colorbar(cax, ticks=[0, 1, 2, 3])
cbar.ax.set_yticklabels(['I', 'N', 'D', 'A'])
plt.title("Mapa da Malha Tipo (19x19)")
plt.xlabel("j (coluna)")
plt.ylabel("i (linha)")
plt.grid(False)
plt.show()
