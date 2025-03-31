
import csv
import os
import matplotlib.pyplot as plt
import pandas as pd

# Criar pasta de saída
os.makedirs("img", exist_ok=True)

# Carregar dados
df = pd.read_csv("sedimentacao_validacao.csv")

# Função auxiliar
def filtrar(df, St=None, Re=None, h=None):
    temp = df.copy()
    if St is not None:
        temp = temp[temp["St"] == St]
    if Re is not None:
        temp = temp[temp["Re_s"] == Re]
    if h is not None:
        temp = temp[temp["h"] == h]
    return temp

# Gráfico 1: y(t) analítica vs RK4 Stokes (Re_s = 0.01) para diferentes St
for St in [0.1, 0.5, 1, 2, 5]:
    dados = filtrar(df, St=St, Re=0.01, h=0.01)
    plt.figure()
    plt.plot(dados["t"], dados["y_analitico"], label="Analítica")
    plt.plot(dados["t"], dados["y_rk4_stokes"], '--', label="RK4 Stokes")
    plt.title(f"Gráfico 1 - St = {St} | Re_s = 0.01")
    plt.xlabel("t")
    plt.ylabel("y(t)")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"img/grafico1_st_{St}.png")

# Gráfico 2: erro_stokes(t) para diferentes h
plt.figure()
for h in [0.1, 0.01, 0.001]:
    dados = filtrar(df, St=1, Re=0.01, h=h)
    plt.plot(dados["t"], dados["erro_stokes"], label=f"h = {h}")
plt.title("Gráfico 2 - Erro RK4 Stokes vs t | St = 1, Re_s = 0.01")
plt.xlabel("t")
plt.ylabel("Erro absoluto")
plt.legend()
plt.grid(True)
plt.savefig("img/grafico2_erro_stokes_por_h.png")

# Gráfico 3: y(t) quadrático com diferentes Re_s (St = 1)
plt.figure()
for Re in [0.01, 0.1, 1, 10, 100]:
    dados = filtrar(df, St=1, Re=Re, h=0.01)
    plt.plot(dados["t"], dados["y_rk4_quadratico"], label=f"Re_s = {Re}")
dados_ref = filtrar(df, St=1, Re=0.01, h=0.01)
plt.plot(dados_ref["t"], dados_ref["y_analitico"], 'k--', label="Analítica (Re→0)")
plt.title("Gráfico 3 - y(t) Quadrático | St = 1")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.legend()
plt.grid(True)
plt.savefig("img/grafico3_quadratico_res.png")

# Gráfico 4: erro_quadratico(t) com Re_s variando (St = 1)
plt.figure()
for Re in [0.01, 0.1, 1, 10, 100]:
    dados = filtrar(df, St=1, Re=Re, h=0.01)
    plt.plot(dados["t"], dados["erro_quadratico"], label=f"Re_s = {Re}")
plt.title("Gráfico 4 - Erro Quadrático vs t | St = 1")
plt.xlabel("t")
plt.ylabel("Erro absoluto")
plt.legend()
plt.grid(True)
plt.savefig("img/grafico4_erro_quadratico_res.png")

# Gráfico 5: y(t) quadrático com diferentes h | St = 1, Re_s = 10
plt.figure()
for h in [0.1, 0.01, 0.001]:
    dados = filtrar(df, St=1, Re=10, h=h)
    plt.plot(dados["t"], dados["y_rk4_quadratico"], label=f"h = {h}")
plt.title("Gráfico 5 - y(t) Quadrático com h variando | St = 1, Re_s = 10")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.legend()
plt.grid(True)
plt.savefig("img/grafico5_h_quadratico.png")

# Gráfico 6: Erro final vs h (Stokes)
dados6 = df[(df["St"] == 1) & (df["Re_s"] == 0.01)]
erro_final_vs_h = dados6.groupby("h").apply(lambda g: g[g["t"] == g["t"].max()]["erro_stokes"].values[0])
plt.figure()
plt.plot(erro_final_vs_h.index, erro_final_vs_h.values, marker='o')
plt.title("Gráfico 6 - Erro Final vs h | RK4 Stokes, St=1, Re_s=0.01")
plt.xlabel("h")
plt.ylabel("Erro absoluto final")
plt.grid(True)
plt.savefig("img/grafico6_erro_final_stokes_vs_h.png")

# Gráfico 7: Erro final vs Re_s (Quadrático)
dados7 = df[(df["St"] == 1) & (df["h"] == 0.01)]
erro_final_vs_Re = dados7.groupby("Re_s").apply(lambda g: g[g["t"] == g["t"].max()]["erro_quadratico"].values[0])
plt.figure()
plt.plot(erro_final_vs_Re.index, erro_final_vs_Re.values, marker='o')
plt.title("Gráfico 7 - Erro Final vs Re_s | RK4 Quadrático, St=1, h=0.01")
plt.xlabel("Re_s")
plt.ylabel("Erro absoluto final")
plt.grid(True)
plt.savefig("img/grafico7_erro_final_quadratico_vs_res.png")
