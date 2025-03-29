import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("stokes_curves_final.csv")

for col in ["y_analytical", "y_euler", "error_abs"]:
    plt.figure(figsize=(10, 6))
    for St in sorted(df["St"].unique()):
        df_st = df[df["St"] == St]
        plt.plot(df_st["t"], df_st[col], label=f"St = {St}")
    plt.title(f"Gráfico: {col}")
    plt.xlabel("t (tempo)")
    plt.ylabel(col)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()
