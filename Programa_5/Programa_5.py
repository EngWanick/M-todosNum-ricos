# ============================================================
#   PROGRAMA 5 - SOLUÇÃO DE PROBLEMAS LINEARES - ESTUDO DE CASO
#   a) Concentração dos Reatores - Gauss-Seidel
#   b) Condução 1D Transiente com Geração Interna - Thomas + Crank-Nicolson + Validação com solução exata
#   Módulo integrado com menu interativo
# ============================================================

import numpy as np
import matplotlib.pyplot as plt

# ======== FUNÇÕES DO PROBLEMA DE CONDUÇÃO TRANSIENTE ========

# Função Thomas (inalterada)
def thomas(a, b, c, d):
    n=len(b); c_star=np.empty(n); d_star=np.empty(n); x=np.empty(n)
    c_star[0]=c[0]/b[0]; d_star[0]=d[0]/b[0]
    for i in range(1,n):
        denom=b[i]-a[i]*c_star[i-1]
        c_star[i]=c[i]/denom if i<n-1 else 0.
        d_star[i]=(d[i]-a[i]*d_star[i-1])/denom
    x[-1]=d_star[-1]
    for i in range(n-2,-1,-1): x[i]=d_star[i]-c_star[i]*x[i+1]
    return x

# Raízes e série analítica
def mu_roots(Bi, n_roots=50, tol=1e-12):
    roots=[]; guess=1.5
    for _ in range(n_roots):
        for _ in range(100):
            f  = guess*np.tan(guess)-Bi
            df = np.tan(guess)+guess/np.cos(guess)**2
            step=f/df; guess-=step
            if abs(step)<tol: break
        roots.append(guess); guess+=np.pi
    return np.array(roots)

def theta_series(X, tau, Bi, n_terms=50):
    mu=mu_roots(Bi,n_terms)
    An=4*np.sin(mu)/(2*mu+np.sin(2*mu))
    cos=np.cos(np.outer(mu,X))[:, :, None]
    exp=np.exp(-np.outer(mu**2,tau))[:, None, :]
    return (An[:,None,None]*cos*exp).sum(axis=0)

# Solver de condução
def run_simulation(N, dt, qdot, tempos,
                   Bi_eff, L=10e-3, k=30, alpha=5e-6,
                   h=1100, Tinf=250, T0=25, beta=0.5, tol=1e-6):
    dx=L/(N-1); Fo=alpha*dt/dx**2
    Agen=Fo*(qdot*dx**2/k)
    T=np.full(N,T0); t=0; snaps={}
    next_idx=0; tmax=max(tempos)
    while t<tmax+1e-12:
        Told=T.copy()
        Bi_FV = h*dx/k
        a=np.full(N,-beta*Fo); a[0]=0; a[-1]=-2*beta*Fo
        b=np.full(N,1+2*beta*Fo); b[-1]+=2*beta*Fo*Bi_FV
        c=np.full(N,-beta*Fo); c[0]=-2*beta*Fo; c[-1]=0
        d=Told+(1-beta)*Agen
        d[0]=Told[0]+Agen
        d[-1]=Told[-1]+Agen+2*beta*Fo*Bi_FV*Tinf
        T=thomas(a,b,c,d)
        if np.max(np.abs(T-Told))<tol:
            for τ in tempos[next_idx:]: snaps[τ]=T.copy(); next_idx=len(tempos)
            break
        t+=dt
        while next_idx<len(tempos) and t>=tempos[next_idx]-dt/2:
            snaps[tempos[next_idx]]=T.copy(); next_idx+=1
    return snaps, Fo, dx

# ======== PROBLEMA: CONCENTRAÇÃO DOS REATORES ===============

def run_concentracao():
    print("\n--- RESOLVENDO: CONCENTRAÇÃO DOS REATORES ---")

    # Dados fornecidos
    Q = {
        'Q01': 5, 'Q03': 8, 'Q12': 3, 'Q15': 3, 'Q23': 1, 'Q24': 1, 'Q25': 1, 
        'Q31': 1, 'Q34': 8, 'Q44': 11, 'Q54': 2, 'Q55': 2
    }
    c01, c03 = 10, 20

    # Montagem da matriz
    A = np.array([
        [Q['Q12']+Q['Q15'], 0, -Q['Q31'], 0, 0],
        [-Q['Q12'], Q['Q23']+Q['Q24']+Q['Q25'], 0, 0, 0],
        [0, -Q['Q23'], Q['Q31']+Q['Q34'], 0, 0],
        [0, -Q['Q24'], -Q['Q34'], Q['Q44'], -Q['Q54']],
        [-Q['Q15'], -Q['Q25'], 0, 0, Q['Q54']+Q['Q55']]
    ])

    b = np.array([
        Q['Q01']*c01,
        0,
        Q['Q03']*c03,
        0,
        0
    ])

    tol = 1e-12
    max_iter = 1000
    x = np.zeros(len(b))

    for iteration in range(max_iter):
        x_new = np.copy(x)
        for i in range(A.shape[0]):
            s1 = sum(A[i, j] * x_new[j] for j in range(i))
            s2 = sum(A[i, j] * x[j] for j in range(i + 1, A.shape[1]))
            x_new[i] = (b[i] - s1 - s2) / A[i, i]
        if np.all(np.abs(x_new - x) / np.abs(x_new + 1e-10) < tol):
            break
        x = x_new

    print("\n--- RESULTADOS FINAIS ---")
    for i, conc in enumerate(x, start=1):
        print(f"c{i} = {conc:.6f} mg/m³")
    print(f"\nNúmero de iterações: {iteration + 1}\n")

# ======== PROBLEMA: CONDUÇÃO TRANSIENTE =====================

def run_conducao():
    print("\n--- RESOLVENDO: CONDUÇÃO TRANSIENTE COM GERAÇÃO ---")

    L=10e-3; k=30; h=1100; alpha=5e-6
    N=41; dt=0.003
    tempos_base=[1,5,7,10,15,20,30]#,120,300,600,1000]

    dx=L/(N-1); Bi_FV=h*dx/k
    snap_base, Fo, dx = run_simulation(N, dt, qdot=1e7,
                                       tempos=tempos_base,
                                       Bi_eff=Bi_FV)

    print(f'\nN={N:<3}  Δx={dx:.4e} m   Δt={dt:.4e} s   Fo={Fo:.3f}\n')

    x_mm=np.linspace(0,10,N)
    for t,Tv in snap_base.items():
        print(f'Tempo = {t:>6.1f} s | T(0) = {Tv[0]:.2f} °C | T(L) = {Tv[-1]:.2f} °C')

    plt.figure(figsize=(8,5))
    for t,Tv in snap_base.items():
        plt.plot(x_mm,Tv,label=f'{t:g}s')
    plt.title('Condução Transiente com Geração Interna')
    plt.xlabel('x (mm)')
    plt.ylabel('T (°C)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Validação
    resp = input('\nDeseja realizar validação com solução exata? [s/N] ').strip().lower()
    if resp == 's':
        tempos_val=[5,30,120,500]
        Bi_plate=h*L/k
        snaps_val,_,_=run_simulation(N,dt,qdot=0.0,
                                     tempos=tempos_val,
                                     Bi_eff=Bi_plate)

        X=np.linspace(0,1,N)
        tau=np.array(tempos_val)*alpha/L**2
        theta=theta_series(X,tau,Bi_plate,n_terms=50)
        T_exact=theta*(25-250)+250

        plt.figure(figsize=(8,5))
        for j,t in enumerate(tempos_val):
            plt.plot(x_mm,snaps_val[t],'o-',label=f'num {t:g}s')
            plt.plot(x_mm,T_exact[:,j],'--',label=f'anal {t:g}s')
            err=np.max(np.abs(snaps_val[t]-T_exact[:,j]))
            print(f'Erro máx em {t:g}s → {err:.3f} °C')
        plt.title('Validação: numérico × analítico  (q_dot = 0)')
        plt.xlabel('x (mm)')
        plt.ylabel('T (°C)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
    else:
        print('Validação pulada.')

# ======== BANNER E MENU PRINCIPAL ===========================

def banner():
    print("\n" + "="*62)
    print("||{:^58}||".format(""))
    print("||{:^58}||".format("PROGRAMA 5"))
    print("||{:^58}||".format("SOLUÇÃO DE PROBLEMAS LINEARES"))
    print("||{:^58}||".format("ESTUDO DE CASO"))
    print("||{:^58}||".format("a) Concentração dos Reatores                 "))
    print("||{:^58}||".format("b) Condução 1D Transiente com Geração Interna"))
    print("||{:^58}||".format("ENG.WANICK"))
    print("||{:^58}||".format("MÉTODOS NUMÉRICOS"))
    print("||{:^58}||".format("MESTRADO EM ENGENHARIA MECÂNICA"))
    print("||{:^58}||".format("UnB"))
    print("||{:^58}||".format(""))
    print("="*62 + "\n")

def main():
    banner()
    opcoes = {'1': run_concentracao, '2': run_conducao}

    escolha = input("Escolha a opção para resolver:\n"
                    "1 - Concentração dos Reatores\n"
                    "2 - Condução Transiente com Geração\n"
                    "Opção: ").strip()

    if escolha not in opcoes:
        print("\nOpção inválida. Encerrando programa.\n")
        return

    opcoes[escolha]()  # roda a opção escolhida

    # Perguntar se deseja rodar a outra
    outra = '2' if escolha == '1' else '1'
    resp = input(f"\nDeseja resolver a opção ({outra})? [s/N] ").strip().lower()
    if resp == 's':
        opcoes[outra]()
    else:
        print("\nPrograma finalizado.\n")

# ======== EXECUÇÃO ==========================================

if __name__ == '__main__':
    main()
