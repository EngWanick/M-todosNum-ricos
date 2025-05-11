# Programa 2: Métodos de Bissecção e Falsa Posição para cálculo de raízes


def banner_inicial():
    print("\n" + "*" * 60)
    print("*{:^58}*".format(""))
    print("*{:^58}*".format("MÉTODOS DA BISSECÇÃO E FALSA POSIÇÃO"))
    print("*{:^58}*".format("CÁLCULO DE RAÍZES DE POLINÔMIOS"))
    print("*{:^58}*".format(""))
    print("*{:^58}*".format("f1(x) = -0.5 * x**2 + 2.5 * x + 4.5"))
    print("*{:^58}*".format("f2(x) = 5x³ - 5x² + 6x - 2"))
    print("*{:^58}*".format("f3(x) = -25 + 82x - 90x² + 44x³ - 8x⁴ + 0.7x"))
    print("*{:^58}*".format("f4(x) = sin(x) - x³"))
    print("*{:^58}*".format("f5(x) = ln(x⁴) - 0.7"))
    print("*{:^58}*".format(""))
    print("*{:^58}*".format("ENG.WANICK"))
    print("*{:^58}*".format("MÉTODOS NUMÉRICOS"))
    print("*{:^58}*".format("MESTRADO EM ENGENHARIA MECÂNICA"))
    print("*{:^58}*".format("UnB"))
    print("*{:^58}*".format(""))
    print("*" * 60)
    input("Pressione ENTER para continuar...\n")

banner_inicial()

# === FUNÇÕES NUMÉRICAS FUNDAMENTAIS ===
# Cáculo da função seno usando uma expanção da série de Taylor;
def sin_taylor(x, n=20):
    # Redução de argumento para [-π, π]
    pi = 3.141592653589793
    x = x % (2 * pi)
    if x > pi:
        x -= 2 * pi

    termo = x
    soma = x
    sinal = -1
    for i in range(1, n):
        termo *= x * x / ((2 * i) * (2 * i + 1))
        soma += sinal * termo
        sinal *= -1
    return soma

# Cálculo da função logaritmo natural [para a f_5(x)] usando uma expanção da série de Taylor;
def ln_taylor(x, n=100):
    if x <= 0:
        raise ValueError("Logaritmo indefinido para x <= 0")
    z = (x - 1) / x
    termo = z
    soma = termo
    for k in range(2, n + 1):
        termo *= z
        soma += (termo / k) * ((-1) ** (k + 1))
    return soma


# Cálculo do log2(x), quando x for a diferença entre os limites superior e inferior do intervalo;
def ln_estavel(x, n=100):
    ln2=0.6931471805599453 # Valor da cte ln(2);
    if x <= 0:
        raise ValueError("Log indefinido para x <= 0") #Definição de propiedade do logarítmica;
    k = 0
    while x > 2:
        x /= 2
        k += 1
    while x < 0.5:
        x *= 2
        k -= 1
    return k * ln2 + ln_taylor(x, n)


# Cálculo do número de iterações

def calcular_max_iter(xu, xl, tol):
    razao = abs(xu - xl) / tol
    if razao <= 1:
        return 1
    ln2 = 0.6931471805599453 # Valor da cte ln(2);
    ln_aprox = ln_estavel(razao, n=100)
    return int(ln_aprox / ln2) + 1 #Iterações necessárias +1;


# === DEFINIÇÕES DAS FUNÇÕES DO PROBLEMA ===
def f1(x):
    return -0.5 * x**2 + 2.5 * x + 4.5

def f2(x):
    return 5 * x**3 - 5 * x**2 + 6 * x - 2

def f3(x):
    return -25 + 82 * x - 90 * x**2 + 44 * x**3 - 8 * x**4 + 0.7 * x**5

def f4(x):
    return sin_taylor(x) - x**3

def f5(x):
    return ln_taylor(x**4) - 0.7


# === MÉTODO DA BISSECÇÃO ===
def bisseccao(f, xl, xu, tol):
    iteracoes = []
    max_iter = calcular_max_iter(xu, xl, tol)
    print(f"Número máximo de iterações calculado (bissecção): {max_iter}")

    xr_ant = None
    for i in range(1, max_iter + 1):
        xr = (xl + xu) / 2
        fr = f(xr)
        fl = f(xl)
        ea = abs(xr - xr_ant) if xr_ant is not None else abs(xu - xl)
        iteracoes.append((i, xr, fr, ea))
        if fr == 0 or ea < tol:
            break
        elif fl * fr < 0:
            xu = xr
        else:
            xl = xr
        xr_ant = xr
    return xr, iteracoes


# === MÉTODO DA FALSA POSIÇÃO ===
def falsa_posicao(f, xl, xu, tol, max_iter=100):
    iteracoes = []
    xr_ant = None
    for i in range(1, max_iter + 1):
        fl = f(xl)
        fu = f(xu)
        if fl == fu:
            print(f"Divisão por zero evitada: f({xl}) = f({xu}) = {fl}")
            break
        xr = xu - fu * (xl - xu) / (fl - fu)
        fr = f(xr)
        ea = abs(xr - xr_ant) if xr_ant is not None else abs(xu - xl)
        iteracoes.append((i, xr, fr, ea))
        if fr == 0 or ea < tol:
            break
        elif fl * fr < 0:
            xu = xr
        else:
            xl = xr
        xr_ant = xr
    return xr, iteracoes


# === MENU INTERATIVO ===
def menu():
    # Definição de entrada + intervalos pré-definidos;
    funcoes = {
        '1': (f1, 5, 10),
        '2': (f2, 0, 1),
        '3': (f3, 0.5, 1),
        '4': (f4, 0.5, 1),
        '5': (f5, 0.5, 2)
    }

    while True:
        print("\n=== MÉTODOS NUMÉRICOS ===")
        print("Escolha a função para resolver f(x) = 0:")
        print("1 - f1(x) = -0.5x² + 2.5x + 4.5")
        print("2 - f2(x) = 5x³ - 5x² + 6x - 2")
        print("3 - f3(x) = -25 + 82x - 90x² + 44x³ - 8x⁴ + 0.7x⁵")
        print("4 - f4(x) = sin(x) - x³")
        print("5 - f5(x) = ln(x⁴) - 0.7")
        print("0 - Sair")
        escolha = input("Digite sua opção: ")

        if escolha == '0':
            print("Encerrando o programa.")
            break
        elif escolha in funcoes:
            f, a, b = funcoes[escolha]
            tol = float(input("Informe a tolerância desejada (ex: 1e-5): "))

            raiz_bis, tabela_bis = bisseccao(f, a, b, tol)
            raiz_fp, tabela_fp = falsa_posicao(f, a, b, tol, len(tabela_bis))

            print(f"\nRaiz encontrada (Bissecção):     {raiz_bis:.10f}")
            print(f"Raiz encontrada (Falsa Posição): {raiz_fp:.10f}\n")

            print("=== Iterações - Bissecção ===")
            print(f"{'Iter':<6}{'x':<20}{'f(x)':<20}{'Erro':<20}")
            for i, x, fx, e in tabela_bis:
                print(f"{i:<6}{x:<20.10f}{fx:<20.10f}{e:<20.10f}")

            print("\n=== Iterações - Falsa Posição ===")
            print(f"{'Iter':<6}{'x':<20}{'f(x)':<20}{'Erro':<20}")
            for i, x, fx, e in tabela_fp:
                print(f"{i:<6}{x:<20.10f}{fx:<20.10f}{e:<20.10f}")
        else:
            print("Opção inválida. Tente novamente.")

if __name__ == "__main__":
    menu()
