

% Capa
titlepage
    \centering
    \includegraphics[width=70mm]{imagens/unb_bandeira.png} \\
    1cm
    \Large Universidade de Brasília \\
    Departamento de Ci\^encias Mec\^anicas \\
    Programa de Pós-Graduaç\~ao \\
    \vfill
    {\Large\bfseries Programa 1} \\
    0.5cm
    Disciplina: Métodos Numéricos \\
    Professor: Rafael Gabler Gontijo \\
    Data: \today \\
    \vfill
    Aluno: Eng. Lucas Wanick — Mestrando em Ci\^encias Mec\^anicas \\
    0.5cm
titlepage

# Introduç\~ao ao Problema
O presente estudo tem como objetivo simular a din\^amica de uma partícula esferoidal em um fluido viscoso, sob aç\~ao da força peso, empuxo e força de arrasto. A equaç\~ao do movimento é baseada na segunda lei de Newton e contempla o modelo de arrasto linear (Stokes). A partícula parte do repouso e tende à velocidade terminal ao longo do tempo. 

figure[H]
    \centering
    \includegraphics[width=0.6\textwidth]{imagens/esquema_exercicio.png}
    Esquematizaç\~ao do problema proposto pelo professor.
figure

# Modelagem Matemática e Soluç\~ao Analítica
O ponto de partida para a modelagem do problema \'{e} a aplicaç\~ao da Segunda Lei de Newton à partícula em movimento dentro de um fluido viscoso. Considerando a força peso, a força de empuxo e a força de arrasto linear (regime de Stokes), a equaç\~ao do movimento \'{e} descrita por:

$$
    y(t) = 1 - e^{-t/St}
$$

A funç\~ao exponencial foi implementada manualmente utilizando a série de Taylor:

$$
    e^x = \sum_{n=0}^{\infty} x^n{n!}
$$

Para manter a estabilidade numérica e evitar estouros por overflow ou underflow, utilizamos uma abordagem baseada na identidade:

$$
    e^{-x} = 1{e^x}
$$

A série foi truncada automaticamente com base em um critério de toler\^ancia de $10^{-15}$, o que assegura a converg\^encia da soluç\~ao dentro de uma precis\~ao aceitável.

\newpage
## Trecho de código (soluç\~ao analítica via Taylor)
```python
def exp_taylor(x, n=100, tol=1e-15):
    soma = 0.0
    termo = 1.0
    for i in range(n):
        soma += termo
        termo *= x / (i + 1)
        if abs(termo) < tol:
            break
    return soma

def exp_negativo(x):
    return 1 / exp_taylor(x)

def analytical(t, St):
    return 1 - exp_negativo(t / St)
```

# Implementaç\~ao do Método de Euler
O método de Euler foi implementado explicitamente como aproximaç\~ao da derivada por diferença progressiva:

$$
    y_{n+1} = y_n + h \cdot -y_n + 1{St}
$$

O passo de integraç\~ao adotado foi $h = 0{,}01$ para todos os casos, garantindo boa estabilidade e precis\~ao numérica.

Utilizar uma abordagem numérica em paralelo à soluç\~ao analítica é fundamental para validar a consist\^encia da simulaç\~ao, principalmente em contextos onde a soluç\~ao exata n\~ao está disponível ou onde se deseja verificar a robustez de métodos aproximativos.

## Trecho de código (método de Euler)
```python
def euler_step(y, h, St):
    return y + h * (-y + 1) / St
```

# Análise do Erro
O erro absoluto foi calculado ponto a ponto pela diferença entre a soluç\~ao numérica e a analítica:
$$
    \epsilon(t) = \left| y_{Euler}(t) - y_{analítico}(t) \right|
$$

# Resultados Gráficos
-1em
figure[H]
    \centering
    \includegraphics[width=0.9\textwidth]{imagens/grafico_analitico.png}
    Soluç\~ao analítica via série de Taylor com toler\^ancia e reescrita.
figure

figure[H]
    \centering
    \includegraphics[width=0.9\textwidth]{imagens/grafico_euler.png}
    Soluç\~ao numérica via método de Euler (h = 0.01).
figure

figure[H]
    \centering
    \includegraphics[width=0.9\textwidth]{imagens/grafico_erro.png}
    Erro absoluto entre a soluç\~ao analítica e o método de Euler.
figure

# Conclusão  
A resolução do problema proposto demonstrou ser plenamente vi\'{a}vel por meio de implementações num\'{e}ricas diretas em Python, mesmo sem recorrer a bibliotecas prontas para o c\'{a}lculo da exponencial. A utilização da s\'{e}rie de Taylor truncada com crit\'{e}rio de toler\^{a}ncia mostrou-se eficaz e computacionalmente est\'{a}vel, desde que aliada a uma reescrita inteligente da exponencial para evitar problemas num\'{e}ricos associados \`{a} representação em ponto flutuante (float64). A exportação dos dados simulados para um arquivo .csv permitiu a posterior geração de gr\'{a}ficos, favorecendo a visualização dos resultados e a comparação entre as soluções anal\'{i}tica e num\'{e}rica. O estudo evidencia a robustez da abordagem proposta e reforça a aplicabilidade de m\'{e}todos computacionais simples na modelagem de fen\^omenos f\'{i}sicos com efici\^encia e controle.

