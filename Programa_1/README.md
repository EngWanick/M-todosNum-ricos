![image](imagens/unb_bandeira.png)\
# Universidade de Brasília
**Departamento de Ciências Mecânicas**\
**Programa de Pós-Graduação**\
**Programa 1**\
**Disciplina: Métodos Numéricos**\
Professor: Rafael Gabler Gontijo\
Data: 2025-03-29\
**Aluno: Eng. Lucas Wanick --- Mestrando em Ciências Mecânicas**\


# Introdução ao Problema

O presente estudo tem como objetivo simular a dinâmica de uma partícula
esferoidal em um fluido viscoso, sob ação da força peso, empuxo e força
de arrasto. A equação do movimento é baseada na segunda lei de Newton e
contempla o modelo de arrasto linear (Stokes). A partícula parte do
repouso e tende à velocidade terminal ao longo do tempo.

![Esquematização do problema proposto pelo
professor.](imagens/esquema_exercicio.png)

# Modelagem Matemática e Solução Analítica

O ponto de partida para a modelagem do problema é a aplicação da Segunda
Lei de Newton à partícula em movimento dentro de um fluido viscoso.
Considerando a força peso, a força de empuxo e a força de arrasto linear
(regime de Stokes), a equação do movimento é descrita por:

$$y(t) = 1 - e^{-t/St}$$

A função exponencial foi implementada manualmente utilizando a série de
Taylor:

$$e^x = \sum_{n=0}^{\infty} \frac{x^n}{n!}$$

Para manter a estabilidade numérica e evitar estouros por overflow ou
underflow, utilizamos uma abordagem baseada na identidade:

$$e^{-x} = \frac{1}{e^x}$$

A série foi truncada automaticamente com base em um critério de
tolerância de $10^{-15}$, o que assegura a convergência da solução
dentro de uma precisão aceitável.

## Trecho de código (solução analítica via Taylor)

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

# Implementação do Método de Euler

O método de Euler foi implementado explicitamente como aproximação da
derivada por diferença progressiva:

$$y_{n+1} = y_n + h \cdot \frac{-y_n + 1}{St}$$

O passo de integração adotado foi $h = 0{,}01$ para todos os casos,
garantindo boa estabilidade e precisão numérica.

Utilizar uma abordagem numérica em paralelo à solução analítica é
fundamental para validar a consistência da simulação, principalmente em
contextos onde a solução exata não está disponível ou onde se deseja
verificar a robustez de métodos aproximativos.

## Trecho de código (método de Euler)
    def euler_step(y, h, St):
        return y + h * (-y + 1) / St

# Análise do Erro

O erro absoluto foi calculado ponto a ponto pela diferença entre a
solução numérica e a analítica:

$$\epsilon(t) = \left| y_{\text{Euler}}(t) - y_{\text{anal\'itico}}(t) \right|$$

# Resultados Gráficos

![Solução analítica via série de Taylor com tolerância e
reescrita.](imagens/grafico_analitico.png)

![Solução numérica via método de Euler (h =
0.01).](imagens/grafico_euler.png)

![Erro absoluto entre a solução analítica e o método de
Euler.](imagens/grafico_erro.png)

# Conclusão

A resolução do problema evidenciou a coerência entre a solução analítica
 e a aproximação numérica obtida pelo método de Euler. A forma analítica,
 construída a partir da equação diferencial, apresentou um comportamento
 esperado para sistemas com regime de arrasto linear. A aproximação exponencial
 por série de Taylor demonstrou-se estável, com erro controlado mediante truncamento
 adequado.


