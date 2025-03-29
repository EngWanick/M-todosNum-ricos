::: titlepage
![image](imagens/unb_bandeira.png){width="70mm"}\
[Universidade de Brasília]{.smallcaps}\
[Departamento de Ciências Mecânicas]{.smallcaps}\
[Programa de Pós-Graduação]{.smallcaps}\
**Programa 1**\
**Disciplina: Métodos Numéricos**\
Professor: Rafael Gabler Gontijo\
Data: 2025-03-29\
**Aluno: Eng. Lucas Wanick --- Mestrando em Ciências Mecânicas**\
:::

# Introdução ao Problema {#introdução-ao-problema .unnumbered}

O presente estudo tem como objetivo simular a dinâmica de uma partícula
esferoidal em um fluido viscoso, sob ação da força peso, empuxo e força
de arrasto. A equação do movimento é baseada na segunda lei de Newton e
contempla o modelo de arrasto linear (Stokes). A partícula parte do
repouso e tende à velocidade terminal ao longo do tempo.

![Esquematização do problema proposto pelo
professor.](imagens/esquema_exercicio.png){width="60%"}

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

## Trecho de código (solução analítica via Taylor) {#trecho-de-código-solução-analítica-via-taylor .unnumbered}

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

## Trecho de código (método de Euler) {#trecho-de-código-método-de-euler .unnumbered}

    def euler_step(y, h, St):
        return y + h * (-y + 1) / St

# Análise do Erro

O erro absoluto foi calculado ponto a ponto pela diferença entre a
solução numérica e a analítica:
$$\epsilon(t) = \left| y_{\text{Euler}}(t) - y_{\text{anal\'itico}}(t) \right|$$

# Resultados Gráficos

![Solução analítica via série de Taylor com tolerância e
reescrita.](imagens/grafico_analitico.png){width="90%"}

![Solução numérica via método de Euler (h =
0.01).](imagens/grafico_euler.png){width="90%"}

![Erro absoluto entre a solução analítica e o método de
Euler.](imagens/grafico_erro.png){width="90%"}

# Conclusão {#conclusão .unnumbered}

A resolução do problema proposto demonstrou ser plenamente viável por
meio de implementações numéricas diretas em `Python`, mesmo sem recorrer
a bibliotecas prontas para o cálculo da exponencial. A utilização da
série de Taylor truncada com critério de tolerância mostrou-se eficaz e
computacionalmente estável, desde que aliada a uma reescrita inteligente
da exponencial para evitar problemas numéricos associados à
representação em ponto flutuante (`float64`). A exportação dos dados
simulados para um arquivo `.csv` permitiu a posterior geração de
gráficos, favorecendo a visualização dos resultados e a comparação entre
as soluções analítica e numérica. O estudo evidencia a robustez da
abordagem proposta e reforça a aplicabilidade de métodos computacionais
simples na modelagem de fenômenos físicos com eficiência e controle.
