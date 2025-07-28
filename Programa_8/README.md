::: titlepage
![image](img/unb_bandeira.png){width="12cm"}\
[Universidade de Brasília]{.smallcaps}\
[Departamento de Ciências Mecânicas]{.smallcaps}\
[Programa de Pós-Graduação]{.smallcaps}\
**Programa 8**\
**Equações Diferenciais**\
**Tranferência de Calor em Regime Transiente**\
**Geometria Complexa**\
**Disciplina: Métodos Numéricos**\
Professor: Dr. Rafael Gabler Gontijo\
**Aluno: Eng. Lucas Wanick --- Mestrando em Engenharia Mecânica**\
2025-07-28\
:::

# Introdução

O presente relatório documenta o desenvolvimento do último programa do
curso de *Métodos Numéricos* ministrado pelo Professor Rafael Gabler
Gontijo na Universidade de Brasília. O objetivo central foi implementar
um solver numérico capaz de resolver a equação de condução de calor
estacionária em uma geometria retangular com cantos curvos
(adiabáticos), submetida a condições de contorno mistas:

![Enunciado do problema - Lousa da Aula
35.](img/enunciado.png){width="80%"}

- Temperaturas fixas (Dirichlet) nas faces superior e inferior;

- Fluxos de calor prescritos (Neumann) nas faces laterais;

- Isolamento térmico (condição adiabática) nos cantos arredondados.

O desenvolvimento visou obter uma solução robusta, eficiente e
escalável, capaz de lidar com malhas refinadas , preservando a
fidelidade geométrica e garantindo estabilidade numérica.

![Geometria do problema. Dimenções arbitradas: $L = 1.50\,\mathrm{m}$,
$H = 2.50\,\mathrm{m}$, raio dos cantos
$R = 0.25\,\mathrm{m}$.](img/geometria.png){width="60%"}

# 2. Formulação Matemática {#formulação-matemática .unnumbered}

A equação governante é a equação de Laplace bidimensional para regime
estacionário ($\nabla^2T=0$):

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

## 2.1 Condições de contorno {#condições-de-contorno .unnumbered}

- **Dirichlet (superior e inferior):**
  $$T(y=0) = 45^\circ C, \qquad T(y=H) = 55^\circ C$$

- **Neumann (laterais):**
  $$-k \frac{\partial T}{\partial x}\bigg|_{x=0} = q_{\text{left}}, 
      \qquad 
      -k \frac{\partial T}{\partial x}\bigg|_{x=L} = q_{\text{right}}$$
  onde $k$ é a condutividade térmica do material, definida como 71
  $W/mK$ (Aço) e $q_{left}$, $q_{right}$ são os fluxos de calor
  prescritos de 250 $W/m^2$ e 210 $W/m^2$, respectivamente.

- **Adiabático (cantos):** $$\frac{\partial T}{\partial n} = 0$$

Geometria: $L = 1.50\,\mathrm{m}$, $H = 2.50\,\mathrm{m}$, cantos curvos
de raio $R = 0.25\,\mathrm{m}$.

# 3. Discretização Numérica {#discretização-numérica .unnumbered}

O domínio foi discretizado em uma malha uniforme de $n \times n$ pontos,
utilizando o método das diferenças finitas de 5 pontos:

$$-4T_{i,j} + T_{i+1,j} + T_{i-1,j} + T_{i,j+1} + T_{i,j-1} = 0$$

As condições de contorno foram incorporadas da seguinte forma:

- **Dirichlet:** valores fixos diretamente no vetor fonte $b$.

- **Neumann:** usando a formulação ajustada:
  $$T_{\text{parede}} = T_{\text{vizinho}} + \frac{q'' \, \Delta x}{k}$$
  com compensação na matriz: $$A_{kk} \leftarrow A_{kk} + 1$$

- **Adiabático:** regiões de contorno curvo foram identificadas por
  verificação geométrica $(x,y)$ e excluídas do sistema.

# 4. Implementação Computacional {#implementação-computacional .unnumbered}

- Linguagem: Python 3.12

- Bibliotecas: `numpy`, `scipy.sparse` (LIL $\rightarrow$ CSR),
  `matplotlib`

Estratégias de desempenho:

- Matriz esparsa, evitando alocação densa (redução de memória de TB
  $\rightarrow$ MB).

- Solver `spsolve` para eficiência em sistemas grandes.

- Máscaras (`NaN`) para lidar com cantos adiabáticos.

## 4.1 Geração da malha e classificação dos nós {#geração-da-malha-e-classificação-dos-nós .unnumbered}

A função `gerar_malha_tipo` atribui a cada nó:

- **I:** interno

- **D:** Dirichlet

- **N:** Neumann

- **A:** adiabático

![Malha gerada com classificação dos
nós.](img/Figure_2.png){width="80%"}

## 4.2 Montagem do sistema {#montagem-do-sistema .unnumbered}

- Mapeamento de incógnitas internas.

- Montagem de $A$ e $b$ considerando as diferentes condições.

- Implementação rigorosa da Neumann com compensação na diagonal.

## 4.3 Interpolação para visualização {#interpolação-para-visualização .unnumbered}

Pós-processamento para interpolar valores em regiões Neumann e
adiabáticas, garantindo visualização contínua sem alterar a solução
numérica.

# 5. Resultados Numéricos {#resultados-numéricos .unnumbered}

## 5.1 Distribuição de temperatura {#distribuição-de-temperatura .unnumbered}

O solver convergiu para distribuições estáveis em malhas de até
$729 \times 729$. Temperaturas coerentes com as condições impostas:

- Regiões centrais mais quentes devido ao aporte de calor lateral.

- Gradiente vertical suave entre $45^\circ C$ e $55^\circ C$.

![Distribuição de temperatura calculada. Valores interpolados para
visualização contínua.](img/Figure_1.png){width="80%"}

## 5.2 Comportamento físico {#comportamento-físico .unnumbered}

- Alteração do sinal de $q_{\text{left}}$ ou $q_{\text{right}}$ produz
  gradientes opostos, confirmando implementação correta da Neumann.

- Perfis de fluxo e isolamento térmico nos cantos respeitados.

- Balanço energético coerente: fluxo líquido $\rightarrow$ aumento da
  temperatura média.

# 6. Conclusão {#conclusão .unnumbered}

O programa desenvolvido cumpre plenamente os objetivos:

- Resolver a equação de condução de calor em geometria complexa com
  múltiplos tipos de contorno.

- Suportar malhas muito refinadas com alta eficiência de memória.

- Produzir resultados fisicamente coerentes e visualmente
  interpretáveis.

Este projeto encerra o curso com demonstração de proficiência em métodos
numéricos, implementação computacional avançada e capacidade analítica
de validar soluções numéricas frente à realidade física.

# 7. Trabalhos futuros {#trabalhos-futuros .unnumbered}

- Implementação de parâmetros $\alpha$ e $\beta$ para ajuste fino da
  curvatura em malhas grosseiras.

- Extensão para problemas transientes (dependentes do tempo).

- Acoplamento com modelos de convecção e radiação para simular trocas
  térmicas complexas.
