# Configurações gerais do Gnuplot
    #gnuplot_script = """
set terminal pngcairo size 1000,800 enhanced font 'Arial,10'
set output 'mapa_bairstow.png'

set title "Fractal de Bairstow - Número de Iterações até Convergência"
set xlabel "r_0 (chute inicial)"
set ylabel "s_0 (chute inicial)"
set cblabel "Número de Iterações"

set palette defined ( \
  0 "#440154", \
  1 "#482777", \
  2 "#3E4989", \
  3 "#31688E", \
  4 "#26828E", \
  5 "#1F9E89", \
  6 "#35B779", \
  7 "#6CCE59", \
  8 "#B4DE2C", \
  9 "#FDE725" \
)

set view map
set datafile missing "999"

# Plotagem
splot 'bairstow_iteracoes.txt' using 1:2:3 with image notitle
"""
    with open("fractal_bairstow.gnu", "w") as f:
        f.write(gnuplot_script)
    print("Script Gnuplot salvo em fractal_bairstow.gnu")