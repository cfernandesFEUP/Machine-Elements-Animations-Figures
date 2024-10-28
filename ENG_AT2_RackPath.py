## INPUT ##
alpha = 20#float(input("Ângulo de Pressão:"))
haP = 1.25#float(input("Altura de Cabeça da Cremalheira = Altura de Pé da Roda Dentada:"))
hfP = 1.#float(input("Altura de Pé da Cremalheira = Altura de Cabeça da Roda Dentada:"))
roh = 0.38#float(input("Raio da Ferramenta:"))
m = float(input("Módulo:"))
z = float(input("Número de dentes:"))
xdesv = float(input("Correcção de Dentado:"))
## DISCRETIZAÇÃO ##
nl = 30
n = 1000
## FUNÇÃO ##
import figures.functions.rack as rack
# from matplotlib.animation import FileMovieWriter
ani = rack.animation(alpha, haP, hfP, roh, m, z, xdesv, nl, n)

#import figures.functions.saveframes as sf
#ani.save('figures/pdf/anim/geracao.pdf', writer=sf.BunchOFiles())