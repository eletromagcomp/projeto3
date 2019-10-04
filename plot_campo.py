#%%RECOMENDAÇÕES
#Baixem o celluloid no terminal "pip install celluloid"
#Estou usando o codigo do projeto passado, mas atualizando com o potencial com variação temporal vai ser bem tranquilo
#Qualquer duvida perguntem, sabado eu atualizo com o projeto de verdade
#Ainda falta tbm adicionar titulo e tals

#%% IMPORT

import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera

#%% VARIÁVEIS
def n_malha():
    return 10

def epsilon():
    return 0.00001

def tamanho_solido():
    return 0.2

def pot_solido():
    return 100

#%% DEFINIÇÃO DA MALHA E DAS CONDIÇÕES DE CONTORNO
def condutor(caso):
    n = n_malha()
    pot_sol = pot_solido()
    condutor_bool = np.zeros((n,n), dtype=bool)
    potencial = np.ones((n,n))
    
    #Quadrado
    if caso==0:
        solido_lado = int(n*tamanho_solido())
        solido = np.arange(int((n - solido_lado)/2), int((n + solido_lado)/2))
    
        condutor_bool[0, :] = True
        condutor_bool[n-1, :] = True
        condutor_bool[:, 0] = True
        condutor_bool[:, n-1] = True

        potencial[0,:] = 0
        potencial[n-1,:] = 0
        potencial[:,0] = 0
        potencial[:,n-1] = 0
       
        for i in solido:
           condutor_bool[i, solido] = True
           potencial[i, solido] = pot_sol
    
    #Círculo
    if caso==1:
        raio = n*tamanho_solido()/2 #raio do solido
        raio_ext = n/2            #raio da superficie
        solido = np.arange(int(n/2-raio),int(n/2+raio))
        
        for i in range(n):
            for j in range(n):
                if( (i-n/2)**2+(j-n/2)**2 >= raio_ext**2):
                    condutor_bool[i,j] = True
                    potencial[i,j] = 0
                    
        for i in solido:
            for j in solido:
                if( (i-n/2)**2+(j-n/2)**2 <= raio**2):
                    condutor_bool[i,j] = True
                    potencial[i,j] = pot_sol
                    
    return potencial, condutor_bool  

#%% MALHA DOS VIZINHOS
def vizinhos(potencial):
    #Esquerda
    potencial_esquerda = np.roll(potencial, 1, axis=1)
    #Direita
    potencial_direita = np.roll(potencial, -1, axis=1)
    #Em baixo
    potencial_baixo = np.roll(potencial, -1, axis=0)
    #Em cima
    potencial_cima = np.roll(potencial, 1, axis=0)
    return potencial_esquerda, potencial_direita, potencial_baixo, potencial_cima

#%% POTENCIAL NUMÉRICO - RESOLUÇÃO DA EQUAÇÃO DE LAPLACE
def diferenca_finita(caso):
    potencial, condutor_bool = condutor(caso)
    eps = epsilon()
    diferenca_max = eps + 1
    while diferenca_max>=eps:
        potencial_esquerda, potencial_direita, potencial_baixo, potencial_cima = vizinhos(potencial)
        potencial_novo = np.where(condutor_bool == False, 1/4 * (potencial_esquerda + potencial_direita
                                                             + potencial_cima + potencial_baixo), potencial)
        diferenca = np.absolute(potencial_novo - potencial)
        diferenca_max = np.amax(diferenca)
        potencial = potencial_novo
    return potencial

#%% PLOT - POTENCIAL, CAMPO ELÉTRICO E EQUIPOTENCIAIS
def plot_campo(potencial, levels=10, linewidth=1, density=0.5,
               arrowsize=1.5, surface_label = False, fig='',
               fig1_name='Mapa_de_Cor.png', fig2_name='Campo.png'):
    
    figure, ax = plt.subplots()
    pot = [[],[]]
    pot[0] = potencial
    pot[1] = potencial
    figure.set_size_inches((7,7))
    camera = Camera(figure)
    t = (0,1)
    print(t)
    for i in t:
        
        #Equipotenciais e linhas de campo
        x = np.linspace(-1, 1, n_malha())
        y = np.linspace(-1, 1, n_malha())
        X, Y = np.meshgrid(x, y)
            
        Ex = np.roll(pot[int(i)], 1, axis=0) - np.roll(pot[int(i)], -1, axis=0) 
        Ey = np.roll(pot[int(i)], 1, axis=1) - np.roll(pot[int(i)], -1, axis=1)
    
        #CS = ax.contour(X, Y, pot[int(i)], cmap=plt.cm.inferno, levels = levels)
        #if surface_label: plt.clabel(CS, fontsize=14)
    
        color = 2 * (np.hypot(Ex, Ey))**(1/2)
        ax.streamplot(y, x, Ey, Ex, color=color, linewidth=linewidth, cmap=plt.cm.inferno, 
                      density=density, arrowstyle='->', arrowsize=arrowsize)
    #    ax.set_xlabel('$x$')
    #    ax.set_ylabel('$y$')
        ax.set_xlim(-1.05,1.05)
        ax.set_ylim(-1.05,1.05)
        ax.set_aspect('equal')
        camera.snap()
    animation = camera.animate()
    animation.save('celluloid_minimal.gif', writer = 'imagemagick')

#%% CÁLCULOS

#casos = {'quadrado': 0, 'circulo':1}

potencial = diferenca_finita(0)
print(potencial)

#potencial_1d = potencial[100,100:199]
plot_campo(potencial, surface_label=True, density=1, fig="oi" + '_')
