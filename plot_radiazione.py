import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
import projeto2 as old

def A():
    return 0.4

def w():
    return np.pi

def posicao_y_da_carga(t):
        return A()*np.cos(w()*t)

def posicao_x_da_carga(t):
        return A()*np.cos(w()*t)

def plot_rad(potencial):
    #Campo baseado no potencial antigo    
    Ey = np.roll(potencial, 1, axis=0) - np.roll(potencial, -1, axis=0) 
    Ex = np.roll(potencial, 1, axis=1) - np.roll(potencial, -1, axis=1)
    
    #Coordenadas
    x = np.linspace(-1, 1, old.n_malha())
    y = np.linspace(-1, 1, old.n_malha())
    X, Y = np.meshgrid(x, y)    
    
    figure, ax = plt.subplots()
    figure.set_size_inches((7,7))
    camera = Camera(figure)
    for t in range(20):
        Ry = Y - posicao_y_da_carga(t/10)
        Rx = X - posicao_x_da_carga(t/10)
        radiation = np.abs(Ry*Ey + Rx*Ex)
        
        ax.pcolor(radiation, cmap=plt.cm.inferno_r)
        camera.snap()    
    animation = camera.animate()
    animation.save('radiazione.gif', writer = 'imagemagick')
    
casos = {'quadrado': 0, 'circulo':1, 'capacitor': 2, 'barra': 3}
key = 'circulo'
caso = casos[key]

potencial = old.diferenca_finita(caso)

#No caso essa função está afora recebendo o potencial do projeto antigo
#Precisamos atualizar pro caso novo dispois
plot_rad(potencial)