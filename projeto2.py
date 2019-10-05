#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:39:05 2019
@author: hiro
"""
import numpy as np
import matplotlib.pyplot as plt
import time

#testando a bagaca no terminal

#%% VARIÁVEIS
def n_malha():
    return 100

def epsilon():
    return 0.00001

def tamanho_solido():
    return 0.2
    
def largura_capacitor():
    return 0.2
    
def distancia_placas():
    return 0.1
    
def pot_placa():
    return 100

def pot_solido():
    return 100

def caminhantes_():
    return 10

#%% DEFINIÇÃO DA MALHA E DAS CONDIÇÕES DE CONTORNO
def condutor(caso):
    n = n_malha()
    pot_sol = pot_solido()
    V = pot_placa()
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
                    
    #Capacitor
    if caso==2:
        largura = int(n*largura_capacitor())
        altura = int(n*distancia_placas()/2)
        placa = np.arange(int((n - largura)/2), int((n + largura)/2))
    
        condutor_bool[0, :] = True
        condutor_bool[n-1, :] = True
        condutor_bool[:, 0] = True
        condutor_bool[:, n-1] = True

        potencial[0,:] = 0
        potencial[n-1,:] = 0
        potencial[:,0] = 0
        potencial[:,n-1] = 0
        
        condutor_bool[int(n/2+altura), placa] = True
        condutor_bool[int(n/2-altura), placa] = True
        potencial[int(n/2+altura), placa] = V
        potencial[int(n/2-altura), placa] = -1*V
    
    #Barra
    if caso==3:
        condutor_bool[0, :] = True
        condutor_bool[n-1, :] = True
        condutor_bool[:, 0] = True
        condutor_bool[:, n-1] = True
        
        potencial[0,:] = 0
        potencial[n-1,:] = pot_sol
        potencial[:,0] = 0
        potencial[:,n-1] = 0
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

#%% POTENCIAL ANALÍTICO (BARRA)
def potencial_analitico_barra(caso):
    eps = epsilon()
    n = n_malha()
    x = np.arange(n)
    y = np.arange(n)
    X, Y = np.meshgrid(x,y)
    dif_max = eps +1
    potencial = np.zeros(X.shape)
    if caso==3:
        pot_sol = pot_solido()
        i = 1
        while dif_max>=eps/1000:
            potencial_novo = 4*pot_sol/((2*i - 1)*np.pi) *np.sin(((2*i - 1) *np.pi * X)/n)*(np.exp((2*i -1)*np.pi*(Y/n - 1)) - np.exp(-(2*i -1)*np.pi*(Y/n + 1)))/(1- np.exp(-2**(2*n-1)*np.pi)) + potencial
            dif = np.absolute(potencial_novo - potencial)
            dif_max = np.amax(dif)
            potencial = potencial_novo
            i = i +1
    return potencial

#%% POTENCIAL ANALÍTICO (circulo)
def potencial_analitico_circulo():
    radius = int(n_malha()/2)
    condutor = int(n_malha()*tamanho_solido()/2)-1
    
    potencial = np.zeros(radius)
    
    for i in range(condutor):
        potencial[i] = pot_solido()
    for i in range(condutor,radius):
        potencial[i] = pot_solido()*np.log(i/radius)/np.log(condutor/radius)
    return potencial

#%% PLOT - POTENCIAL, CAMPO ELÉTRICO E EQUIPOTENCIAIS
def plot_campo(potencial, levels=10, linewidth=1, density=0.5,
               arrowsize=1.5, surface_label = False, fig='',
               fig1_name='Mapa_de_Cor.png', fig2_name='Campo.png'):
    
    #Potencial
    plt.figure(figsize=(7, 6))
    plt.pcolor(potencial)
    plt.colorbar()
#    plt.show()
    plt.savefig(fig+fig1_name, dpi=200, bbox_inches='tight')
    
    #Equipotenciais e linhas de campo
    x = np.linspace(-1, 1, n_malha())
    y = np.linspace(-1, 1, n_malha())
    X, Y = np.meshgrid(x, y)
        
    Ex = np.roll(potencial, 1, axis=0) - np.roll(potencial, -1, axis=0) 
    Ey = np.roll(potencial, 1, axis=1) - np.roll(potencial, -1, axis=1)
    
    figure, ax = plt.subplots()
    figure.set_size_inches((7,7))

    CS = ax.contour(X, Y, potencial, cmap=plt.cm.inferno, levels = levels)
    if surface_label: plt.clabel(CS, fontsize=14)

    color = 2 * (np.hypot(Ex, Ey))**(1/2)
    ax.streamplot(y, x, Ey, Ex, color=color, linewidth=linewidth, cmap=plt.cm.inferno, 
                  density=density, arrowstyle='->', arrowsize=arrowsize)
#    ax.set_xlabel('$x$')
#    ax.set_ylabel('$y$')
    ax.set_xlim(-1.05,1.05)
    ax.set_ylim(-1.05,1.05)
    ax.set_aspect('equal')
#    plt.show()
    plt.savefig(fig+fig2_name, dpi=200, bbox_inches='tight')

#%% POTENCIAL EM COORDENADAS POLARES
def potencial_polar(caso=1):
    n = n_malha()/2
    tamanho_condutor = int(n*tamanho_solido())
    resto_da_malha = int(n*(1-tamanho_solido()))
    
    potencial = np.append( pot_solido() * np.ones(tamanho_condutor), np.zeros(resto_da_malha) )
    
    eps = epsilon()
    diferenca_max = eps + 1
    
    if caso == 0: # Gauss-Seidel
        potencial_novo = np.copy(potencial)
        while diferenca_max >= eps:
            for i in range(tamanho_condutor, int(n)-1):
                potencial_novo[i] = ((i-1/2)*potencial_novo[i-1] + (i+1/2)*potencial_novo[i+1]) / (2*i)
            
            diferenca_max = np.max( np.absolute(potencial_novo - potencial) )
            potencial = np.copy(potencial_novo)
            
    if caso == 1: # Busca direta
        potencial[tamanho_condutor] = 0.9*pot_solido() #chute
        passo = 0.1*pot_solido()
        
        while diferenca_max >= eps:
            for i in range(tamanho_condutor, int(n)-1):
                potencial[i+1] = (2*i*potencial[i] - (i-1/2)*potencial[i-1]) / (i+1/2)
            
            diferenca_max = np.absolute(potencial[int(n-1)])   
            if potencial[int(n-1)] < 0:
                passo = passo/2
                potencial[tamanho_condutor] += passo
            elif potencial[int(n-1)] > 0:
                passo = passo/2
                potencial[tamanho_condutor] += -passo
            else:         
                diferenca_max = 0
    
    return potencial

            
#%% PLOT POLAR
def plot_polar(potencial, levels=8, linewidth=1, density=0.5,
               arrowsize=1.5, surface_label = True, fig='',
               fig1_name='Mapa_de_Cor_Polar.png', fig2_name='Campo_Polar.png'):
    #Variáveis theta e raio
    t = np.linspace(0, 2*np.pi, n_malha()/2)
    r = np.linspace(0, 1, n_malha()/2)
    
    T, potencial_2d = np.meshgrid(t, potencial)
    
    #Plot Mapa de Cor
    figure, ax = plt.subplots(1, subplot_kw=dict(projection='polar'))
    figure.set_size_inches((7,6))
    plot = ax.contourf(t, r, potencial_2d, 100, cmap=plt.cm.viridis)
    plt.colorbar(plot, ax=ax)
    #ax.grid(False)
    #plt.show()
    plt.savefig(fig+fig1_name, dpi=200, bbox_inches='tight')
    
    #Cálculo do Campo
    Er = np.roll(potencial, 1) - potencial
    Er[0] = 0
    
    T, Er = np.meshgrid(t, Er)
    a = np.zeros_like(Er)
    
    figure, ax = plt.subplots(1, subplot_kw=dict(projection='polar'))
    figure.set_size_inches((7,6))
    
    CS = ax.contour(t, r, potencial_2d, cmap=plt.cm.inferno, levels = levels)
    if surface_label: plt.clabel(CS, fontsize=14)
    
    color = 2 * np.sqrt(Er)
    plot = ax.streamplot(t, r, a, Er, color=color, linewidth=linewidth, cmap=plt.cm.inferno, 
                  density=density, arrowstyle='->', arrowsize=arrowsize)
    
    ax.set_ylim(0, 1.04)
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.spines['polar'].set_visible(False)
    #plt.show()
    plt.savefig(fig+fig2_name, dpi=200, bbox_inches='tight')
    
#%% RANDOM WALK

def random_walk(caso):
    #Passando informações de contorno
    potencial, extremidades = condutor(caso)
    n = n_malha()
    n_walkers = caminhantes_()
    #Andando:
    variacoes = [[1,0],[-1,0],[0,1],[0,-1]]
    for i in range(int(n/4),int(3*n/4+1)):
        for j in range(int(n/4),int(3*n/4+1)):
            p = 0
            potencial_b = np.zeros(1)
            while p < n_walkers:    
                lugar = np.array([i,j])
                if extremidades[lugar[0],lugar[1]] != True:
                    while extremidades[lugar[0],lugar[1]] != True:
                        lugar = lugar + variacoes[int(np.random.randint(4, size=1))]    
                    potencial_b = np.append(potencial_b, potencial[lugar[0]][lugar[1]])
                if extremidades[lugar[0],lugar[1]] == True:
                    potencial_b = np.append(potencial_b, potencial[lugar[0]][lugar[1]])
                p = p + 1
            potencial[i,j] = (np.sum(potencial_b))/n_walkers
        print('Acabou linha ' + str(i))
    return potencial