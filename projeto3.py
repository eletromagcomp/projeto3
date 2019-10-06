#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 20:03:28 2019
@author: hiro
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from celluloid import Camera
from scipy.optimize import newton
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pdb

#%% VARIÁVEIS
def n_malha():
    return 30

def amplitude():
    return 2

def v_max(): #Entre 0 e 1
    return 0.8

def time_step():
    return 1

def frequencia():
    return v_max()/amplitude()

def periodo():
    return 2*np.pi/frequencia()
#%% FUNÇÃO DA POSIÇÃO DA CASA
def charge_position(t_ret):
    x = 0.5
    y = 0
    z = amplitude()*np.cos(t_ret*frequencia())
    xyz_charge = np.array((x, y, z))
    return xyz_charge

#%% FUNÇÃO DA VELOCIDADE DA CARGA
def charge_velocity(t_ret):
    v_x = 0
    v_y = 0
    v_z = -v_max() * np.sin(t_ret*frequencia())
    v_charge = np.array((v_x, v_y, v_z))
    return v_charge
#%% FUNÇÃO DA ACELERAÇÃO DA CARGA
def charge_acceleration(t_ret):
    a_x = 0
    a_y = 0
    a_z = -(v_max()**2)/amplitude() * np.cos(t_ret*frequencia())
    a_charge = np.array((a_x, a_y, a_z))
    return a_charge

#%% CAMPO ELETRICO
def electric_field(t_ret, xyz):
    xyz_charge = charge_position(t_ret)
    R = xyz - xyz_charge
    R_norm = np.linalg.norm(R)
    R_hat = R/R_norm
        
    v_charge = charge_velocity(t_ret)
    v_charge_norm = np.linalg.norm(v_charge)
    a_charge = charge_acceleration(t_ret)
    u = R_hat - v_charge
    
    electric = R_norm/(np.dot(R, u))**3 * ((1 - v_charge_norm**2)*u + np.cross(R, np.cross(u, a_charge)))
    return electric
    
#%% SIMULAÇÃO
def simulation(t):
    def t_retarded(t_ret,*args):
        x,y,z = xyz
        dist = np.sqrt(((x - 0.5)**2 + y**2 + (z - amplitude()*np.cos(frequencia()*t_ret))**2))
        return dist - (t - t_ret)

    def t_retarded_prime(t_ret, *args):
        x,y,z = xyz
        dist = np.sqrt(((x - 0.5)**2 + y**2 + (z - amplitude()*np.cos(frequencia()*t_ret))**2))
        return (v_max()*np.sin(t_ret*frequencia())*(z - amplitude()*np.cos(t_ret*frequencia())))/dist + 1

    x = np.arange(-n_malha()/2, n_malha()/2)
    y = 0
    z = np.arange(-n_malha()/2, n_malha()/2)
    
    electric_x = np.zeros((n_malha(), n_malha()))
    electric_y = np.zeros((n_malha(), n_malha()))
    electric_z = np.zeros((n_malha(), n_malha()))
    
    guess = 1
    for i in range(n_malha()):
        for j in range(n_malha()):
            xyz = np.array((x[i], y, z[j]))
            t_ret = newton(t_retarded, guess, fprime=t_retarded_prime, args=(xyz, t), tol=10**(-3), maxiter = 100)
            guess = t_ret
            electric_x[i, j], electric_y[i, j], electric_z[i, j] = electric_field(t_ret, xyz)
    electric_xyz = np.array([electric_x, electric_y, electric_z])
    return  electric_xyz
#%% PLOT
def plot(electric_xyz, t):
    x = np.arange(-n_malha()/2, n_malha()/2)
    z = np.arange(-n_malha()/2, n_malha()/2)
    X, Z = np.meshgrid(x, z)
    x_charge, y_charge, z_charge = charge_position(t)
    electric_x, electric_y, electric_z = electric_xyz
    
    #Carga
    ax_stream.plot(x_charge, z_charge, 'bo')
    ax_quiver.plot(x_charge, z_charge, 'bo')
    ax_both.plot(x_charge, z_charge, 'bo')
    ax_radiation.plot(x_charge+ n_malha()/2 + 0.5, z_charge + n_malha()/2 , 'bo')
    
    #Streamplot
    ax_stream.streamplot(z, x, np.transpose(electric_x), np.transpose(electric_z), color = 'black', arrowstyle='-', density = 1.5)
    ax_both.streamplot(z, x, np.transpose(electric_x), np.transpose(electric_z), color = 'black', arrowstyle='-', density = 1.5)
    
    #Quiver
    U = 3*electric_z / np.sqrt(electric_z**2 + electric_x**2);
    V = 3*electric_x / np.sqrt(electric_z**2 + electric_x**2);
    ax_quiver.quiver(Z[::3,::3], X[::3,::3], V[::3,::3], U[::3,::3], pivot='mid', color='r', units='inches')
    ax_both.quiver(Z[::2,::2], X[::2,::2], V[::2,::2], U[::2,::2], pivot='mid', color='r', units='inches')
    
    #Radiação
    x_charge, y_charge, z_charge = charge_position(t)
    Rz = Z - z_charge
    Rx = X - x_charge
    radiation = np.abs(Rz*electric_z + Rx*electric_x)
    im = ax_radiation.pcolor(radiation.transpose(), cmap=plt.cm.inferno, vmin=0, vmax=1)
    divider = make_axes_locatable(ax_radiation)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    figure_radiation.colorbar(im, cax=cax)
    
    ax_stream.set_aspect('equal')
    ax_quiver.set_aspect('equal')
    ax_both.set_aspect('equal')
    
    camera_stream.snap()
    camera_quiver.snap()
    camera_both.snap()
    camera_radiation.snap()

#%%RUN
def run(t_interval):
    for t in t_interval:
        electric_xyz = simulation(t)
        plot(electric_xyz, t)
        print('t = ' + str(t) + ' de ' + str(t_interval[-1]))

#%% CÁLCULOS
#Criar as figuras para os plots
figure_stream, ax_stream = plt.subplots()
figure_stream.set_size_inches((7,7))
camera_stream = Camera(figure_stream)

figure_quiver, ax_quiver = plt.subplots()
figure_quiver.set_size_inches((7,7))
camera_quiver = Camera(figure_quiver)

figure_both, ax_both = plt.subplots()
figure_both.set_size_inches((7,7))
camera_both = Camera(figure_both)

figure_radiation, ax_radiation = plt.subplots()
figure_radiation.set_size_inches((7,7))
camera_radiation = Camera(figure_radiation)

#Definição do intervalo de tempo que vamos considerar, rodando pra um período é suficiente
t_interval = np.arange(0, periodo(), time_step())

start = time.time()

#Rodando para os diferentes intervalos de tempo
run(t_interval)

#for t in t_interval:
#    electric_xyz = simulation(t)
#    plot(electric_xyz)
#    print('t = ' + str(t) + ' de ' + str(t_interval[-1]))
    
end = time.time()
tempo = end - start
print(tempo)
#%%
#Cria o gif
animation_stream = camera_stream.animate()
animation_quiver = camera_quiver.animate()
animation_both = camera_both.animate()
animation_radiation = camera_radiation.animate()

animation_stream.save('n' + str(n_malha()) + '_a' + str(amplitude()) + '_v' + str(v_max()) + '_stream.gif' , writer = 'imagemagick')
animation_quiver.save('n' + str(n_malha()) + '_a' + str(amplitude()) + '_v' + str(v_max()) + '_quiver.gif', writer = 'imagemagick')
animation_both.save('n' + str(n_malha()) + '_a' + str(amplitude()) + '_v' + str(v_max()) + '_both.gif' , writer = 'imagemagick')
animation_radiation.save('n' + str(n_malha()) + '_a' + str(amplitude()) + '_v' + str(v_max()) + '_radiation.gif' , writer = 'imagemagick')
