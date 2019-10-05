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

#%% VARIÁVEIS

def n_malha():
    return 400

def amplitude():
    return 5

def v_max():
    return 0.5

def time_step():
    return 1

def n_ciclos():
    return 1  
    
#%% TEMPO RETARDADO
def t_retarded(t_ret,*args):
    x,y,z = xyz
    return ((x - 0.1)**2 + y**2 + (z - amplitude()*np.cos(v_max()*t_ret/amplitude()))**2)**(1/2) - (t - t_ret)

#%% DERIVADA TEMPO RETARDADO
def t_retarded_prime(t_ret, *args):
    x,y,z = xyz
    return (v_max()*np.sin(t_ret*v_max()/amplitude())*(z - amplitude()*np.cos(t_ret*v_max()/amplitude())))*((z - amplitude()*np.cos(t_ret*v_max()/amplitude()))**2 + (x- 0.1)**2 + y**2)**(-1/2) + 1

#%%
def charge_position(t_ret):
    x = 0.1
    y = 0
    z = amplitude()*np.cos(t_ret*v_max()/amplitude())
    xyz_charge = np.array((x, y, z))
    return xyz_charge

#%%
def charge_velocity(t_ret):
    v_x = 0
    v_y = 0
    v_z = -v_max() * np.sin(t_ret*v_max()/amplitude())
    v_charge = np.array((v_x, v_y, v_z))
    return v_charge
#%%
def charge_acceleration(t_ret):
    a_x = 0
    a_y = 0
    a_z = -(v_max()**2)/amplitude() * np.cos(t_ret*v_max()/amplitude())
    a_charge = np.array((a_x, a_y, a_z))
    return a_charge

#%%
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
    
#%%
'''
x = np.arange(-n_malha()/2, n_malha()/2)
y = 0
z = np.arange(-n_malha()/2, n_malha()/2)
malha_t_ret = np.zeros((n_malha(), n_malha()))
periodo = 2*np.pi/(v_max()/amplitude())
t_interval = np.arange(0, n_ciclos()*periodo, time_step())
guess = 1
electric_x = np.zeros((n_malha(), n_malha()))
electric_y = np.zeros((n_malha(), n_malha()))
electric_z = np.zeros((n_malha(), n_malha()))


start = time.time()

figure, ax = plt.subplots()
figure.set_size_inches((7,7))
camera = Camera(figure)
for t in t_interval:
    for i in range(n_malha()):
        for j in range(n_malha()):
            xyz = np.array((x[i], y, z[j]))
            t_ret = newton(t_retarded, guess, fprime=t_retarded_prime, args=(xyz, t))
            malha_t_ret[i,j] = t_ret
            guess = t_ret
            electric_x[i, j], electric_y[i, j], electric_z[i, j] = electric_field(t_ret, xyz)
    print('t =' + str(t) + ' terminado')
    seed_points_x = np.linspace(-10, 10, 20)
    seed_points_z = np.sqrt(10**2 - seed_points_x**2)
    seed_points = np.array([seed_points_x, seed_points_z])
    ax.streamplot(x, z, electric_x, electric_z, color='black', linewidth=1, cmap=plt.cm.inferno, 
                      density=1, arrowstyle='->', start_points = seed_points.T, arrowsize=1.5)
    x_charge, y_charge, z_charge = charge_position(t)
    ax.plot(x_charge, z_charge, 'bo')
    ax.set_aspect('equal')
    camera.snap()
end = time.time()
tempo = end - start
print(tempo)
animation = camera.animate()
animation.save('celluloid_minimal_2.gif', writer = 'imagemagick')
'''
#%%
x = np.arange(-n_malha()/2, n_malha()/2)
y = 0
z = np.arange(-n_malha()/2, n_malha()/2)
malha_t_ret = np.zeros((n_malha(), n_malha()))
periodo = 2*np.pi/(v_max()/amplitude())
t = 50
guess = 1
electric_x = np.zeros((n_malha(), n_malha()))
electric_y = np.zeros((n_malha(), n_malha()))
electric_z = np.zeros((n_malha(), n_malha()))


start = time.time()

for i in range(n_malha()):
    print(i)
    for j in range(n_malha()):
        xyz = np.array((x[i], y, z[j]))
        t_ret = newton(t_retarded, guess, fprime=t_retarded_prime, args=(xyz, t))
        malha_t_ret[i,j] = t_ret
        guess = t_ret
        electric_x[i, j], electric_y[i, j], electric_z[i, j] = electric_field(t_ret, xyz)
#%%
        
figure, ax = plt.subplots()
figure.set_size_inches((7,7))
camera = Camera(figure)
ax.quiver(x, z, electric_x, electric_z)
x_charge, y_charge, z_charge = charge_position(t)
ax.plot(x_charge, z_charge, 'bo')
ax.set_aspect('equal')
camera.snap()
end = time.time()
tempo = end - start
print(tempo)
animation = camera.animate()
animation.save('celluloid_minimal_3.gif', writer = 'imagemagick')
