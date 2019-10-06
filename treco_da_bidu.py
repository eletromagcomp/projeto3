#Coloquei separado pra nao dar erro, pq nao tinha certeza se alguem ja tinha feito alguma alteração

import numpy as np
import matplotlib.pyplot as plt
import time
from celluloid import Camera
from scipy.optimize import newton

#%% VARIÁVEIS

def n_malha():
    return 400

def amplitude():
    return 2

def v_max():
    return 0.8

def time_step():
    return 0.5

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


theta = np.linspace(0,2*np.pi, 20)


X,Z = np.meshgrid(x, z)

start = time.time()

figure1, ax1 = plt.subplots()
figure1.set_size_inches((7,7))
camera_hiro = Camera(figure1)

figure2, ax2 = plt.subplots()
figure2.set_size_inches((7,7))
camera_bidu = Camera(figure2)

for t in t_interval:
    for i in range(n_malha()):
        for j in range(n_malha()):
            xyz = np.array((x[i], y, z[j]))
            t_ret = newton(t_retarded, guess, fprime=t_retarded_prime, args=(xyz, t), tol=10**(-3))
            malha_t_ret[i,j] = t_ret
            guess = t_ret
            electric_x[i, j], electric_y[i, j], electric_z[i, j] = electric_field(t_ret, xyz)
    print('t =' + str(t) + ' terminado')
    
    
    ax1.streamplot(X, Z, electric_z, electric_x, color = 'lightgray', arrowstyle='-', density = 2)
    
    U = electric_z
    V = electric_x
    
    #NORMALIZANDO
    U = U / np.sqrt(U**2 + V**2);
    V = V / np.sqrt(U**2 + V**2);
    Q = ax2.quiver(X, Z, U, V, pivot='mid', color='r', units='inches')
    '''
    Q = ax1.quiver(X, Z, U, V, pivot='mid', color='r', units='inches')
    '''
    x_charge, y_charge, z_charge = charge_position(t)
    ax1.plot(z_charge, x_charge, 'bo')
    
    ax2.plot(z_charge, x_charge, 'bo')
    
    ax1.set_aspect('equal')
    
    camera_hiro.snap()
    
    camera_bidu.snap()
    
end = time.time()
tempo = end - start
print(tempo)

start = time.time()

animation_hiro = camera_hiro.animate()
animation_hiro.save('linhas_de_campo.gif', writer = 'imagemagick')

end = time.time()
tempo = end - start
print(tempo)

start = time.time()

animation_bidu = camera_bidu.animate()
animation_bidu.save('campo_normalizado.gif', writer = 'imagemagick')

end = time.time()
tempo = end - start
print(tempo)
