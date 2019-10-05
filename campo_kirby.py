import numpy as np
import matplotlib as plt

def malha():
    m = np.zeros((201,201,2), dtype= 'int')
    for i in range(201):
        for j in range(201):
            m[i][j] = np.array([i-100,j-100])
    return m

tempo = np.arange(0,100)

#chute inicial para o tempo retardado
def chute0():
    return 0.0 

def tolerancia():
    return: tolerance = 0.001

#fazemos unidades com q=c=1.
#amplitude do movimento
def amplitude():
    return 10.0

#frequencia angular
def omega():
    return 0.01

#PS: para termos sempre v<c devemos ter omega*amplitude < c=1


#posicao da carga em funcao do tempo
def pos_carga(t):
    a = amplitude()
    w = omega()
    
    x=0.0
    z = a*np.cos(w*t)
    
    r = np.array([x,z])
    return r
    
#velocidade da carga
def vel_carga(t):
    a = amplitude()
    w = omega()
    
    vx = 0.0
    vz = -a*w*np.sin(w*t)
    
    v = np.array([vx,vz])
    return v

#calcula a funcao f(t_ret) da qual queremos as raizes
def f(r,t,t_ret):
    dist = np.sqrt( np.sum((r-pos_carga(t_ret))**2 ) )
    fun = t - t_r - dist
    return fun
    
#calcula a derivada da funcao f(t_ret) da qual queremos as raizes
def df(r,t,t_ret):
    dist = np.sqrt( np.sum**2(r-pos_carga(t_ret) ) )
    dfun = (r-pos_carga(t_ret))*vel_carga(t_ret)/dist - 1
    return dfun
    

#calcula o tempo retardado para um dado ponto (x,z,t)
def tempo_retardado(r,t, chute, tolerance):
    while f(r,t,chute) > tolerance:
        chute = chute - f(r,t,chute)/df(r,t,chute)
    return chute    
    
#calcula o campo eletrico para um dado instante de tempo
def campo_eletrico(r,t):
    ponto = malha()
    campo = np.zeros((201,201,2), dtype= 'float')
    for i in range(101,200):
        for j in range(0,200):
            #calcula
    
    return campo













