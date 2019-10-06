import numpy as np
import math
import matplotlib.pyplot as plt
import time

#cria malha de dimensao 201 por 201
def malha():
    m = np.zeros((201,201,2), dtype= 'int')
    for i in range(201):
        for j in range(201):
            m[i][j] = np.array([i-100,j-100])
    return m


#chute inicial para o tempo retardado
def chute0():
    return 0.0 

def tolerancia():
    return 0.001

#fazemos unidades com q=c=1.
#amplitude do movimento
def amplitude():
    return 10.0

#frequencia angular
def omega():
    return 0.025

#PS: para termos sempre v<c devemos ter omega*amplitude < c=1


#posicao da carga em funcao do tempo
def pos_carga(t):
    b = amplitude()
    w = omega()
    
    x=0.0
    z = b*np.cos(w*t)
    
    r = np.array([x,z])
    return r
    
#velocidade da carga
def vel_carga(t):
    b = amplitude()
    w = omega()
    
    vx = 0.0
    vz = -b*w*np.sin(w*t)
    
    v = np.array([vx,vz])
    return v

#calcula a funcao f(t_ret) da qual queremos as raizes
def f(r,t,t_ret):
    dist = np.sqrt( np.sum((r-pos_carga(t_ret))**2 ) )
    fun = t - t_ret - dist
    return fun
    
#calcula a derivada da funcao f(t_ret) da qual queremos as raizes
def df(r,t,t_ret):
    dist = np.sqrt( np.sum((r-pos_carga(t_ret))**2 ) )
    dfun = np.sum( (r-pos_carga(t_ret))*vel_carga(t_ret) )/dist - 1
    return dfun
    

#calcula o tempo retardado para um dado ponto (x,z,t)
def tempo_retardado(r,t, chute):
    tol = tolerancia()
    efe = f(r,t,chute)
    while efe > tol:
        chute = chute - f(r,t,chute)/df(r,t,chute)
        efe = f(r,t,chute)
    return chute    
    
#calcula o campo eletrico para um dado instante de tempo
def campo_eletrico(t):
    ponto = malha()
    campo = np.zeros((201,201,2), dtype= 'float')
    
    w = omega()
    
    chute = chute0() #chute inicial para t_ret
    for i in range(101,200):
        for j in range(0,200):
            t_ret = tempo_retardado(ponto[i][j], t, chute)
            
            r = ponto[i][j]
            dR = r - pos_carga(t_ret) #separacao no tempo retardado
            dist = np.sum(dR**2)     #distancia ao tempo retardado
            v = vel_carga(t_ret)     #velocidade no tempo retardado
            a = -w**2*(r-dR)         #aceleracao no tempo retardado
                        
            u = dR/dist - v
            
            #calcula o campo por aquela formula dado o tempo retardado. Eu uso alguma peculiaridades do nosso sistema (movimento unidimensional no eixo z "[1]" pra otimizar umas contas, abrindo umas componentes explicitamente) 
            campo[i][j] = ( (1-v[1]**2)*u*dist + np.array([dR[1],-dR[0]])*a[1]*r[0] )/np.sum(dR*u)**3
            
            #espelha o campo no eixo x
            campo[200-i][j][0] = -campo[i][j][0]
            campo[200-i][j][1] = campo[i][j][1]
    
    return campo 
    
t = 0.0    
t0 = time.time()
field = campo_eletrico(t)
tf = time.time()

print("tempo =", tf-t0) #tempo pra calcular em 1 instante

x = np.arange(-100,101)
z = np.arange(-100,101)


pi = math.pi
tempo = np.linspace(0.0,6*pi/omega(),90) #cria a sequencia temporal onde sao calculados os campos varrendo 3 periodos

k=0
for t in tempo:
    
    k=k+1
    print(k,t)
    
    field = campo_eletrico(1.0*t)

    electric_x = np.transpose(field[:,:,0])
    electric_z = np.transpose(field[:,:,1])

    plt.close(k-1)
    plt.figure(k)
    plt.streamplot(x,z , electric_x, electric_z, color='black', linewidth=1, cmap=plt.cm.inferno, 
                      density=1, arrowstyle='->', arrowsize=1.5)
    x_charge, z_charge = pos_carga(t)
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    plt.plot(x_charge, z_charge, 'bo')
    #plt.set_aspect('equal')
    plt.savefig(str(k)+".png")
