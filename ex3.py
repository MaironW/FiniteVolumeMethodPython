#! /usr/bin/python
# -*- coding: utf-8 -*-

# Mairon de Souza Wolniewicz    (16250094)

# Universidade Federal de Santa Catarina - CTJ
# Engenharia Aeroespacial
# EMB5413 Mecânica dos Fluidos Computacional

# -----------------------------------------------
# EXERCÍCIO 3
# -----------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

# -----------------------------------------------
# PROPRIEDADES

# propriedades do escoamento
h = 500.0           # [W/m^2 K] coeficiente convectivo
T_inf = 100.0       # [deg C]   temperatura do escoamento

# propriedades da parede
T_i = 0.0           # [deg C]   temperatura inicial
rho = 2500.0        # [kg/m^3]  densidade do material
k = 25.0            # [W/mK]    condutividade térmica
cp = 1.0            # [J/kgK]   calor específico do material à pressão constante
L = 0.1             # [m]       espessura da parede

# -----------------------------------------------
# FUNÇÕES

# Perfil de temperatura de acordo com a equação analítica que governa o problema
# Input:    x       [m]   posição ao longo da placa
#           t       [s]   tempo analisado
# Output:   T       [K]   temperatura nos dados posição e instante
def TemperatureExact(x,t):
    return 1.1795*np.exp(-1.0769**2*k*t/(rho*cp*L**2))*np.cos(1.0769*x/L)*(T_i-T_inf)+T_inf


# Calcula  o perfil de temperatura ao longo da parede dividida em N volumes finitos em um passo de tempo discreto dt.
# Input:    N           número de volumes para discretização da parede
#           dt      [s] passo de tempo da solução numérica
#           t       [s] instante de  interesse
# Output:   T       [K] vetor com as temperaturas em cada volume
#           x       [m] vetor com a posição central de cada volume
def TransientSolution(N,dt,t):
    dx = L/float(N) # tamanho dos volumes

    # listas de coeficientes
    aw = np.zeros(N)
    ae = np.zeros(N)
    Su = np.zeros(N)
    Sp = np.zeros(N)

    T = np.ones(N)*T_i # vetor com temperaturas

    ap = rho*cp*dx/dt

    # condição de contorno no primeiro volume:
    aw[0] = 0
    ae[0] = k/dx

    # condição de contorno no último volume:
    T_B = T_inf + (T[-1]-T_inf)/(h*dx/(2*k)+1)  # temperatura na parede leste
    aw[-1] = k/dx
    ae[-1] = 0
    Su[-1] = 2*k/dx * T_B
    Sp[-1] = -2*k/dx

    # condição de contorno para volumes intermediários
    for i in range(1,N-1):
        aw[i] = k/dx
        ae[i] = k/dx

    # solução transiente
    ts = 0
    To = np.copy(T)

    while ts<t:
        T[0] = (ae[0]*To[1] + (ap-ae[0]-aw[0])*To[0])/ap

        for i in range(1,N-1):
            T[i] = (aw[i]*To[i-1] + ae[i]*To[i+1] + (ap-ae[i]-aw[i])*To[i])/ap

        T_B = T_inf + (To[-1]-T_inf)/(h*dx/(2*k)+1)  # temperatura na parede leste
        Su[-1] = 2*k/dx * T_B
        T[-1] = (aw[-1]*To[-2] + (ap - aw[-1] + Sp[-1])*To[-1] + Su[-1])/ap

        print rho*cp*dx**2/(3*k),(- aw[-1] + Sp[-1]),(ap - aw[-1] + Sp[-1])
        ts+=dt
        To = np.copy(T)

    # discretização do intervalo numérico para o gráfico:
    X = np.zeros(N)
    X[0] = dx/2
    for i in range(1,N):
        X[i] = X[i-1]+dx

    return T, X

# -----------------------------------------------
# RESULTADOS

x = np.linspace(0,L,100)    # [m]   intervalo analítico

plt.figure()

# 1)
# plt.subplot(1,2,1)
plt.plot(x,TemperatureExact(x,0.2),'k',label='Analitico')

T_num, X = TransientSolution(10,0.002,0.2)
plt.plot(X,T_num,'bo-',label='$\Delta t=0.002$ s')
T_num, X = TransientSolution(10,0.005,0.2)
plt.plot(X,T_num,'ro-',label='$\Delta t=0.005$ s')
T_num, X = TransientSolution(10,0.0053,0.2)
plt.plot(X,T_num,'go-',label='$\Delta t=0.0053$ s')

plt.xlabel("$x$ (m)")
plt.ylabel("$T$ (\degree C)")
plt.grid()
plt.legend()
tikz_save("q1a.tex")
plt.show()

# 2)
# plt.subplot(1,2,2)
plt.plot(x,TemperatureExact(x,0.2),'b',label='Analitico')
plt.plot(x,TemperatureExact(x,1.0),'r',label='Analitico')

T_num, X = TransientSolution(10,0.002,0.2)
plt.plot(X,T_num,'bo',label='$t=0.2$ s')
T_num, X = TransientSolution(10,0.002,1.0)
plt.plot(X,T_num,'ro',label='$t=1.0$ s')

plt.xlabel("$x$ (m)")
plt.ylabel("$T$ (\degree C)")
plt.grid()
plt.legend()
tikz_save("q2a.tex")
plt.show()





