#! /usr/bin/python
# -*- coding: utf-8 -*-

# Mairon de Souza Wolniewicz (16250094)
# Universidade Federal de Santa Catarina - CTJ
# Engenharia Aeroespacial
# -----------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

# ===============================================
# PARÂMETROS DO PROBLEMA

q = 0.3e6       # W/m^3
k = 25.0        # W/mK
h = 500.0       # W/m^2 K
T_inf = 92.0    # °C
L = 0.1         # m
T_B = q*L/h + T_inf # K

# ===============================================
# ANALÍTICO

def T(x):
    return q/(2*k) * (L**2 - x**2) + q*L/h + T_inf

# ===============================================
# NUMÉRICO

# função que realiza a resolução do problema para N volumes
# retorna um vetor T_num com as temperaturas calculadas
# e um vetor X com o intervalo discretizado

def Solucao(N):

    #N = 5           # número de volumes
    dx = L/float(N) # tamanho dos volumes

    # listas de coeficientes:
    aw = np.zeros(N)
    ae = np.zeros(N)
    Su = np.zeros(N)
    Sp = np.zeros(N)
    ap = np.zeros(N)

    T = np.zeros(N) # vetor com temperaturas

    # condição de contorno no primeiro volume:
    aw[0] = 0     
    ae[0] = k/dx
    Sp[0] = 0
    Su[0] = q*dx
    ap[0] = ae[0] + aw[0] - Sp[0]

    # condição de contorno no N-ésimo volume:
    aw[-1] = k/dx
    ae[-1] = 0
    Sp[-1] = -2.0*k/dx
    Su[-1] = 2*k/dx * T_B + q*dx
    ap[-1] = ae[-1] + aw[-1] - Sp[-1]

    # condição para volumes intermediários:
    for i in range(1,N-1,1):
        aw[i] = k/dx
        ae[i] = k/dx
        Sp[i] = 0
        Su[i] = q*dx
        ap[i] = ae[i] + aw[i] - Sp[i]
            
    # matriz de coeficientes [A]:
    A = np.zeros([N,N])

    A[0,0] = ap[0]
    A[0,1] = -ae[0]

    A[-1,-1] = ap[-1]
    A[-1,-2] = -aw[-1]

    for i in range(1,N-1):
        A[i,i] = ap[i]
        A[i,i-1] = -aw[i]
        A[i,i+1] = -ae[i]
        
    # vetor de temperaturas resultante do procedimento numérico:
    T_num = np.linalg.solve(A,Su)

    # discretização do intervalo numérico para o gráfico:
    X = np.zeros(N)
    X[0] = dx/2
    for i in range(1,N):
        X[i] = X[i-1]+dx
        
    return T_num, X

# ===============================================
# GRÁFICOS

# gráfico de temperaturas
plt.xlabel("x (m)")
plt.ylabel("T (Celcius)")
plt.grid()

x = np.linspace(0,L,100) # intervalo para plot analítico
plt.plot(x,T(x),'k',label="Analitico" )     # plot analítico
T_num, X = Solucao(5)
plt.plot(X,T_num,'or',label="5 volumes")
T_num, X = Solucao(20)
plt.plot(X,T_num,'+b',label="20 volumes")

plt.legend()
tikz_save("q1.tex")

plt.show()


# gráfico de erros
plt.xlabel("x (m)")
plt.ylabel("DT (Celcius)")
plt.grid()

T_num, X = Solucao(5)

# cálculo do erro relativo
DT = np.zeros(len(X))
for i in range(len(X)):
    DT[i] = abs(T(X[i])-T_num[i])
    print "%6.4f\t%6.4f\t%6.4f\t" %(T(X[i]),T_num[i],DT[i])

plt.plot(X,DT,'o-r',label="5 volumes")

print 

T_num, X = Solucao(20)

# cálculo do erro relativo
DT = np.zeros(len(X))
for i in range(len(X)):
    DT[i] = abs(T(X[i])-T_num[i])
    print "%6.4f\t%6.4f\t%6.4f\t" %(T(X[i]),T_num[i],DT[i])
    
plt.plot(X,DT,'+-b',label="20 volumes")

plt.legend()
tikz_save("q2.tex")

plt.show()

