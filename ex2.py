#! /usr/bin/python
# -*- coding: utf-8 -*-

# Mairon de Souza Wolniewicz    (16250094)

# Universidade Federal de Santa Catarina - CTJ
# Engenharia Aeroespacial
# EMB5413 Mecânica dos Fluidos Computacional

# -----------------------------------------------
# EXERCÍCIO 2
# -----------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save

# -----------------------------------------------
# PROPRIEDADES

# propriedades do escoamento
h = 4.0             # [W/m^2 K] coeficiente convectivo
T_inf = 30.0+273    # [K]       temperatura do escoamento
u = 0.1             # [m/s]     velocidade de laminação

# propriedades da placa
T_0 = 600.0+273     # [K]       temperatura em x=0
rho = 2000.0        # [kg/m^3]  densidade do material
k = 400.0           # [W/mK]    condutividade térmica
cp = 1.0            # [J/kgK]   calor específico do material à pressão constante
L = 10.0            # [m]       comprimento avaliado da placa
w = 1.0             # [m]       largura da placa
t = 10e-3           # [m]       espessura da placa
P = 2*t + 2*w       # [m]       perímetro da seção transversal da placa
A = t*w             # [m]       área da seção transversal da placa

# -----------------------------------------------
# FUNÇÕES

# Perfil de temperatura de acordo com a equação analítica que governa o problema
# Input:    x       [m]   posição ao longo da placa
# Output:   T       [K]   temperatura na dada posição
def TemperatureExact(x):
    m = (h*P/(k*A))**0.5
    return T_inf + (T_0-T_inf) * np.sinh(m*(L-x))/np.sinh(m*L)

# Central Differencng Scheme
# Calcula o perfil de temperatura ao longo da placa dividida em N volumes finitos
# As propriedades nas paredes dos volumes são definidas pela aproximação de diferenças centrais
# Input:    N             número de volumes para discretização da placa
# Output:   T_num   [K]   vetor com as temperaturas em cada volume
#           x       [m]   vetor com a posição central de cada volume
def CDS(N):

    dx = L/float(N) # tamanho dos volumes

    C = h*P/cp * dx
    G = k*A/cp
    D = G/dx
    F = rho*u*A

    # listas de coeficientes:
    aw = np.zeros(N)
    ae = np.zeros(N)
    Su = np.zeros(N)
    Sp = np.zeros(N)
    ap = np.zeros(N)

    T = np.zeros(N) # vetor com propriedades

    # condição de contorno no primeiro volume:
    aw[0] = 0
    ae[0] = D - F*0.5
    Sp[0] = -(2*D + F + C)
    Su[0] = (2*D + F)*T_0 + C*T_inf
    ap[0] = ae[0] + aw[0] - Sp[0]

    # condição de contorno no N-ésimo volume:
    aw[-1] = D + F*0.5
    ae[-1] = 0
    Sp[-1] = -(- F + C)
    Su[-1] = (-F + C)*T_inf
    ap[-1] = ae[-1] + aw[-1] - Sp[-1]

    # condição para volumes intermediários:
    for i in range(1,N-1,1):
        aw[i] = D+F*0.5
        ae[i] = D-F*0.5
        Sp[i] = -C
        Su[i] = C*T_inf
        ap[i] = ae[i] + aw[i] - Sp[i]

    # matriz de coeficientes [Am]:
    Am = np.zeros([N,N])

    Am[0,0] = ap[0]
    Am[0,1] = -ae[0]

    Am[-1,-1] = ap[-1]
    Am[-1,-2] = -aw[-1]

    for i in range(1,N-1):
        Am[i,i] = ap[i]
        Am[i,i-1] = -aw[i]
        Am[i,i+1] = -ae[i]

    # vetor de propriedades resultante do procedimento numérico:
    T_num = np.linalg.solve(Am,Su)

    # discretização do intervalo numérico para o gráfico:
    X = np.zeros(N)
    X[0] = dx/2
    for i in range(1,N):
        X[i] = X[i-1]+dx

    return T_num, X

# Upwind Differencng Scheme
# Calcula o perfil de temperatura ao longo da placa dividida em N volumes finitos
# As propriedades nas paredes dos volumes são definidas considerando a direção do escoamento
# Input:    N             número de volumes para discretização da placa
# Output:   T_num   [K]   vetor com as temperaturas em cada volume
#           x       [m]   vetor com a posição central de cada volume
def UDS(N):

    dx = L/float(N) # tamanho dos volumes

    C = h*P/cp * dx
    G = k*A/cp
    D = G/dx
    F = rho*u*A

    # listas de coeficientes:
    aw = np.zeros(N)
    ae = np.zeros(N)
    Su = np.zeros(N)
    Sp = np.zeros(N)
    ap = np.zeros(N)

    T = np.zeros(N) # vetor com propriedades

    # condição de contorno no primeiro volume:
    aw[0] = 0
    ae[0] = D
    Sp[0] = -(2*D + F + C)
    Su[0] = (2*D + F)*T_0 + C*T_inf
    ap[0] = ae[0] + aw[0] - Sp[0]

    # condição de contorno no N-ésimo volume:
    aw[-1] = D + F
    ae[-1] = 0
    Sp[-1] = -C
    Su[-1] = C*T_inf
    ap[-1] = ae[-1] + aw[-1] - Sp[-1]

    # condição para volumes intermediários:
    for i in range(1,N-1,1):
        aw[i] = D+F
        ae[i] = D
        Sp[i] = -C
        Su[i] = C*T_inf
        ap[i] = ae[i] + aw[i] - Sp[i]

    # matriz de coeficientes [A]:
    Am = np.zeros([N,N])

    Am[0,0] = ap[0]
    Am[0,1] = -ae[0]

    Am[-1,-1] = ap[-1]
    Am[-1,-2] = -aw[-1]

    for i in range(1,N-1):
        Am[i,i] = ap[i]
        Am[i,i-1] = -aw[i]
        Am[i,i+1] = -ae[i]

    # vetor de propriedades resultante do procedimento numérico:
    T_num = np.linalg.solve(Am,Su)

    # discretização do intervalo numérico para o gráfico:
    X = np.zeros(N)
    X[0] = dx/2
    for i in range(1,N):
        X[i] = X[i-1]+dx

    return T_num, X

# -----------------------------------------------
# RESULTADOS

# 1)
u = 0                               # [m/s] velocidade de laminação da placa
x = np.linspace(0,10,100)           # intervalo analítico
plt.plot(x,TemperatureExact(x)-273,'k',linewidth=2,
         label='Analitico') # plot da solução analítica

T_num, X = CDS(50)
plt.plot(X,T_num-273,'--o',label='CDS 50 volumes')
T_num, X = CDS(5)
plt.plot(X,T_num-273,'--o',label='CDS 5 volumes')

T_num, X = UDS(50)
plt.plot(X,T_num-273,'--+',label='UDS 50 volumes')
T_num, X = UDS(5)
plt.plot(X,T_num-273,'--+',label='UDS 5 volumes')

plt.grid()
plt.xlabel("x [m]")
plt.ylabel("T [K]")
plt.legend()
# tikz_save("q1a.tex")
plt.show()

# 2)
u = 0.1 # [m/s] velocidade de laminação da placa

T_num, X = CDS(50)
plt.plot(X,T_num-273,'--o',label='CDS 50 volumes u = 0.1 m/s')
T_num, X = CDS(5)
plt.plot(X,T_num-273,'--o',label='CDS 5 volumes u = 0.1 m/s')

u = 1.0 # [m/s] velocidade de laminação da placa

T_num, X = CDS(50)
plt.plot(X,T_num-273,'--+',label='CDS 50 volumes u = 1.0 m/s')
T_num, X = CDS(5)
plt.plot(X,T_num-273,'--+',label='CDS 5 volumes u = 1.0 m/s')

plt.grid()
plt.xlabel("x [m]")
plt.ylabel("T [K]")
plt.legend()
# tikz_save("q2a.tex")
plt.show()

# 3)
u = 0.1 # [m/s] velocidade de laminação da placa

T_num, X = UDS(50)
plt.plot(X,T_num-273,'--o',label='UDS 50 volumes u = 0.1 m/s')
T_num, X = UDS(5)
plt.plot(X,T_num-273,'--o',label='UDS 5 volumes u = 0.1 m/s')

u = 1.0 # [m/s] velocidade de laminação da placa

T_num, X = UDS(50)
plt.plot(X,T_num-273,'--+',label='UDS 50 volumes u = 1.0 m/s')
T_num, X = UDS(5)
plt.plot(X,T_num-273,'--+',label='UDS 5 volumes u = 1.0 m/s')

plt.grid()
plt.xlabel("x [m]")
plt.ylabel("T [K]")
plt.legend()
# tikz_save("q3a.tex")
plt.show()
