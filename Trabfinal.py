#!/usr/bin/env python
# coding: utf-8

# In[1]:


from numpy import exp, abs, angle, conj
import numpy as np
from scipy.constants import mu_0, epsilon_0
from scipy.special import k1, k0, i1, i0 #fucoes que representam as funcoes de bessel
from mpmath import *
mp.dps = 25
mp.pretty = True


class Linha_transmissao:
#Classe que simplifica o calculo dos parametros de uma LT (matriz de impedancia, matriz de admitancia...)
#calculando os parametros com base nos inputs fornecidos (ver definicao de cada um dentro do __init__)
    def __init__(self, r_int, r_ext, nfase, npr, xc, yc, rhoc, rhoc_pr, rf, rpr):
        #r_int: raio interno do condutor de fase
        #r_ext: raio externo do condutor de fase
        #nfase: numero de condutores de fase
        #npr: numero de para-raios
        #xc: posicoes x dos centros de todos os condutores (fase e para-raio)
        #yc: posicoes x dos centros de todos os condutores (fase e para-raio)
        #rhoc: rho do condutor de fase
        #rhoc_pr: rho do para-raio
        #rf: raio externo do condutor de fase
        
        self.f = 60
        self.omega = 2*np.pi*self.f
        self.epsilon_r = 10 
        self.sigma_s = 1*10**(-3)
        self.r_int = r_int
        self.r_ext = r_ext
        self.nfase = nfase
        self.npr = npr
        self.xc = xc
        self.yc = yc
        self.rhoc = rhoc
        self.rhoc_pr = rhoc_pr
        self.rf = rf
        self.rpr = rpr
        self.ri = r_int + 10**(-6)
        #sqrt so funciona se for do modulo mpmath, se usar o do numpy nao funciona
        self.gama_S = sqrt(1j*self.omega*mu_0*(self.sigma_s + 1j*self.omega*self.epsilon_r*epsilon_0))
        self.gama_ar = 1j*self.omega*np.sqrt(mu_0*epsilon_0)
        #sqrt so funciona se for do modulo mpmath
        self.eta = sqrt(self.gama_S**2 - self.gama_ar**2)
        self.ncond = len(self.xc)


    # Calcula a matriz de impedancia de retorno do solo
    def S1(self):
        ncond = len(self.xc)
        s1 = np.eye(ncond)*1j
        for i in range(ncond):
            for j in range(ncond):
                if(i != j):
                    s1[i, j] = log(
                        1 + (2/(self.eta * sqrt((self.xc[i] - self.xc[j])**2 + (self.yc[i]+self.yc[j])**2))))
                elif (i+1) <= (ncond-self.npr):
                    s1[i, j] = log(
                        (1 + 2/(self.eta * sqrt(4*self.yc[i]**2 + self.rf**2))))
                else:
                    s1[i, j] = log(
                        (1 + 2/(self.eta * sqrt(4*self.yc[i]**2 + self.rpr**2))))
        return s1


    def Mpot(self):
        # Calcula a matriz dos potenciais
        ncond = len(self.xc)
        pot = np.eye(ncond)*1j
        for i in range(ncond):
            for j in range(ncond):
                if i != j:  # entre fases
                    num = (self.xc[i]-self.xc[j])**2+(self.yc[i]+self.yc[j])**2
                    den = (self.xc[i]-self.xc[j])**2+(self.yc[i]-self.yc[j])**2
                    pot[i, j] = 0.5*log(num/den)
                elif (i+1) <= (ncond-self.npr):  # condutor de fase
                    pot[i, j] = log((2*(self.yc[i])/self.rf))
                else:  # para-raio
                    pot[i, j] = log((2*(self.yc[i])/self.rpr))
        return pot
        

    def Zint(self):
        # Calcula a impedancia interna de um condutor cilindrico
        etac = np.sqrt((1j*self.omega*mu_0)/self.rhoc)
        zint = self.rhoc*(etac/(2*pi*self.rf)) *             (i0(abs(etac)*self.rf)/i1(abs(etac)*self.rf))
        return zint


    def Zinttub(self):
        # Calcula a impedancia interna de um condutor tubular
        etac = np.sqrt((1j*self.omega*mu_0)/self.rhoc)
        num = i0(abs(etac*self.rf))*k1(abs(etac*self.ri)) +             k0(abs(etac*self.rf))*i1(abs(etac*self.ri))
        den = i1(abs(etac*self.rf))*k1(abs(etac*self.ri)) -             i1(abs(etac*self.ri))*k1(abs(etac*self.rf))
        zin = self.rhoc*(etac/(2*pi*self.rf))*(num/den)
        return zin

    def Zin(self):
        ncond = len(self.xc)
        zin = np.eye(ncond)*1j
        for i in range(ncond):
            for j in range(ncond):
                if i != j:  # entre fases
                    zin[i, j] = 0
                elif (i+1) <= (ncond-self.npr):  # cabos de fase
                    zin[i, j] = self.Zinttub()
                else:  # pararaio
                    zin[i, j] = self.Zint()
        return zin

    def impedancia(self):
        #Calcula a matriz de impedancia
        Z = self.Zin() + (((1j*self.omega*mu_0)/2/pi) * (self.Mpot() + self.S1()))
        return Z*1000 #Saida da matriz de impedancia em ohms/km

    def admitancia(self):
        #Calcula a matriz de admitancia
        Y =  1j*self.omega*2*np.pi*epsilon_0*(np.linalg.inv(self.Mpot())) + 3.0*10**(-11)*np.eye(self.ncond)
        return Y*1000 #Saida da matriz de admitancia em ohms/km


# In[2]:


from numpy import sqrt

# função para calcular o raio equivalente
def raio_eq(n, rext, R):
    #R = raio do cicrulo que sera o novo condutor
    #rext = raio do feixe de condutor
    re = (rext*n*(R**(n-1)))**(1.0/n)
    return re

# potencia caracteristica

def Pnat(Vs, Z, Y):
    Zc = sqrt(Z/Y)
    Pnat = (Vs*Vs)/Zc.real
    return Pnat


# In[7]:


from numpy import sqrt
import numpy as np
from numpy.linalg import inv

np.set_printoptions(linewidth=500)

#------------------------------------------------------------------------------------------------------
# Dados dos condutores de fase
# Name: (r0,r1,pfase) ohms/m
CondutoresEspecs = {"Bluejay": (8.702*10**-3, 15.977*10**-3, 29.544*10**-9),
                    "Rail": (8.702*10**-3, 14.796*10**-3, 29.538*10**-9),
                    "3/8 EHS": (0, 4.570*10**-3, 276.470*10**-9)}

#------------------------------------------------------------------------------------------------------

#Name: (rext,rint,(Xc-1,Xc0,Xc1),(Yc-1,Yc,Yc+1),(Xpr,Ypr,Xpr2,Ypr2),n)
# n é o numero de condutores por fase
# Xc cordenadas dos centros dos 4 condutores
# Yc cordenadas dos centros dos 4 condutores PS: na formula subtrai por 2/3 da flecha
#Pros dois primeiros botei o raio como sendo o db*sqrt(2)/2
CondutoresPos = {"Bluejay": (raio_eq(4, CondutoresEspecs["Bluejay"][1], .475*sqrt(2)/2),raio_eq(4, CondutoresEspecs["Bluejay"][0], .475*sqrt(2)/2) ,(-15.85, 0, 15.85), (35.9-2*20.9/3, 35.9-2*20.9/3, 35.9-2*20.9/3), (-14.45, 45.9-2*14.7/3, 14.45, 45.9-2*14.7/3), 4),
                 "Rail Normal": (raio_eq(4, CondutoresEspecs["Rail"][1],.475*sqrt(2)/2), raio_eq(4, CondutoresEspecs["Rail"][0],.475*sqrt(2)/2), (-15, -11, -6, 6, 11, 15), (23.2, 33.2, 23.2, 23.2, 33.2, 23.2), (-8.8, 42.7, 8.8, 42.7), 4)
                 }
#------------------------------------------------------------------------------------------------------


# In[22]:


##distancia minima entre o eixo dos circuitos
# eixo seria o centro da torre
# referenciar norma tecnica ABNT NBR 5422 NBR5422 Projeto de linhas aéreas de
L=0
L=L+0.22+0.01*750000 #distancia minima entre os condutores localizados em circuitos de transmissao distintos


# In[39]:


#temos que reconstruir a matriz Xc. Yc por sua vez permanece igual
#Bluejay permanece igual, rail a gnt soma a distancia minima
delta1=abs(CondutoresPos["Rail Normal"][2][1]-CondutoresPos["Rail Normal"][2][0])
delta2=abs(CondutoresPos["Rail Normal"][2][1]-CondutoresPos["Rail Normal"][2][3])
delta3=abs(CondutoresPos["Rail Normal"][2][4]-CondutoresPos["Rail Normal"][2][3])
rail1=L+CondutoresPos["Bluejay"][2][2]
rail2=rail1+delta1
rail3=rail2+delta2
rail4=rail3+delta3
rail5=rail4+delta2
rail6=rail5+delta1

Xc=(CondutoresPos["Bluejay"][2][0],CondutoresPos["Bluejay"][2][1],CondutoresPos["Bluejay"][2][2],rail1,rail2,rail3,rail4,rail5,rail6)

Yc=[]
for i in (CondutoresPos["Bluejay"][3]):
    Yc.append(i)
for i in (CondutoresPos["Rail Normal"][3]):
    Yc.append(i)
    
    
print(Xc)
print(Yc)


# In[ ]:




