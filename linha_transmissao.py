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
        zint = self.rhoc*(etac/(2*pi*self.rf)) * \
            (i0(abs(etac)*self.rf)/i1(abs(etac)*self.rf))
        return zint


    def Zinttub(self):
        # Calcula a impedancia interna de um condutor tubular
        etac = np.sqrt((1j*self.omega*mu_0)/self.rhoc)
        num = i0(abs(etac*self.rf))*k1(abs(etac*self.ri)) + \
            k0(abs(etac*self.rf))*i1(abs(etac*self.ri))
        den = i1(abs(etac*self.rf))*k1(abs(etac*self.ri)) - \
            i1(abs(etac*self.ri))*k1(abs(etac*self.rf))
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
