from numpy import exp, abs, angle, conj, sqrt
import numpy as np
from scipy.constants import mu_0, epsilon_0
from scipy.special import k1, k0, i1, i0
from mpmath import *
mp.dps = 25
mp.pretty = True


class Linha_transmissao:
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
        self.omega = 2*pi*self.f
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
        self.ri = rint + 10**(-6)
        self.gama_S = np.sqrt(
            1j*omega*mu_0*(sigma_s + 1j*omega*epsilon_r*epsilon_0))
        self.gama_ar = 1j*omega*np.sqrt(mu_0*epsilon_0)
        self.eta = np.sqrt(gama_S**2 - gama_ar**2)

    # Calcula a matriz de impedancia de retorno do solo

    def S1(self, xc, npr, yc, rf, rpr):
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

    def Mpot(self, xc, yc, npr, rf, rpr):
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

    def zexternal(self, rho, npr, xc, yc, rf, rpr):
        p = sqrt(self.rho/(1j*self.omega*mu_0))
        ncond = len(self.xc)
        zout = np.eye(ncond)*1j
        for i in range(ncond):
            for j in range(ncond):
                if i != j:  # entre fases
                    num = (self.xc[i]-self.xc[j]**2) + \
                        (2*p+self.yc[i]+self.yc[j]**2)
                    den = (self.xc[i]-self.xc[j]**2)+(self.yc[i]-self.yc[j]**2)
                    zout[i, j] = 1j*self.omega*mu_0/(4*pi)*log(num/den)
                elif (i+1) <= (ncond-self.npr):
                    zout[i, j] = 1j*self.omega*mu_0 / \
                        (2*pi)*log((2*(self.yc[i]+p))/self.rf)
                else:
                    zout[i, j] = 1j*self.omega*mu_0 / \
                        (2*pi)*log((2*(self.yc[i]+p))/self.rpr)
        return zout

    def Zint(self, omega, rhoc, rf):
        # Calcula a impedancia interna de um condutor cilindrico
        etac = sqrt((1j*self.omega*mu_0)/self.rhoc)
        zint = rhoc*(etac/(2*pi*self.rf)) * \
            (i0(etac*self.rf)/i1(etac*self.rf))
        return zint

    def Zinttub(self, omega, rhoc, rf, ri):
        # Calcula a impedancia interna de um condutor tubular
        etac = sqrt((1j*omega*mu_0)/self.rhoc)
        num = i0(etac*self.rf)*k1(etac*self.ri) + \
            k0(etac*self.rf)*i1(etac*self.ri)
        den = i1(etac*self.rf)*k1(etac*self.ri) - \
            i1(etac*self.ri)*k1(etac*self.rf)
        zin = self.rhoc*(etac/(2*pi*self.rf))*(num/den)
        return zint

    def Zin(self, xc, npr, rhoc, rhoc_pr, rf, ri):
        ncond = len(self.xc)
        zin = np.eye(ncond)*1j
        rpr = 0  # qual a variavel?

        for i in range(ncond):
            for j in range(ncond):
                if i != j:  # entre fases
                    zin[i, j] = 0
                elif (i+1) <= (ncond-self.npr):  # cabos de fase
                    zin[i, j] = self.Zinttub(
                        self, omega, self.rhoc, self.rf, self.ri)
                else:  # pararaio
                    zin[i, j] = self.Zint(self.omega, self.rhoc_pr, self.rpr)
        return zin

    def impedancia(self, epsilon_r, sigma_s, r_int, r_ext, nfase, npr, xc, yc, rhoc, rhoc_pr, rf, rpr):
        Z = self.Zin(self.xc, self.npr, self.rhoc, self.rhoc_pr, self.rf, self.ri) + (((1j*omega*mu_0)/2/pi)
                                                                                      * (self.Mpot(self.xc, self.yc, self.npr, self.rf, self.rpr) + self.S1(self.xc, self.npr, self.yc, self.rf, self.rpr)))
        return Z

    def impedanciaY(self, epsilon_r, sigma_s, r_int, r_ext, nfase, npr, xc, yc, rhoc, rhoc_pr, rf, rpr):
        Y =  1j*omega*2*math.pi*epsilon_0*(np.linalg.inv(Mpot(xc,yc,npr,rf,rpr))) + 3.0*10**(-11)*np.eye(ncond)
        return Y

