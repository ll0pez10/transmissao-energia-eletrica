#!/usr/bin/env python
# coding: utf-8


from linha_transmissao import Linha_transmissao
from metodos_linhas import raio_eq, Pnat

from numpy import sqrt
import numpy as np
from numpy.linalg import inv

np.set_printoptions(linewidth=500)

#------------------------------------------------------------------------------------------------------
# Dados dos condutores
# Name: (r0,r1,pfase) (raio interno, raio externo, resistividade) ohms/m
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


##distancia minima entre o eixo dos circuitos
# eixo seria o centro da torre
# referenciar norma tecnica ABNT NBR 5422 NBR5422 Projeto de linhas aéreas de
L=0
L=L+0.22+0.01*750000 #distancia minima entre os condutores localizados em circuitos de transmissao distintos


#temos que reconstruir a matriz Xc. Yc por sua vez permanece igual
#Bluejay permanece igual, rail a gnt soma a distancia minima
delta1=abs(CondutoresPos["Rail Normal"][2][1]-CondutoresPos["Rail Normal"][2][0])
delta2=abs(CondutoresPos["Rail Normal"][2][1]-CondutoresPos["Rail Normal"][2][2])
delta3=abs(CondutoresPos["Rail Normal"][2][2]-CondutoresPos["Rail Normal"][2][3])
rail1=L+CondutoresPos["Bluejay"][2][2]
rail2=rail1+delta1
rail3=rail2+delta2
rail4=rail3+delta3
rail5=rail4+delta2
rail6=rail5+delta1

#Condutores
Xc=[]
for i in (CondutoresPos["Bluejay"][2]):
    Xc.append(i) 

Xc=Xc+[rail1,rail2,rail3,rail4,rail5,rail6]

Yc=[]
for i in (CondutoresPos["Bluejay"][3]):
    Yc.append(i)
for i in (CondutoresPos["Rail Normal"][3]):
    Yc.append(i)
    
    
print(Xc)
print(Yc)

#Para-raio
Xpr=[]
Ypr=[]
for i in range(0,len(CondutoresPos["Bluejay"][4]),2):
    Xpr.append(CondutoresPos["Bluejay"][4][i]) 
    
for i in range(1,len(CondutoresPos["Bluejay"][4]),2):
    Ypr.append(CondutoresPos["Bluejay"][4][i]) 

delta=abs(CondutoresPos["Rail Normal"][2][0]-CondutoresPos["Rail Normal"][4][0]) #diferenca da distancia do pararaio pro primeiro condutor do rail
Xpr.append(rail1+delta)
Xpr.append(rail6-delta)

for i in range(1,len(CondutoresPos["Rail Normal"][4]),2):
    Ypr.append(CondutoresPos["Rail Normal"][4][i]) 

print(Xpr)
print(Ypr)
