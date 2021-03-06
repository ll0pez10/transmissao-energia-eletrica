#!/usr/bin/env python
# coding: utf-8

from linha_transmissao import Linha_transmissao
from metodos_linhas import raio_eq, Pnat, derivacao, quadlinha, compenslinha
from numpy import savetxt

from numpy import sqrt
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=500)

#======================================================================================================
#==================================== Dados do Projeto ================================================
#======================================================================================================

# Especificacao dos condutores dos condutores
# Name: (r0,r1,pfase) (raio interno, raio externo, resistividade) ohms/m
CondutoresEspecs = {"Bluejay": (8.702*10**-3, 15.977*10**-3, 29.544*10**-9),
                    "Rail": (8.702*10**-3, 14.796*10**-3, 29.538*10**-9),
                    "3/8 EHS": (0, 4.570*10**-3, 276.470*10**-9)}

#------------------------------------------------------------------------------------------------------

#Posicao dos condutores
#Name: (rext,rint,(Xc-1,Xc0,Xc1),(Yc-1,Yc,Yc+1),(Xpr,Ypr,Xpr2,Ypr2),n)
# n é o numero de condutores por fase
# Xc cordenadas dos centros dos condutores equivalentes
# Yc cordenadas dos centros dos condutores equivalentes PS: na formula subtrai por 2/3 da flecha
#Pros dois primeiros botei o raio como sendo o db*sqrt(2)/2
CondutoresPos = {"Bluejay": (raio_eq(4, CondutoresEspecs["Bluejay"][1], .475*sqrt(2)/2),
                             raio_eq(4, CondutoresEspecs["Bluejay"][0], .475*sqrt(2)/2),
                             (-15.85, 0, 15.85),
                             (35.9-2*20.9/3, 35.9-2*20.9/3, 35.9-2*20.9/3),
                             (-14.45, 45.9-2*14.7/3, 14.45, 45.9-2*14.7/3), 4),
                 "Rail Normal": (raio_eq(4, CondutoresEspecs["Rail"][1],.475*sqrt(2)/2),
                                 raio_eq(4, CondutoresEspecs["Rail"][0],.475*sqrt(2)/2),
                                 (-15, -11, -6, 6, 11, 15),
                                 (23.2, 33.2, 23.2, 23.2, 33.2, 23.2),
                                 (-8.8, 42.7, 8.8, 42.7), 4)
                 }



#====================================================================================
#===================== Calculo da Distancia entre as torres =========================
#====================================================================================

##distancia minima entre o eixo dos circuitos
# o eixo seria a linha vertical que corta a torre ao meio
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


#Xc e Yc sao listas com as coordenadas dos centros dos condutores equivalentes das duas torres
#Xpr e Ypr sao lista com as coordenadas dos centros dos pararaios dos das duas torres


#plt.plot(XcondRail,YcondRail,'bx',Xc[3:],Yc[3:],'ro')
#plt.show()


#=================================================================================================
#============================= Matrizes de Impedancia e Admitancia ===============================
#=================================================================================================

#Construcao da matriz de sequencias
a = np.exp(1j * np.deg2rad(120))
a2 = a**2
A = np.array([[1, 1, 1], [1, a2, a], [1, a, a2]]) #matriz de sequencias

npr = 2 #numero de pararaios


#------------------------------------------ Bluejay -------------------------------------------
tipo = "Bluejay"
r_ext = CondutoresEspecs[tipo][0]
r_int = CondutoresEspecs[tipo][1]
nfase = 3
xc = np.array(Xc[0:3] + Xpr[0:2]) #vetor com a posicao X de todos os condutores do Bluejay(fase e pararaio)
yc = np.array(Yc[0:3] + Ypr[0:2]) #vetor com a posicao Y de todos os condutores (fase e pararaio)
rhoc = CondutoresEspecs[tipo][2]
rhoc_pr = CondutoresEspecs["3/8 EHS"][2]
rf = r_ext
rpr = CondutoresEspecs["3/8 EHS"][1]

#intanciacao do objeto que representa a linha bluejay
LinhaBluejay = Linha_transmissao(r_int, r_ext, nfase, npr, xc, yc, rhoc, rhoc_pr, r_ext, rpr)
mat_pot = LinhaBluejay.Mpot()
#Retiramos as info do pararaio das matrizes atraves da reducao de kron
pot_abc_bluejay = 0j + np.zeros((nfase,nfase))

pot_abc_bluejay = mat_pot[0:nfase,0:nfase] - mat_pot[0:nfase,nfase:] @ inv(mat_pot[nfase:,nfase:]) @ mat_pot[nfase:,0:nfase]
savetxt('mat_pot.csv', pot_abc_bluejay, delimiter=',')

Z_bluejay = LinhaBluejay.impedancia()
Y_bluejay = LinhaBluejay.admitancia()
Z_bluejay = Z_bluejay.astype(np.csingle)
Y_bluejay = Y_bluejay.astype(np.csingle)

#Retiramos as info do pararaio das matrizes atraves da reducao de kron
Zabc_bluejay = 0j + np.zeros((nfase,nfase))
Yabc_bluejay = 0j + np.zeros((nfase,nfase))

Zabc_bluejay = Z_bluejay[0:nfase,0:nfase] - Z_bluejay[0:nfase,nfase:] @ inv(Z_bluejay[nfase:,nfase:]) @ Z_bluejay[nfase:,0:nfase]
Yabc_bluejay = Y_bluejay[0:nfase,0:nfase] - Y_bluejay[0:nfase,nfase:] @ inv(Y_bluejay[nfase:,nfase:]) @ Y_bluejay[nfase:,0:nfase]

z012_bluejay = inv(A)@Zabc_bluejay@A
y012_bluejay = inv(A)@Yabc_bluejay@A

#------------------------------------- Rail normal ---------------------------------------

#Atribuimos os valores que serão utilizados para chamar a classe. Raios, numero de condutores, posição espacial dos condutores, rhos
tipo = "Rail Normal"
r_ext = CondutoresPos[tipo][0]
r_int = CondutoresPos[tipo][1]
nfase = 6
xc = np.array(Xc[3:] + Xpr[2:]) #vetor com a posicao X de todos os condutores do Rail Normal(fase e pararaio)
yc = np.array(Yc[3:] + Ypr[2:]) #vetor com a posicao Y de todos os condutores do Rail Normal (fase e pararaio)
rhoc = CondutoresEspecs["Rail"][2]
rhoc_pr = CondutoresEspecs["3/8 EHS"][2]
rf = r_ext
rpr = CondutoresEspecs["3/8 EHS"][1]

#intanciacao do objeto que representa a linha rail normal
LinhaRail = Linha_transmissao(r_int, r_ext, nfase, npr, xc, yc, rhoc, rhoc_pr, r_ext, rpr)

#Uma vez construida a linha podemos determinar seus parametros de impedancia e admitancia
Z_rail = LinhaRail.impedancia()
Y_rail = LinhaRail.admitancia()
Z_rail = Z_rail.astype(np.csingle) #Codigo necessario para o python ler o tipo de variavel 
Y_rail = Y_rail.astype(np.csingle)
#Redução de Kron para eliminar o para-raio
Zabc_rail = 0j + np.zeros((6,6))
Yabc_rail = 0j + np.zeros((6,6))
Zabc_rail = Z_rail[0:nfase,0:nfase] - Z_rail[0:nfase,nfase:] @ inv(Z_rail[nfase:,nfase:]) @ Z_rail[nfase:,0:nfase]
Yabc_rail = Y_rail[0:nfase,0:nfase] - Y_rail[0:nfase,nfase:] @ inv(Y_rail[nfase:,nfase:]) @ Y_rail[nfase:,0:nfase]

#A = [[A, 0] matriz de sequencia precisa ser 6x6 agr
#    [0, A]]
zero = np.zeros((3,3)) #matriz de zeros auxiliar
aux1 = np.concatenate((A,zero))
aux2 = np.concatenate((zero,A))
A = np.concatenate((aux1,aux2),axis=1) #matriz de sequencia para o caso 6x6


z012_rail = inv(A)@Zabc_rail@A
y012_rail = inv(A)@Yabc_rail@A


#=================================================================================================
#============================================= Quadripolos =======================================
#=================================================================================================
L = 750 #km
Sbase = 300e6 #VA
Vbaseblue = 750e3 #V
Vbaserail = 500e3 #V
Vbasegera = 250e3 #V
Zbaseblue = Sbase/Vbaseblue #ohms
Zbaserail = Sbase/Vbaserail #ohms
#Ybase = 1/Zbase
Xg = 2.4 #impedancia em pu do gerador
Qg = np.array( [[1,Xg],[0, 1]] ) #quadripolo do gerador
Xt = 3 #impedancia em pu do transformador
Qt = np.array( [[1,Xt],[0, 1]] ) #quadripolo do transformador

#======================================== Bluejay ==============================================

#------ Linha Bluejay (modelo pi) -------
Z1 = z012_bluejay[1][1]/Zbaseblue*L #impedancia considerando a sequencia positiva
Y1 = y012_bluejay[1][1]*Zbaseblue*L #admitancia considerando a sequencia negativa

#Componentes da matriz de quadripolo
A1 = 1 + Z1*Y1/2
A2 = Z1
A3 = Y1*(1 + Z1*Y1/4)
A4 = 1 + Z1*Y1/2

QL_b = np.array( [[A1, A2],[A3, A4]] ) #Quadripolo da linha bluejay


#-------------- Testes --------------------------
#fator de potencia unitario na carga

Vg = 750e3 #tensao de entrada na linha
Vg = Vg/Vbaseblue
#Ig = Sg/(np.sqrt(3) * Vg) #corrente no gerador
#VI_carga = inv(Qg @ QL_b) @ np.array([Vg, Ig])

#Caso 1: Vazio sem compensacao -> Corrente na carga e nula, a segunda coluna do quadripolo e desconsiderada
print("- Linha em vazio -BLUEJAY ")
Vr = Vg/A1
Ig = A3*Vr
print("Vr = %.2e < %.2f   V =   %.2f pu" % (abs(Vr),np.rad2deg(np.angle(Vr)),abs(Vr)))
print("Ig = %.2e < %.2f  A" % (abs(Ig),np.rad2deg(np.angle(Ig))))


#Caso 2: Vazio com compensacao indutiva em derivacao na fonte e na carga
Y = (1 - A1)/A2 #valor da admitancia de cada reator indutivo adicionado (p. 240 do Fuchs Vol. 1)

A1c = A1 + A2*Y
A2c = A2
A3c = A3 + A1*Y + A4*Y + A2*Y*Y
A4c = A4 + A2*Y

Qcomp = np.array( [[A1c, A2c],[A3c, A4c]] ) #Quadripolo equivalente considerando a compensacao

print("\n\n- Linha com compensacao -")
Vr = Vg/A1c
Ig = A3c*Vr
print("Vr = %.2e < %.2f   V =   %.2f pu" % (abs(Vr),np.rad2deg(np.angle(Vr)),abs(Vr/Vg)))
print("Ig = %.2e < %.2f  A" % (abs(Ig),np.rad2deg(np.angle(Ig))))






#======================================== Rail Normal ===========================================


#------ Linha Rail normal (modelo pi) -------
Z1 = z012_rail[1][1]/Zbaserail*L #impedancia considerando a sequencia positiva
Y1 = y012_rail[1][1]*Zbaserail*L #admitancia considerando a sequencia negativa

#montagem do quadripolo considerando apenas o modelo pi da linha
A1 = 1 + Z1*Y1/2
A2 = Z1
A3 = Y1*(1 + Z1*Y1/4)
A4 = 1 + Z1*Y1/2

QL_r = np.array( [[A1, A2],[A3, A4]] ) #quadripolo da linha rail normal

#Caso 1: Vazio sem compensacao -> Corrente na carga e nula, o segunda coluna do quadripolo e desconsiderada
Vg = 500e3 #tensao de entrada na linha
Vg = Vg/Vbaserail
print("\n\n- Linha em vazio - RAIL")
Vr = Vg/A1
Ig = A3*Vr
print("Vr = %.2e ang( %.2f )   V =   %.2f pu" % (abs(Vr),np.rad2deg(np.angle(Vr)),abs(Vr)))
print("Ig = %.2e ang( %.2f ) A" % (abs(Ig),np.rad2deg(np.angle(Ig))))

#======================================== Grafico da variacao da tensão com a distancia do RAIL ===========================================

#%matplotlib inline
#O objetivo aqui e ir avancando do ponto inicial da linha ate a outra ponta. Nesse trajeto, quando encontrarmos um ponto em que a tensao
#ultrapassar 1.05, marcaremos esse ponto para depois colocar uma subestacao. Continuamos as iteracoes, considerando que no ponto que foi
#marcado a tensao volta ao normal, ate acharmos o proximo ponto que a tensao ultrapassar o limite de 1.05 pu.
Vg = 500e3 #tensao de entrada na linha
Vg = Vg/Vbaserail
dim=0
for l in range(750):
    dist=l-dim
    Z1 = z012_rail[1][1]/Zbaserail*dist #impedancia considerando a sequencia positiva
    Y1 = y012_rail[1][1]*Zbaserail*dist #admitancia considerando a sequencia negativa

    A1 = 1 + Z1*Y1/2
    A2 = Z1
    A3 = Y1*(1 + Z1*Y1/4)
    A4 = 1 + Z1*Y1/2

    QL_r = np.array( [[A1, A2],[A3, A4]] ) #quadripolo da linha rail normal
       
    #Caso 1: Vazio sem compensacao -> Corrente na carga e nula, o segundda coluna do quadripolo e desconsiderada

    Vr = Vg/A1
    Ig = A3*Vr
    
    
    plt.plot(l, Vr, 'o', color='black');

plt.xlabel('Distancia (km)')
plt.ylabel('tensão (PU)')
plt.title('Linha sem compensação shunt')

comp=np.array( [[1, 0],[0.7, 1]] )



fig2=plt.figure()
#primeiro a gnt descobre a matriz de quadripolos desse pedaço todo da linha
l=244
Z1 = z012_rail[1][1]/Zbaserail*l #impedancia considerando a sequencia positiva
Y1 = y012_rail[1][1]*Zbaserail*l #admitancia considerando a sequencia negativa

A1 = 1 + Z1*Y1/2
A2 = Z1
A3 = Y1*(1 + Z1*Y1/4)
A4 = 1 + Z1*Y1/2

#compensação shunt, pag 240 fuchs    
k=1#quanto queremos compensar, 1 = tudo. Uentrada/Usaida se nao    
Y=(k-A1)/A2

for l in range(750):    
    if l < 244:
        Z1 = z012_rail[1][1]/Zbaserail*l #impedancia considerando a sequencia positiva
        Y1 = y012_rail[1][1]*Zbaserail*l #admitancia considerando a sequencia negativa
        
        A1 = 1 + Z1*Y1/2
        A2 = Z1
        A3 = Y1*(1 + Z1*Y1/4)
        A4 = 1 + Z1*Y1/2
      
        #nova matriz de quadripolos originaria da [shunt][quadripolos][shunt]
        A = A1 +A2*Y
        B = A2
        C = A3+A1*Y+A4*Y+A2*Y*Y
        D = A4 + A2*Y
        
        Vr = Vg*A
        plt.plot(l, abs(Vr), 'o', color='black');
        plt.xlabel('Distancia (km)')
        plt.ylabel('tensão (PU)')
        plt.title('Linha com compensação shunt de 100%')


#=====================================================================================================
#============================================= Linha com carga =======================================
#=====================================================================================================
Vg = 1 #pu
Ig = 1.07 #pu para garantirmos assim a potencia na fonte a principio de 3000 MW, provavelmente teremos que aumentar dps devido as perdas na linha
P0 = 3*(Vg*np.conj(Ig))
print("A potencia ativa total do gerador é : {:.2f} GW \n" .format(P0.real))
L1 = 244 #km ate a 1 subestação
L2 = 267 #km ate a 2 subestação
L3 = 239 #km ate a 3 subestação

#======================================== Duas Redes em paralelo ==============================================
x1=0.06 #fator de com pensacao da impedancia
x2=0.3 #fator de compensacao da admitancia

R = np.array( [ [0,1,0],[0,0,1],[1,0,0] ] ) #matriz de rotacao usada para fazer a transposicao

def quadparalelo(L,compz,compy):
#Calcula o quadripolo equivalente dos dois circuitos em paralelo (bluejay e rail)

    Z1 = z012_bluejay[1][1]/Zbaseblue*L
    Z1 = complex(Z1.real,Z1.imag*compz) #incluindo a compensacao de impedancia
    Y1 = y012_bluejay[1][1]*Zbaseblue*L
    Y1 = complex(Y1.real,Y1.imag*compy) #incluindo a compensacao de admitancia
    
    #quadripolo bluejay
    A1 = 1 + Z1*Y1/2
    B1 = Z1
    C1 = Y1*(1 + Z1*Y1/4)
    D1 = 1 + Z1*Y1/2
    
    Z1 = z012_rail[1][1]/Zbaserail*L
    Z1 = complex(Z1.real,Z1.imag*compz) #incluindo a compensacao de impedancia
    Y1 = y012_rail[1][1]*Zbaserail*L 
    Y1 = complex(Y1.real,Y1.imag*compy) #incluindo a compensacao de admitancia

    #quadripolo rail
    A2 = 1 + Z1*Y1/2
    B2 = Z1
    C2 = Y1*(1 + Z1*Y1/4)
    D2 = 1 + Z1*Y1/2
    
    #quadripolo resultado do paralelo entre o bluejay e o rail
    A = (A1*B2+B1*A2)/(B1+B2)
    B = B1*B2/(B1+B2)
    C = C1+C2+((A1-A2)*(D2-D1)/(B1+B2))
    D = (B1*D2+D1*B2)/(B1+B2)
    
    #Calculando o paralelo do circuito equivalente acima com o outro circuito que sobrou da linha dupla
    
    A = (A*B2+B*A2)/(B+B2)
    B = B*B2/(B+B2)
    C = C+C2+((A-A2)*(D2-D)/(B+B2))
    D = (B*D2+D*B2)/(B+B2)
    
    
    return np.array( [[A, B],[C, D]] )

#Quad = trecho de 244 km @ trecho de 267 km @ trecho de 239 km

Quad = quadparalelo(L1,x1,x2) @ quadparalelo(L2,x1,x2)@ quadparalelo(L3,x1,x2)
print("Matriz de quadripolos")
print(Quad)

#[Vr,Ir] = np.linalg.solve(Quad,[Vg,Ig])
Ir=1
Vr=1
#Vr=(Vg-Quad[0,1]*Ir)/Quad[0,0]
#Ig=Quad[1,0]*Vr+Quad[1,1]*Ir
[Vg,Ig]=np.linalg.solve(inv(Quad),[Vr,Ir])


P0= (3*(Vg*np.conj(Ig))).real
print("Vg = "+str(abs(Vg))+" L"+str(np.rad2deg(np.angle(Vg))))
print("Ig = "+str(abs(Ig))+" L"+str(np.rad2deg(np.angle(Ig))))
print("\nPotencia na carga: %.2f GW" % (P0))

print("\n- Linhas em paralelo -")
print("Vr = %.2f" % abs(Vr))
print("Ir = %.2f" % abs(Ir))
print("Vr * Ir = %.2f" % (Vr*Ir).real)
P1=(3*(Vr*np.conj(Ir))).real
print("\nPotencia na carga: %.2f GW" % (P1))
print("\nPerdas na linha: %.2f Porcentos" % (100*(-P1+P0)/P0))