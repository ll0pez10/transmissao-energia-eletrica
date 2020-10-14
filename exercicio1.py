from linha_transmissao import Linha_transmissao
from metodos_linhas import raio_eq, Pnat
from numpy import sqrt
import numpy as np
from numpy.linalg import inv

np.set_printoptions(linewidth=500)

# Dados dos condutores de fase
# Name: (r0,r1,pfase) ohms/m
CondutoresEspecs = {"Bluejay": (8.702*10**-3, 15.977*10**-3, 29.544*10**-9),
                    "Rail": (8.702*10**-3, 14.796*10**-3, 29.538*10**-9),
                    "Ruddy": (7.821*10**-3, 14.364*10**-3, 29.559*10**-9),
                    "Grossbeak": (7.456*10**-3, 12.573*10**-3, 29.538*10**-9),
                    "Dove": (6.988*10**-3, 11.773*10**-3, 29.567*10**-9),
                    "Penguim": (4.123*10**-3, 7.150*10**-3, 28.554*10**-9),
                    "Leghorn": (4.840*10**-3, 6.7180*10**-3, 29.586*10**-9),
                    "Minorca": (4.390*10**-3, 6.0960*10**-3, 29.579*10**-9),
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

# ============ 500kV rail convencional ===============
#2 ultimos elementos sao pararraio, da forma (x, y)
n_rail_comp = 3
R_rail_compExt = (CondutoresEspecs["Rail"][1] *
               np.sqrt(3))/3 + CondutoresEspecs["Rail"][1]

R_rail_compInt = (CondutoresEspecs["Rail"][0] *
               np.sqrt(3))/3 + CondutoresEspecs["Rail"][0]

xy = ((-11.000, 10.736), (-10.771, 11.132),
      (-11.229, 11.132), (0.000, 10.736), (0.229,
                                           11.132), (-0.229, 11.132), (11.000, 10.736),
      (10.771, 11.132), (11.229, 11.132), (-6.000, 22.000), (6.000, 22.000))

#centro dos condutores equivalentes (os dois ultimos pares sao dos para-raios)
xyc = (((-10.771 + -11.229)/2, (2*11.132 + 10.736)/3),
       ((0.229 + -0.229)/2, (2*11.132 + 10.736)/3),
       ((10.771 + 11.229)/2, (2*11.132 + 10.736)/3),
       (-6.000, 22.000), (6.000, 22.000))

R = sqrt((xyc[0][0] - xy[0][0])**2 + (xyc[0][1] - xy[0][1])**2)

CondutoresPos["Rail Convencional"] = (raio_eq(3, CondutoresEspecs["Rail"][1], R),raio_eq(4, CondutoresEspecs["Rail"][0], R),
                                    (xyc[0][0], xyc[1][0], xyc[2][0]),
                                      (xyc[0][1], xyc[1][1], xyc[2][1]),
                                      (xyc[3][0], xyc[3][1], xyc[4][0], xyc[4][1]), 3)


# ================= 500kV rail compacto ====================
#2 ultimos elementos sao pararraio, da forma (x, y)

n_rail_comp = 4
R_rail_comp = 0

xy2 = ((-4.271, 10.771), (-4.271, 11.229), (-4.729, 11.229), (-4.729, 10.771), (0.229, 15.271),
       (0.229, 15.729), (-0.229, 15.729), (-0.229, 15.271), (4.271,10.771),
       (4.271, 11.229), (4.729, 11.229), (4.729, 10.771), (-3.500, 26.000), (3.500, 26.000))
       
# centro dos condutores equivalentes (os dois ultimos pares sao dos para-raios)
xy2c = (((-4.271 + -4.729)/2, (10.771+11.229)/2),
        ((0.229 + -0.229)/2, (15.729 + 15.271)/2),
        ((4.271 + 4.729)/2, (11.229 + 10.771)/2),
        (-3.500, 26.000), (3.500, 26.000))
        
R2 = sqrt((xy2c[0][0] - xy2[0][0])**2 + (xy2c[0][1] - xy2[0][1])**2)

Rext = sqrt((xy2c[0][0] - xy2[0][0])**2 + (xy2c[0][1] - xy2[0][1])**2)/2
Rint = Rext*CondutoresEspecs["Rail"][0]/CondutoresEspecs["Rail"][1]

CondutoresPos["Rail Compacto"] = (raio_eq(4, CondutoresEspecs["Rail"][1], R2),raio_eq(4, CondutoresEspecs["Rail"][1], R2),
                                 (xy2c[0][0], xy2c[1][0], xy2c[2][0]),
                                      (xy2c[0][1], xy2c[1][1], xy2c[2][1]),
                                      (xy2c[3][0], xy2c[3][1], xy2c[4][0], xy2c[4][1]), 4)

# =============== 500kV rail recapacitado =====================
#2 ultimos elementos sao pararraio, da forma (x, y)
n_rail_cap = 4
R_rail_cap1 = 0
R_rail_cap2 = 0

# raio do condutor equivalente
xy3 = ((-7.229, 8.639), (-7.229, 10.500), (-6.771, 10.500), (-6.771,
                                                             8.639), (-0.229, 10.500),
       (-0.229, 9.876), (0.229, 9.876), (0.229, 10.500), (7.229,
                                                          8.639), (7.229, 10.500), (6.771, 10.500),
       (6.771, 8.639), (-5.000, 20.500), (5.000, 20.500))

xy3c = (((-7.229 + -6.771)/2, (8.639 + 10.500)/2),
        (((-0.229 + 0.229)/2, (9.876 + 10.500)/2)),
        (((6.771 + 7.229)/2, (8.639 + 10.500)/2)),
        (-5.000, 20.500), (5.000, 20.500))

R3 = sqrt((xy3c[0][0] - xy3[0][0])**2 + (xy3c[0][1] - xy3[0][1])**2)

Rext = sqrt((xy3c[0][0] - xy3[0][0])**2 + (xy3c[0][1] - xy3[0][1])**2)/2
Rint = Rext*CondutoresEspecs["Rail"][0]/CondutoresEspecs["Rail"][1]


CondutoresPos["Rail Recapacitado"] = (raio_eq(4, CondutoresEspecs["Rail"][1], R3),raio_eq(4, CondutoresEspecs["Rail"][1], R3),
                                     (xy3c[0][0], xy3c[1][0], xy3c[2][0]),
                                      (xy3c[0][1], xy3c[1][1], xy3c[2][1]),
                                      (xy3c[3][0], xy3c[3][1], xy3c[4][0], xy3c[4][1]), 4)

# =================== 345kv ruddy ==========================
#2 ultimos elementos sao pararraio, da forma (x, y)
n_ruddy = 2
# raio do condutor equivalente é igual ao diametro
R_ruddy = 2*CondutoresEspecs["Ruddy"][1]
xy1 = ((-7.229, 10.5), (-6.771, 10.500), (-0.229, 10.500), (0.229,
                                                            10.500), (7.229, 10.500),
       (6.771, 10.500), (-5.000, 20.500), (5.000, 20.500))

xy1c = (((-7.229 + -6.771)/2, 10.5),
        ((-0.229 + 0.229)/2, 10.5),
        ((7.229 + 6.771)/2, 10.5),
        (-5.000, 20.500), (5.000, 20.500))


Rext = sqrt((xy1c[0][0] - xy1[0][0])**2 + (xy1c[0][1] - xy1[0][1])**2)/2
Rint = Rext*CondutoresEspecs["Ruddy"][0]/CondutoresEspecs["Ruddy"][1]
R1 = sqrt((xy1c[0][0] - xy1[0][0])**2 + (xy1c[0][1] - xy1[0][1])**2)

CondutoresPos["Ruddy"] = (raio_eq(2, CondutoresEspecs["Ruddy"][1], R1),raio_eq(2, CondutoresEspecs["Ruddy"][1], R1),
                                     (xy1c[0][0], xy1c[1][0], xy1c[2][0]),
                                      (xy1c[0][1], xy1c[1][1], xy1c[2][1]),
                                      (xy1c[3][0], xy1c[3][1], xy1c[4][0], xy1c[4][1]), 4)


#===================================================================================================
#========================== CALCULO DOS PARAMETROS PARA CADA CONDUTOR ==============================
#===================================================================================================


#(self, r_int, r_ext, nfase, npr, xc, yc, rhoc, rhoc_pr, rf, rpr):
nfase = 3
npr = 2
Vs = {"Bluejay": 750*1000,"Rail normal": 500*1000, "Rail convencional":500*1000, "Rail compacto":500*1000, "Rail recapacitado":500*1000,"Ruddy": 345*1000}
for i in CondutoresPos:
    
    if i == "Bluejay":
        nameSpec = "Bluejay"
    else:
        nameSpec = "Rail"
    #O Rail normal e diferente. Aparentemente e uma LT dupla
    if i == "Rail Normal":
        continue
        
    #(rext,rint,(Xc-1,Xc0,Xc1),(Yc-1,Yc,Yc+1),(Xpr,Ypr,Xpr2,Ypr2),n)
    r_ext = CondutoresPos[i][0]
    r_int = CondutoresPos[i][1]
    nfase = 3
    xc = np.concatenate((np.array(CondutoresPos[i][2]), np.array([CondutoresPos[i][4][0]]), np.array([CondutoresPos[i][4][2]])))
    yc = np.concatenate(( np.array(CondutoresPos[i][3]), np.array([CondutoresPos[i][4][1]]) , np.array([CondutoresPos[i][4][3]]) ))
    rhoc = CondutoresEspecs[nameSpec][2]
    rhoc_pr = CondutoresEspecs["3/8 EHS"][2]
    rf = r_ext
    rpr = CondutoresEspecs["3/8 EHS"][1]
    
    
    print("Linha: " + i)
    Linha = Linha_transmissao(r_int, r_ext, nfase, npr, xc, yc, rhoc, rhoc_pr, r_ext, rpr)
    Z = Linha.impedancia()
    Y = Linha.admitancia()
    Z = Z.astype(np.csingle)
    Y = Y.astype(np.csingle)
    #Retiramos as info do pararaio das matrizes atraves da reducao de kron
    Zabc = 0j + np.zeros((3,3))
    Yabc = 0j + np.zeros((3,3))
    
    Zabc = Z[0:nfase,0:nfase] - Z[0:nfase,nfase:] @ inv(Z[nfase:,nfase:]) @ Z[nfase:,0:nfase]
    Yabc = Y[0:nfase,0:nfase] - Y[0:nfase,nfase:] @ inv(Y[nfase:,nfase:]) @ Y[nfase:,0:nfase]


    #matrizes de sequencia
    a = np.exp(1j * np.deg2rad(120))
    a2 = a**2
    A = np.array([[1, 1, 1], [1, a2, a], [1, a, a2]])
    z012 = inv(A)@Zabc@A
    y012 = inv(A)@Yabc@A
    Pnat = Pnat(Vs[i], Vr, Ir)
    print("Impedancia Z+ = " + str(z012[1][1]))
    print("Impedancia Y+ = " + str(y012[1][1]))

