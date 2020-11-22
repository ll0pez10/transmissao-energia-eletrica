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

def derivacao(FC, Y, Bc, ltot):
    Bc = Y*ltot
    Bl = -j*FC*Bc
    Ql = [1, 0, Bl, 1]
    return Ql

#teste

#funcoes uteis para achar o quadripolo da linha em função da distancia e em seguida a compensação da linha em função da distancia e da taxa que queremos compensar
def quadlinha(name,L):
    if name == "Bluejay":
        Z1 = z012_bluejay[1][1]/Zbaseblue*L #impedancia considerando a sequencia positiva
        Y1 = y012_bluejay[1][1]*Zbaseblue*L #admitancia considerando a sequencia negativa
    elif name == "Rail":
        Z1 = z012_rail[1][1]/Zbaserail*L #impedancia considerando a sequencia positiva
        Y1 = y012_rail[1][1]*Zbaserail*L #admitancia considerando a sequencia negativa

    #Componentes da matriz de quadripolo
    A1 = 1 + Z1*Y1/2
    A2 = Z1
    A3 = Y1*(1 + Z1*Y1/4)
    A4 = 1 + Z1*Y1/2

    return np.array( [[A1, A2],[A3, A4]] ) #Quadripolo da linha bluejay

def compenslinha(name,L,taxa):
    if name == "Bluejay":
        Y1 = y012_bluejay[1][1]*Zbaseblue*L
    elif name == "Rail":
        Y1 = y012_rail[1][1]*Zbaserail*L
        
    return np.array( [[1, 0],[-taxa*Y1, 1]] ) #se pa temos que ajustar esse 50% dps, a verificar de acordo com os resultados
