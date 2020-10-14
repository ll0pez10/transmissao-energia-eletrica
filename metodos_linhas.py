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
