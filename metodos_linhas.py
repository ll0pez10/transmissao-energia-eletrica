# função para calcular o raio equivalente
def raio_eq(n, rext, R):

    #R = raio do cicrulo que sera o novo condutor
    #rext = raio do feixe de condutor
    re = (rext*n*(R**(n-1)))**(n)
    return re

    re = (rext*n*(R**(n-1)))**(n)
    return re

# potencia caracteristica
def Pnat(Vs, Vr, Ir):
    Zc = Vr/Ir
    Pnat = (Vs**2)/Zc.real
    return Pnat
