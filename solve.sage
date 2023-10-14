from Crypto.Util.number import *

def small_roots(f, bound,r,s,N,m):
    '''
    small private d RSA with moduli N=p^r*q^s,that d < 1-(3*r+s)/(r+s)^2 - eps
    the eps is 
    ((15*s+1)*r^4-(2*s^2-10*s)*r^3-(s^3-6*s^2+8*s)*r^2+(2*s^3-12*s^2+6*s)*r+s^4-4*s^3+s^2)
    /(4*m*(r-s)*(r+s)^3) in theory,by this we can choose the proper m which is lower than theory.
    return : one factor of N.you can factor N by this
    for example:
    p,q : 256
    r,s : 5,3
    d : 1300,m = 9
    d : 1350,m = 14
    d : 1360,m = 20
    d : 1370,m = 30
    d > 1370,m = 40 # spend long time to do this(2800s)
    you can choose the larger m to approach the theory solution
    '''
    t1 = int((r*(r+s-2))/((r-1)*(r+s))*m)
    t2 = m
    bounds = [bound ,1]
    f = f.change_ring(ZZ)
    G = Sequence([], f.parent())
    x = f.variables()[0]
    for k in range(t2+1):
        for i in range(t2+1-k):
            d=max([0,ceil((r-1)*(t1-k)/r),ceil((s-1)*(t2-k)/s)])
            base=N ^ d * f ^ k * x ^ i
            G.append(base)
    B, monomials = G.coefficient_matrix()
    monomials = vector(monomials)
    factors = [monomial(*bounds) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)
    B = B.dense_matrix().LLL()
#     B = flatter(B)
    '''
    another question is can't use flatter because flatter not support 
    the matrix that its row far greater than its cols
    '''
    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1 / factor)
    H = Sequence([], f.parent().change_ring(QQ))
    for h in filter(None, B * monomials):
        for i in h.coefficients():
            if gcd(i,N)!=1 and gcd(i,N)!=N:
                return gcd(h.coefficients()[0],N)
    return 0

def test():
    r,s=5,3
    p,q=getPrime(256),getPrime(256)
    N = p^r*q^s
    phi=p^(r-1)*q^(s-1)*(p-1)*(q-1)
    edge=1350
    d=getPrime(edge)
    e=ZZ(inverse(d,phi))
    a= -int(inverse(e,N)) %N
    PR.<x,y> = PolynomialRing(Zmod(N))
    f=a-x
    m=13
    res=small_roots(f,2^edge,r,s,N,m)
    print(res)
c = 5029723337007818808832754776759893887444879662801349712576810470262312396754408718170770707925015017385645584442512000797339170647864326108311590562844171
N = 3759675650502563695153151584917824507854926661732907939028389279774124406013042838546953581576073208871094857148508406802633489701504134249499597914542124463779252098837657394974239335884726202413862183404321654463922171057401752386333809584179109491691851156579178452248534783401568199131101637952819524702798663560557483668765136026315271419314073036116879338555727777506342904334411379055692051548735524394491162124016936020960278594934832905479002480160251814228947417553571414430429706482587579876807513469399180281619097937606917846073419200176183340717568702116796677513739961828976919478654457637451826294933
e = 229008556409573164371978380452021106688588097729158338468248213917888175178492067982525826582582124769664067336853325793830736524069028120877995293076501685071845478613891814063550943140421907014666043731432899819999451833948917119826997988890110827467830371163637151737633259525347140559518939332719377544459045827417076412346696030494989054139185610912700663962872185889165738547511312375915148946228556442812249150675925829256431970480221003802721711007352304700373376185262322264329636611961422545471143841872105226762691690138343686371843086125054384530349905361116692007292869579146352767879179902925722367519
r,s=5,3
a= -int(inverse(e,N)) %N
PR.<x,y> = PolynomialRing(Zmod(N))
f=a-x
m=40
res=small_roots(f, 2^edge,r,s,N,m)
res # 65669608693032794298247369306859402326551240348853293641688013177984363140773
import gmpy2
q=65669608693032794298247369306859402326551240348853293641688013177984363140773
p=int(gmpy2.iroot(N//q^3,5)[0])
assert N==p^5*q^3
d=inverse(65537,(p-1)*(q-1))
long_to_bytes(ZZ(pow(c,d,p*q)))
#b'DASCTF{8f4ff7c3-4017-455b-9773-371ab6f42a63}'