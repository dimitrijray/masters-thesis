
# coding: utf-8

# In[91]:

### Computing Deuring's correspondence ###
### Requires SageMath ###

# coding: utf-8
#set_random_seed(1) #For debugging
proof.all(False) #For debugging (2)

##Some formatting stuff.
centermod = lambda X,N: ZZ(X%N - N*floor((X%N)/((N+1)/2)))

import itertools
import operator as ops
import time
import random

#Next Powersmooth Number Generator
def smoothies(s,start = list()):
    psize = 1 #Index of Primes()
    isWayTooMuch = False
    #Initialize
    number = start
    if len(number) == 0:
        number.append([Primes()[psize],1])
        ChangeThis = 0
        MaximumToChange = 0
        yield number
    else:
        ChangeThis = 0
        psize = len(number)
        MaximumToChange = psize - 1
    while not(isWayTooMuch):
        if number[ChangeThis][0]^(number[ChangeThis][1]+1) > s:
            if ChangeThis == MaximumToChange:
                if Primes()[psize + 1] <= s:
                    psize += 1
                    MaximumToChange += 1
                    number.append([Primes()[psize],1])
                    #Reset everything!
                    for i in range(MaximumToChange):
                        number[i][1] = 0
                    ChangeThis = 0
                    yield number
                else:
                    isWayTooMuch = True
                    raise IndexError("Out of range! Raise your smoothness!")
            else:
                ChangeThis += 1
        else:
            number[ChangeThis][1] += 1
            if ChangeThis == 0:
                pass
            else:
                for i in range(ChangeThis):
                    number[i][1] = 0
                ChangeThis = 0
            yield number

def blend(smoothies):
    slices = list(itertools.starmap(ops.pow,smoothies))
    blitz = reduce(ops.mul,slices,1)
    return blitz


# In[92]:

#Inputs and precomputed things

#Primes and Numbers
p = 431
S1 = 4515 
S2 = 8948537162565

#Powersmooth bound, for tracking.
smoothness = (7/2)*log(p)


# In[93]:

#Get the quaternion algebra going first.
if p%4 != 3:
    print("Nope, dude! It's gotta be 3 mod 4!")
    exit()
B.<i,j,k> = QuaternionAlgebra(QQ,-1,-p);
OBasis = [1,i,(i+j)/2,(1+k)/2]
O0 = B.quaternion_order(OBasis)

#####Random Coefficient Generation#####
#Generate matrices with a specified determinant.
IWantThisDeterminant = 9
isIGotAnIdeal = False
count = 1
limit = IWantThisDeterminant^5
while not(isIGotAnIdeal) and (count <= limit):
    isIGotThisDeterminant = False
    while not(isIGotThisDeterminant):
        rd = [ZZ.random_element(-4,4) for _ in range(10)]
        M = matrix(ZZ,[
            [rd[0], rd[1], rd[2], rd[3]],
            [0, rd[4], rd[5], rd[6]],
            [0, 0, rd[7], rd[8]],
            [0, 0, 0, rd[9]]
        ]) #Assuming matrix is upper triangular ==> Otherwise can use elem. row ops and get an equiv. ideal.
        if M.determinant() == IWantThisDeterminant:
            isIGotThisDeterminant = True
    IBasis = [M[idx,0]*OBasis[0] + M[idx,1]*OBasis[1] + M[idx,2]*OBasis[2] + M[idx,3]*OBasis[3] for idx in range(4)]
    I = O0.left_ideal(IBasis)
    try:
        normI = I.norm()
        isIGotAnIdeal = True #Flag True first, check in the next lines.
        exitlist = list(itertools.product(OBasis,IBasis))
        for stuff in exitlist:
            check = stuff[0]*stuff[1]
            if not(check in I):
                isIGotAnIdeal = False
                count = count + 1
                break
            else:
                pass
    except AssertionError:
        count = count + 1
    if(isIGotAnIdeal):
        print("We have an ideal.")
assert I.left_order() == O0
#Start timing now!
starttime1 = time.time()
I = O0.left_ideal(IBasis)
assert I.left_order() == O0
normI = I.norm()

#Find the delta such that things get smoother. A priori set m = log(p), says the paper.
#Initialize counter
it = 0
m = -ceil(log(p))-1;
upperbound = ceil(log(p))
x = [m,0,m,0]
#N = p;
N = 1;
while not N in Primes():
    while x[0] <= upperbound:
        it = it + 1
        delta = x[0]*IBasis[0] + x[1]*IBasis[1] + x[2]*IBasis[2] + x[3]*IBasis[3]
        assert delta in O0.unit_ideal()
        N = ZZ(delta.reduced_norm() / normI)
        if N in Primes():
            break
        else:
            x[3] += 1 #A roundabout way to avoid 4 "for"s. Has a bonus of not checking the previous cases.
            if x[3] > upperbound:
                x[3] = m
                x[2] += 1
                if x[2] > upperbound:
                    x[2] = m
                    x[1] += 1
                    if x[1] > upperbound:
                        x[1] = m
                        x[0] += 1
        if it%100001 == 0:
            print("Iter:",it,"x:",x)
    upperbound = upperbound+10
N = ZZ(N) #Typecast to prevent errors.

#For checking purposes. Can comment out, along with the NO0 basis check.
#NO0Basis = [stuff * N for stuff in OBasis]
#NO0 = O0.left_ideal(NO0Basis)

assert delta in O0.unit_ideal()
#Compute I'
IIBasis = [(stuff * delta.conjugate())/normI for stuff in IBasis]
II = O0.left_ideal(IIBasis)
assert II.left_order() == O0
IIBasis = II.basis()

#Cornacchia time!
#Initialize c and d...
m = floor(sqrt((N*S1)/(2*p)));
c = -m-1
d = -m
#Initialize loop for the linear algebra.
isSolvable = False
solve_counter = 0
while not(isSolvable):
    solve_counter = solve_counter + 1
    isFound = False
    if c<m:
        c = c+1
    else:
        d = d+1
        c = -m
    while d<=m:
        rhs = ZZ(N*S1 - p*(c^2 + d^2))
        if rhs >=0:
            try:
                a,b = two_squares(rhs)
                isFound = True
            except ValueError:
                pass
        if not(isFound):
            if c<m:
                c = c+1
            else:
                d = d+1
                c = -m
        else:
            break
    if isFound:
        beta1 = a + b*i + c*j + d*k
        beta1 = O0(beta1)
        #Begin search for alpha such that gcd(norm(alpha),N^2)=N.
        isAlphaFound = False
        while not(isAlphaFound):
            RandomCoeff = [ZZ.random_element(-10,10) for _ in range(4)]
            alpha = sum([IIBasis[idx]*RandomCoeff[idx] for idx in range(4)])
            alpha = O0(alpha)
            if gcd(alpha.reduced_norm(),N^2) == N:
                isAlphaFound = True
        #Solve for beta2.
        cA = vector(alpha.coefficient_tuple())
        cB = beta1.coefficient_tuple()
        coeff_mat = matrix(GF(N),[
            [-cB[2]*p, -cB[3]*p, -cA[0], cA[1], cA[2]*p, cA[3]*p],
            [-cB[3]*p, cB[2]*p, -cA[1], -cA[0], -cA[3]*p, cA[2]*p],
            [cB[0], -cB[1], -cA[2], cA[3], -cA[0], -cA[1]],
            [cB[1], cB[0], -cA[3], -cA[2], cA[1], -cA[0]],
        ])
        solution = coeff_mat.right_kernel().basis()
        isFinallyFound = False
        vv = ZZ.random_element(1,N-1)*solution[0]
        CC = centermod(ZZ(vv[0]),N)
        DD = centermod(ZZ(vv[1]),N)
        beta2 = ZZ(CC)*j + ZZ(DD)*k
        unit = ZZ(vv[2]) + ZZ(vv[3])*i + ZZ(vv[4])*j + ZZ(vv[5])*k
        if (beta2.reduced_norm() % N == 0):
            pass
        else:
            isSolvable = True
            beta2 = ZZ(CC)*j + ZZ(DD)*k
            unit = ZZ(vv[2]) + ZZ(vv[3])*i + ZZ(vv[4])*j + ZZ(vv[5])*k
            #For checking purposes, uncomment if needed.
            #check = beta1 * beta2 - unit*alpha
            #print(check)
            #check_value = check in NO0
            #check_2 = unit.reduced_norm() % N
            #print("Is beta1 * beta2 = u*alpha (mod NO0?) {}".format(check_value))
            #print("check2: {}".format(check_2))
            #print(solve_counter)
        #Find beta2'such that beta2' = lambda * beta2.
        primedex = 1
        temp = p*(CC^2+DD^2)*ZZ(S2)
        legendre = kronecker(temp,N)
        if legendre == 1:
            isQuadraticResidue = True
            r = 1
        else:
            isQuadraticResidue = False
        while not(isQuadraticResidue):
            r = Primes()[primedex]
            temp = p*(CC^2+DD^2)*ZZ(S2)*ZZ(r)
            legendre = kronecker(temp,N)
            if legendre == 1:
                isQuadraticResidue = True
                print("Hey! Found the prime: {}".format(r))
                break
            else:
                primedex = primedex + 1
        S2 = S2*r
        lambsquared = ZZ(S2 * inverse_mod((CC^2+DD^2)*p,N))
        lamb = ZZ(mod(lambsquared,N).sqrt())
        #Cornacchia time! (2)
        check_count = 0
        limit = (2*N)/3 #You can change this, this is just to avoid the program running forever due to small S2.
        while not(isFinallyFound) and check_count<limit:
            check_count += 1
            #Grab the (c,d):
            c = centermod(ZZ.random_element(1,N-1),N)
            rhs_d = (S2 - p*(lamb^2)*(CC^2+DD^2) - 2*lamb*N*CC*c*p)*inverse_mod(2*lamb*DD*p,N^2)
            d = ZZ((rhs_d%(N^2))/N)
            d = centermod(d,N)
            #Set RHS for Cornacchia
            rhs = (S2 - p*(lamb*CC + c*N)^2 - p*(lamb*DD + d*N)^2)/(N^2)
            if rhs >= 0:
                try:
                    a,b = two_squares(rhs)
                    print("Finally! Found it!")
                    isFinallyFound = True
                except ValueError:
                    pass
            else:
                continue
        if isFinallyFound:
            a = N*a
            b = N*b
            c = lamb*CC + c*N
            d = lamb*DD + d*N
            beta2prime = a + b*i + c*j + d*k
            betas = beta1 * beta2prime
            assert betas in II
            JJBasis = [(stuff*betas.conjugate())/N for stuff in IIBasis]
            JJ = O0.left_ideal(JJBasis)
            assert JJ.left_order() == O0
            print(JJ)
        else:
            raise ValueError("Try increasing S2.")
    else:
        print("Tough luck, mate.  Try again!")
        break
#That's it! Stop timing this part!
endtime1 = time.time()


# In[94]:

#Curve declaration
Fq = GF(p).extension(2)
sqmin = Fq(-1).square_root() #The square root of -1 in the finite field
get1 = lambda X: X.curve()(X[0], X[1], X[2])
getI = lambda X: X.curve()(-X[0],sqmin*X[1],X[2])
getJ = lambda X: X.curve()(X[0]**p,X[1]**p,X[2]**p)
getK = lambda X: getI(getJ(X))
E0 = EllipticCurve(Fq,[1,0]) #Lift the curve to GF(p^2) to simplify things.


# In[95]:

def myownkernel(A):
#    M = M.change_ring(ZZ).stack(MatrixSpace(ZZ, 2)(M.base_ring().characteristic()))
#    return M.row_space().basis()
    #Sage doesn't handle kernels of matrices over a ring. This is a workaround...
    print A
    m, n = A.dimensions()
    l = A.base_ring().characteristic()
    M = ZZ^m / (l*ZZ^m)
    phi = M.hom([(ZZ^n/(l*ZZ^n))(a) for a in A])
    K = phi.kernel()
    return K, list(map(M, K.gens()))

def dlog2dimprimepow(P0, Q0, R):
    #2-dimensional discrete log.
    (l, e), = factor(P0.order())
    assert P0.order() == Q0.order() == l**e
    assert P0.order() % R.order() == 0
    try: (l**(e-1)*P0).discrete_log(l**(e-1)*Q0); raise Exception('not a basis')
    except ValueError: pass

    tab = dict()
    for A in range(l**e):
        tab[R - A*P0] = A
    for B in range(l**e):
        A = tab.get(B*Q0, None)
        if A is not None:
            assert A*P0 + B*Q0 == R
            return A, B
    raise ValueError("oops, couldn't find coeffs")

#Begin computing the kernel.
#Start timing this part now!
starttime2 = time.time()
JJfacs = factor(JJ.norm())
print JJfacs
#Prepare the list of isogeny kernels.
iso_kernels = list()
for fac in JJfacs:
    the_exp = fac[0]**fac[1]
    extdeg = next(i for i, d in enumerate(E0.count_points(the_exp), 1) if d % (2*the_exp**2) == 0)
    print("Computing a {}-torsion basis over GF({})...".format(the_exp, factor(p^(2*extdeg))))
    E0ext = E0.change_ring(Fq.extension(extdeg))
    E0ext.set_order(E0.count_points() if extdeg == 1 else E0.count_points(extdeg)[-1])

    #Compute basis points: P0 and Q0.
    while True:
        while True:
            P0 = E0ext.random_element()
            if P0.order() % the_exp: continue
            P0 *= P0.order() // the_exp
            assert P0.order() == the_exp
            break
        while True:
            Q0 = E0ext.random_element()
            if Q0.order() % the_exp: continue
            Q0 *= Q0.order() // the_exp
            assert Q0.order() == the_exp
            break
        print("got points of order {}.".format(the_exp))
        try:
            tmul = fac[0]**(fac[1]-1)
            discrete_log(tmul*P0,tmul*Q0,ZZ(fac[0]),operation='+')
            print("They're not independent, though.")
        except ValueError:
            break

    print("They're independent. Found it!")
    
    #Find the 2-division points.
    divP = P0.division_points(2)[0]
    divQ = Q0.division_points(2)[0]
    assert 2*divP == P0 and 2*divQ == Q0

    imsP = [f(divP) for f in (get1, getI, getJ, getK)]
    imsQ = [f(divQ) for f in (get1, getI, getJ, getK)]
    maizemat = matrix(Zmod(the_exp), 2, 0)
    #Mapping time!
    print("Looking for the matrix...")
    for basicstuff in JJBasis:
        assert basicstuff in O0.unit_ideal()
        base_without_denominators = 2 * basicstuff #Clear the denominators
        stuff = base_without_denominators.coefficient_tuple()
        multips = [ZZ(things) for things in stuff]
        multips = [ZZ(things) % (2*the_exp) for things in stuff]
        #Begin the mapping process.
        P1 = sum(m*R for m, R in zip(multips, imsP))
        Q1 = sum(m*R for m, R in zip(multips, imsQ))
        assert the_exp % P1.order() == 0 and the_exp % Q1.order() == 0
        #Find the coefficient matrix
        A, B = dlog2dimprimepow(P0, Q0, P1)
        C, D = dlog2dimprimepow(P0, Q0, Q1)
        #Looking for kernels.
        cobmat = matrix(Zmod(the_exp),[[A,B],[C,D]])
        maizemat = maizemat.augment(cobmat)
    kr0n, maize = myownkernel(maizemat)
    print("Kernel: {} (size {})".format(maize, kr0n.cardinality()))
    assert kr0n.cardinality() == the_exp
    #Plug all vectors in the kernel.
    for corn in maize:
        cereal = ZZ(corn[0])*P0 + ZZ(corn[1])*Q0
        iso_kernels.append(cereal)
    print("I'm done with this corn. Next...")
    #print("===============================================================")

for K in iso_kernels:
    print K.order()

#Isogeny time!
n_isogs = len(iso_kernels)
Etemp = E0
for idx in range(n_isogs): 
    grain = iso_kernels[idx]
    #Generate the kernel polynomial, whose roots are the x-coordinates of the isogeny kernel.
    kerpoly = prod(polygen(grain.curve().base_field()) - (i*grain).xy()[0] for i in range(1, 1 + grain.order()//2)).change_ring(Fq)

    for jdx in range(idx, n_isogs):
        #Temporarily consider the isogeny over an extension to push remaining kernel points through.
        F = iso_kernels[jdx].curve().base_field()
        iso_kernels[jdx] = Etemp.change_ring(F).isogeny(kerpoly.change_ring(F))(iso_kernels[jdx])

    Etemp = Etemp.isogeny_codomain(kerpoly)

#That's it! Stop timing!
endtime2 = time.time()
time_elapsed1 = abs(endtime1-starttime1)
time_elapsed2 = abs(endtime2-starttime2)
print(Etemp)
print("Time elapsed part 1: {} secs".format(time_elapsed1))
print("Time elapsed part 2: {} secs".format(time_elapsed2))

