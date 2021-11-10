

def Factor(n: int):
    # poly f
    # natural m, f(m)=0 (mod n)

    # complex theta, f(theta)=0
    # poly ring in theta over integers: Ztheta
    # MAPPING phi from Ztheta to Z/nZ,
    #       + and * retained
    #       phi(1)=1 (mod n)
    #       phi(theta)=m (mod n)

    # beta^2 in Ztheta
    # y^2 in Z
    # CONGRUENCE OF SQUARES
    #       find beta^2 = product(a+b(theta)) for a bunch of (a,b) pairs
    #       find y^2 = product(a+b(m)) for those same pairs
    #       then (phi(beta))^2 = y^2 (mod n)
    
    
    # PROCESS
    # R: rational fb (quad residues mod n?? idk)
    # A: algebraic fb, bunch of (a+b(theta)) in Ztheta
    #       primality: same thing, each has no poly divisors in Ztheta
    #       REPRESENTATION:
    #               (r, p)
    #                   r in Z/nZ, f(r)=0 (mod p)
    #                   prime
    #               pairs of (r,p) is bijective with (a+b(theta)) in Ztheta
    #                   so like you can use this instead I guess
    
    # gotta find (a,b) pairs such that their poly is smooth over BOTH fbs
    #       (a+b(theta)) over A, (a+b(m)) over R
    # SOME THEOREMS
    # 1. (r,p) for (c+d(theta)), (c+d(theta))|(a+b(theta)) <-> a=-br (mod p)
    # 2. what: the set U of pairs (r,p) completely factorizes (a+b(theta)) <-> product(all p in U) = (-b)^d * f(-a/b), d = degree of f
    # 3. prime q | (a+bm) <-> a=-bm (mod q)

    # OKAY
    # so this time we're seiving over TWO lists, crazy right
    #       fix b, a varies from -N to N
    #       list1 = (a+b(theta)), list2 = (a+b(m))
    #       

    pass