
import ModCalc as Mod
import PrimeCalc as Prime
import FactorCalc as Factor
import MatrixCalc as Matrix
import QuadraticSieve as QS
import QSHelper
import sympy
import time


DEFAULT_ITERATIONS = 5

def euclid_extended(firstkey, totient):
    inverse = 0
    tempinverse = 1
    tot = totient
    key = firstkey
    while key != 0:
        quotient = int(tot/key)

        newtempinverse = inverse - (quotient*tempinverse)
        (inverse, tempinverse) = (tempinverse, newtempinverse)

        newkey = tot - (quotient*key)
        (tot, key) = (key, newkey)
    if tot > 1:
        return "no inverse found"
    if inverse < 0:
        inverse += totient
    return inverse


def timer(start=time.time()):
    print("------", time.time()-start, "------")
    return time.time()


def run_primality_tests(number):
    print("Number:", number)
    t = timer()
    print("Sympy: ", Prime.isprime_sympy(number))
    t = timer(t)
    print("Fermat: ", Prime.isprime_Fermat(number, DEFAULT_ITERATIONS))
    t = timer(t)
    # print(Prime.isprime_Fermat(number, DEFAULT_ITERATIONS))
    # t = timer(t)
    print("\r")


def get_private_key(pub, factors):
    totient = (factors[0]-1)*(factors[1]-1)
    priv = QSHelper.extended_euclidean_algorithm(pub, totient)
    return priv


def remove_commas(numstr):
    return int(numstr.replace(",", ""))


if __name__ == '__main__':
    t = timer()
    qs = QS.QuadraticSieve()

    splist = [
        974387 * 493351,
        756023 * 509227,
        20809 * 809023,
        950959783 * 532136131,
        731004667 * 95282753,
        1696845001 * 3254172751,
        956424366058369634576511997901,
        845687515256677794057145483142257,
        940562294173322057894192280217862479,
        1050499240906169867321081277725674006861,
        1010327588654592479489518715626174153572611,
        954408344258307963167386933095470926757906833,
        1206641416126754032441387494108570591309330651103,
        973607460392457781571715947566555017699035751061097,
        1259714232554985299730573571992530137288967940148051,
        1312146826492950870824970931914634678701075368845546689397,
        1196512127911231582816317581340375183596021484875123829415441,
    ]

    for sprime in splist:
        qs.factor(sprime)
        t = timer(t)

