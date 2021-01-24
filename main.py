
import ModCalc as Mod
import PrimeCalc as Prime
import FactorCalc as Factor
import MatrixCalc as Matrix
import QuadraticSieve as QS
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
    priv = Factor.extended_euclidean_algorithm(pub, totient)
    return priv


def remove_commas(numstr):
    return int(numstr.replace(",", ""))


class One:
    def __init__(self):
        self.my_list = []
        self.two = Two(self.my_list)
        self.two.add_stuff()
        self.print_stuff()

    def print_stuff(self):
        print(self.my_list)


class Two:
    def __init__(self, my_list):
        self.__my_list = my_list

    def add_stuff(self):
        for i in range(20):
            self.__my_list.append(i)


if __name__ == '__main__':
    # calc = ModCalc.ModCalc(8)
    # calc.create_table("**", (2, 21), (0, 100))
    # print(euclid_extended(17, 47500))
    # run_primality_tests(6597072809)
    run_primality_tests(remove_commas("1,872,160,488,359,313,005,195,657,104,852,822,711,178,010,796,939,601,138,105,193,545,347,015,515,588,641,251,945,858,024,657,444,530,188,818,971,438,716,066,174,295,004,419,071,702,079,031,242,851,266,630,488,578,181,241,742,724,921,449,063"))
    run_primality_tests(101)
    num = 118627064770095767481642687937131474625503554938783849161211
    num = remove_commas("69,773,257,921,009,844,045,424,780,911")
    num = 13
    t = timer()
    print("Number:", num)
    print("Bit Length:", num.bit_length())
    print("Factors:", sympy.factorint(num))
    t = timer(t)
    print(Factor.basic_factoring(6784))
    print(Prime.generate_primes(min_limit=600, amount=7))

    # print(FactorCalc.rational_sieve(187))
    t = timer(t)
    print(Prime.sieve_of_Eratosthenes(1000)[-10:-1])
    t = timer(t)
    print(Prime.generate_primes(1000)[-10:-1])
    t = timer(t)
    num = remove_commas("3,256,628,491")
    num = 53*61
    # print("Semiprime:", num)
    # print("Bit Length:", num.bit_length())
    # print("Factors:", Factor.rational_sieve(num))
    # t = timer(t)
    matrix = [[-7, -3, 3, 12], [2, 2, 2, 0], [-1, -4, 3, -9]]
    matrix = [[5, 2, 0, 2], [2, 1, -1, 0], [2, 3, -1, 3]]
    # matrix = [[0, 4, -2, -2], [3, 4, -1, -6], [0, -2, 10, -8]]
    # matrix = [[3, 4, -1, -6], [0, -2, 10, -8], [0, 4, -2, -2]]
    matrix = [[2, 1, -3, -10], [0, -2, 1, -2], [0, 0, 1, 6]]
    matrix = [[1, -2, 1], [2, -3, 5], [1, 1, 0]]

    # this one has infinitely many solutions
    # matrix = [[1, -1, -1, 4], [2, -2, -2, 8], [5, -5, -5, 20]]

    matrix = [[1, 1, 1, 0], [0, 0, 0, 1], [1, 1, 1, 1]]
    Matrix.print_matrix(matrix)
    # Matrix.Gauss_Jordan_elimination(matrix)
    new_matrices = Matrix.left_null_space(matrix)
    Matrix.print_matrix(Matrix.mod_matrix(new_matrices[0], 2))
    Matrix.print_matrix(Matrix.mod_matrix(new_matrices[1], 2))
    t = timer(t)
    # print(FactorCalc.rational_sieve(remove_commas("36,019")))
    # print(FactorCalc.rational_sieve(11*29))
    print(Factor.generate_quadratic_residue_primes(4, 1007))
    t = timer(t)
    print(Factor.quadratic_sieve(15347))
    t = timer(t)
    # print(Factor.exp_quadratic_sieve(remove_commas("638,373,133"), 10))
    # matrix = [[0, 1, 1, 0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 1, 1, 0, 1, 0], [0, 0, 0, 1, 1, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0, 0, 0, 1], [0, 1, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0, 1, 1, 0, 0, 1, 0], [0, 0, 0, 0, 1, 1, 0, 1, 0], [0, 1, 1, 0, 1, 1, 0, 0, 0], [0, 1, 1, 0, 1, 0, 1, 0, 1]]
    # print(Factor.all_possible_matrix_combinations(matrix))
    print(Factor.PQS(87463))
    t = timer(t)
    #print(Factor.MPQS(remove_commas("316,805,816,423,168,663")))
    t = timer(t)
    #print(Factor.quadratic_sieve(remove_commas("316,805,816,423,168,663")))
    t = timer(t)
    splist = [
        # 953*263,
        # 211*863,
        # 31*139,
        # 887*557,
        974387*493351,
        756023*509227,
        20809*809023
    ]
    # for sprime in splist:
    #      #sp = remove_commas(sprime)
    #      print(Factor.PQS(sprime))
    #      # t = timer(t)
    #      #print(len(Factor.generate_polynomial(sprime)))
    #      t = timer(t)
    #      #print(Factor.quadratic_sieve(sprime))
    #      # t = timer(t)
    #     print(Factor.step_by_step(sprime))
    #     t = timer(t)
    #print(len(Factor.generate_polynomial(887*557)))
    #print(Factor.PQS(502613*360181))
    #print(Factor.PQS(887*557))

    # print(Factor.step_by_step(502613*360181))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 282))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 306))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 310))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 312))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 330))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 346))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 348))
    # t = timer(t)
    # print(Factor.step_by_step(502613*360181, 358))

    # print(Factor.step_by_step(887*557))
    # print(Factor.step_by_step(887*557))
    # print(Factor.step_by_step(remove_commas("421,357,885,598,003")*remove_commas("597,374,069,461,879")))
    # t = timer(t)
    # print(Factor.step_by_step(remove_commas("380,978,779,300,724,935,452,034,511,917")*remove_commas("1,133,840,581,612,979,458,722,383,048,567")))
    # t = timer(t)

    splist = [
        974387 * 493351,
        756023 * 509227,
        20809 * 809023,
        950959783 * 532136131,
        731004667 * 95282753,
        1696845001 * 3254172751
    ]

    qs = QS.QuadraticSieve()
    for sprime in splist:
        print(qs.factor(sprime))
        t = timer(t)

    splist = [
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
        print(qs.factor(sprime))
        t = timer(t)

