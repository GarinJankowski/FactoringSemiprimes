import ModCalc as Mod
import PrimeCalc as Prime
import FactorCalc as Factor
import MatrixCalc as Matrix
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


def remove_commas(numstr):
    return int(numstr.replace(",", ""))


if __name__ == '__main__':
    # calc = ModCalc.ModCalc(8)
    # calc.create_table("**", (2, 21), (0, 100))
    # print(euclid_extended(17, 47500))
    #run_primality_tests(6597072809)
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

    #print(FactorCalc.rational_sieve(187))
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
    #Matrix.Gauss_Jordan_elimination(matrix)
    new_matrices = Matrix.left_null_space(matrix)
    Matrix.print_matrix(Matrix.mod_matrix(new_matrices[0], 2))
    Matrix.print_matrix(Matrix.mod_matrix(new_matrices[1], 2))
    t = timer(t)
    #print(FactorCalc.rational_sieve(remove_commas("36,019")))
    #print(FactorCalc.rational_sieve(11*29))
    print(Factor.generate_quadratic_residue_primes(4, 1007))
    t = timer(t)
    #print(Factor.exp_quadratic_sieve(15347, 5))
    t = timer(t)
    #print(Factor.exp_quadratic_sieve(remove_commas("638,373,133"), 10))
    matrix = [[0, 1, 1, 0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 1, 1, 0, 1, 0], [0, 0, 0, 1, 1, 0, 0, 1, 0], [0, 1, 0, 0, 1, 0, 0, 0, 1], [0, 1, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0, 0, 1], [0, 0, 0, 1, 1, 0, 0, 1, 0], [0, 0, 0, 0, 1, 1, 0, 1, 0], [0, 1, 1, 0, 1, 1, 0, 0, 0], [0, 1, 1, 0, 1, 0, 1, 0, 1]]
    # print(Factor.all_possible_matrix_combinations(matrix))
    t = timer(t)
    Factor.Tonelli_Shanks(5, 41)
    t = timer(t)
    Factor.Tonelli_Shanks(71, 1009)