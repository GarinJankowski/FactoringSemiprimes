import copy
import math

import sympy

import PrimeCalc
import MatrixCalc as Matrix

"""
format for prime factors in a dictionary:
key = prime factor
value = power of that factor

so, 2^3 5^1 7^2 = {
2: 3,
5: 1,
7: 2
}
"""


def sympy_factor_semiprime(n):
    """
    documentation link for sympy's factor/factorint:
    http://man.hubwiz.com/docset/SymPy.docset/Contents/Resources/Documents/_modules/sympy/ntheory/factor_.html
    :param n: number to be factored
    :return: a tuple containing the two prime factors of n, or 0 if n is not a semi prime
    """
    dict = sympy.factorint(n)

    if len(dict) == 2:
        tup = []
        for key in dict:
            if key == 1:
                tup.append(dict[key])
            elif dict[key] == 1:
                tup.append(key)
            else:
                return 0, 0
        return tuple(tup)
    return 0, 0


def basic_factoring(n):
    """
    factors a number by going through a list of primes and checking them
    :param n: number to be factored
    :return: dictionary of factors (format at the top of this file)
    """
    primes = PrimeCalc.generate_primes(n/2)
    factors = {}
    i = 0
    num = n
    while i < len(primes) and num > 1:
        prime = primes[i]
        while num % prime == 0:
            num /= prime
            if prime in factors:
                factors[prime] += 1
            else:
                factors[prime] = 1
        i += 1
    if len(factors) == 0:
        factors[n] = 1
    return factors


def rational_sieve(n, bound=7):
    """
    uses the rational sieve to find the factors of a semiprime
    this thing kinda sucks, maybe its my fault though
    :param n: semiprime to be factored
    :param bound: the limit used to generated the factor base
    :return: pair of prime factors as a list
    """
# STEP 1
    # choose a bound
    # no clue what to choose for a bound, 7 as default for now I guess
    # maybe increase the bound for larger numbers
    print("NEW BOUND:", bound)

    # don't forget to check for prime powers here

# STEP 2
    # generates a list of primes up to the bound, inclusively, called the factor base
    # good to note that generate_primes uses sympy since their primality test is a super fast Rabin-Miller/B-PSW
    factor_base = PrimeCalc.generate_primes(max_limit=bound)
    # print(factor_base)

# STEP 3 (actual sieve part)
    # look for an integer z where z and z+n are both smooth to the bound
    # the amount of relations, (len(z_list)*(len(z_list)-1))/2, is the amount of comparisons we can make from the list.
    # this doesn't include self comparisons, which are a nice bonus,
    # although maybe I should count them as comparisons so this doesn't take forever.
    # len(factor_base)+1 is basically just a one more than the size of the factor base,
    # which should be sufficient for how many z combinations we want since it ensures a linear dependency.
    # I also set a limit to how many can be searched ((bound**2)*1000) because sometimes you just have to give up
    # after giving up, the sieve reaches the end and tries a new bound
    # I've learned this speeds up the process somewhat
    i = 2
    z_list = []
    while (len(z_list)*(len(z_list)-1))/2 <= len(factor_base)+1 and i < bound**2*1000:
        # you can use your own basic_factoring instead of sympy's rho/p-1 stuff if you wanna ramp up the inefficiency
        z = i
        z_factors = sympy.factorint(z)
        # z_factors = basic_factoring(z)

        nz = z + n
        nz_factors = sympy.factorint(nz)
        # nz_factors = basic_factoring(nz)

        # if both z and z+n are (bound)-smooth, they are congruent in mod (bound) and so they are added as a valid pair
        # z, nz, and dictionaries of their prime factors exist as a tuple of tuples: ((z, z_factors), (nz, nz_factors))
        if list(z_factors.keys())[-1] <= bound and list(nz_factors.keys())[-1] <= bound:
            z_list.append(((z, z_factors), (nz, nz_factors)))
            # print(i)
            # print(z_list)

        i += 1

# STEP 4
    # multiply all the combinations of two pairs with each other
    # check congruence of squares between the products and added them if it works
    # this also checks if a pair itself is a congruence of squares
    congruent_squares = find_congruence_of_squares(z_list)
    # print("Congruences:", congruent_squares)

# STEP 5
    # check gcd factorization to get the factors
    # then check if they are non-trivial
    prime_factors = get_factors_from_congruence_of_squares(n, congruent_squares)

    # if no nontrivial factors were found, increase the bound and restart because why not
    if len(prime_factors) == 0:
        return rational_sieve(n, bound+10)
    return prime_factors


def quadratic_sieve(n, base_size=7):
    """
    MUCH faster than the rational sieve
    but just as inconsistent (I need to optimize it)

    This time, the sieve is a bit more particular.
    The factor base only consists of primes that make n a quadratic residue (mod prime) to remove unnecessary primes.
    Checks a bunch of integers starting at 0 by putting them into a polynomial and getting all the y(x) as a list.
    Next it uses a sieve of Eratosthenes to get the numbers from that list that factor only into primes from the base.
    Then it puts all those y(x) into a matrix, each row being the exponents of the prime factors of a y(x).
    Gauss-Jordan elimination is then performed on the matrix augmented with an identity to find its left null space.
    The resulting identity from the elimination tells us which combinations of rows produce a perfect square.
    We then get a congruence of squares from the relation of the square to a product of the x values of the square's factors.
    Finally, we can pop those roots into Euclid's algorithm and hopefully get some non-trivial factors.
    If not, increment the size of the factor base and try again.

    Beware: its so unoptimized that it can rerun too many times, go past the max recursion depth, and throw and error.
    In reality this should actually only have to run once. We'll see if we can get there.

    :param n: semiprime to be factored
    :param base_size: size of the factor base
    :return: non-trivial factors of semiprime n
    """
    # nice little loading indicator for ya
    loading = "."
    if (base_size-3) % 50 == 0:
        loading += "\r"
    print(loading, sep=' ', end='', flush=True)
# STEP 1
    # choose a base_size instead of a bound this time
    # this is because we want a number of primes, not all primes below a maximum
    # this is due to the fact that we want quadratic residue primes,
    # so an amount is probably better than a maximum since idk residue distribution
    # default base_size is in the parameters

# STEP 2
    # for this sieve, you can't just use any old primes for factor base
    # for each prime p, n must be a quadratic residue mod p
    # aka, there must exist an integer x where x^2 = n (mod p)
    # this cuts out any unnecessary primes so we can search on a more relevant set
    factor_base = generate_quadratic_residue_primes(base_size, n)
    # print("Bound:", base_size)
    print(factor_base)

# STEP 3
    # go through a certain range of X values to put into a polynomial
    # idk how to choose the right polynomial so I got a default one, some people even use multiple at once
    # list the resulting y values to get a base for your sieve
    # range for the input is another variable that we can mess with, but I have no clue what to actually make it
    # consult the Fermat-Kraitchik method, maybe some clues there
    sqrt_n = math.sqrt(n)
    if not sqrt_n.is_integer():
        sqrt_n = int(sqrt_n+1)
    y_list = []
    for i in range(500):
        # POLYNOMIAL
        # try to make it better and tailor it to the number being factored
        value = int(math.pow(i + sqrt_n, 2) - n)
        y_list.append(value)
    # print("y list:", y_list)

# STEP 4
    # sieve your y values for numbers where all their prime factors are in the factor base
    # this gives us a list of numbers where we can try all of their product combinations to find a perfect square
    # maybe I could utilize Eratosthenes better here because I'm not really using it
    y_list_copy = copy.deepcopy(y_list)
    z_list = []
    parity_matrix = []
    for i in range(len(y_list)):
        # maybe also optimize this with Shanks-Tonelli
        for prime in factor_base:
            while y_list_copy[i] % prime == 0:
                y_list_copy[i] /= prime
        if y_list_copy[i] == 1:
            index = i + sqrt_n
            # value = int(math.pow(i + sqrt_n, 2) - n)
            value = y_list[i]
            z_list.append((index, value))
            parity_matrix.append(prime_factors_to_parity(sympy.factorint(value), factor_base))
    i_index = [0 for i in factor_base]
    i_increment_index = [1 for i in factor_base]
    y_size = len(y_list)
    # while len(i_index) > 0:
    #     i = 0
    #     while i < len(i_index):
    #         if y_list_copy[i_index[i]] % factor_base[i] == 0:
    #             if i_increment_index[i] == 1:
    #                 i_increment_index[i] = factor_base[i]
    #         while y_list_copy[i_index[i]] % factor_base[i] == 0:
    #             y_list_copy[i_index[i]] /= factor_base[i]
    #
    #         if y_list_copy[i_index[i]] == 1:
    #             index = i_index[i] + sqrt_n
    #             # value = int(math.pow(i + sqrt_n, 2) - n)
    #             value = y_list[i_index[i]]
    #             tup = index, value
    #             if tup not in z_list:
    #                 z_list.append(tup)
    #                 parity_matrix.append(prime_factors_to_parity(sympy.factorint(value), factor_base))
    #
    #         i_index[i] += i_increment_index[i]
    #         if i_index[i] >= y_size:
    #             i_index.pop(i)
    #         else:
    #             i += 1

    print("y list:", y_list)
    print("y list:", y_list_copy)
    print("z list:", z_list)
    Matrix.print_matrix(parity_matrix)
    if len(z_list) == 0:
        # return quadratic_sieve(n, base_size+1)
        return None

# STEP 5
    # gets all possible combinations of the matrix's rows that form a zero matrix
    # does this by getting the left null space of the matrix
    # also for some reason all_possible_matrix_combinations keeps the answers from the last loop of the sieve,
    # so I have to specify that the list of answers should start out empty

    # left_null_space = Matrix.left_null_space(parity_matrix)
    # parity_matrix = Matrix.mod_matrix(left_null_space[0], 2)
    # new_identity_matrix = Matrix.mod_matrix(left_null_space[1], 2)
    # print()
    # Matrix.print_matrix(parity_matrix)
    # print()
    # Matrix.print_matrix(new_identity_matrix)
    # valid_combos = []
    #
    # for r in range(len(parity_matrix)):
    #     if sum(parity_matrix[r]) == 0:
    #         valid_combos.append(new_identity_matrix[r])
    valid_combos = all_possible_matrix_combinations(parity_matrix)

    print("valid combos:", valid_combos)
    if len(valid_combos) == 0 or len(valid_combos[0]) > len(z_list):
        # return quadratic_sieve(n, base_size+1)
        return None

# STEP 6
    # tries all the combos by multiplying them out and putting them into Euclid's algorithm
    # returns any non-trivial factors
    for combo in valid_combos:
        # squares of indices
        square_a = 1
        # squares of polynomial results
        square_b = 1
        for i in range(len(combo)):
            square_a *= z_list[i][0]
            square_b *= z_list[i][1]
        square_b = int(math.sqrt(square_b))
        # print(combo, square_a, square_b, "\r")
        # print(square_a, square_b)
        factor1 = GCD_Stein(math.fabs(square_a - square_b), n)
        factor2 = 0
        if factor1 != 1 and factor1 != n:
            print("Factor base size:", base_size)
            return factor1, n/factor1
    # return quadratic_sieve(n, base_size+1)
    return None


# don't use this anymore
def all_possible_matrix_combinations(og_matrix, index=0, combo=[], answers=[]):
    if index == len(og_matrix):
        if sum(combo) == 0:
            return
        matrix = copy.deepcopy(og_matrix)
        answer = Matrix.multiply_row(matrix[0], combo[0])
        for i in range(1, len(combo)):
            row_to_add = Matrix.multiply_row(matrix[i], combo[i])
            answer = Matrix.add_rows(answer, row_to_add)
        works = True
        for parity in answer:
            if parity % 2 != 0:
                works = False
        if works:
            answers.append(combo)
        return

    combo1 = list(combo)
    combo1.append(0)
    all_possible_matrix_combinations(og_matrix, index+1, combo1, answers)

    combo2 = list(combo)
    combo2.append(1)
    all_possible_matrix_combinations(og_matrix, index+1, combo2, answers)

    return answers


# helper boys
def get_congruence_of_squares(a, b):
    """
    simultaneously checks if each list of prime factors is a perfect square and
    attempts to make them both perfect squares if they aren't
    note these params are meant to be congruent on whatever mod, so any operation on one param has to be done on the other
    this is why it checks of the odds lists are equivalent, because if they aren't, then they won't be congruent any more
    :param factorsA: prime factor dictionary
    :param factorsB: prime factor dictionary
    :return: both resulting squares as a tuple, or returns empty if they aren't a congruence of squares
    """
    factorsA = dict(a)
    factorsB = dict(b)
    # decrements each prime factor with an odd power and records that factor to oddsA
    oddsA = []
    keysA = list(factorsA.keys())
    for n in keysA:
        e = factorsA[n]
        if e % 2 != 0:
            oddsA.append(n)
            factorsA[n] -= 1
            if factorsA[n] == 0:
                del factorsA[n]

    # does the same for b
    oddsB = []
    keysB = list(factorsB.keys())
    for n in keysB:
        e = factorsB[n]
        if e % 2 != 0:
            oddsB.append(n)
            factorsB[n] -= 1
            if factorsB[n] == 0:
                del factorsB[n]

    # if all prime factors have an even power and each side was evenly simplified,
    # then the list of odds should be the same
    # this means we got ourselves a congruence of squares lets go
    # also the resulting numbers from the simplification have to be NOT zero
    new_squares = ()
    if set(oddsA) == set(oddsB) and len(factorsA) > 0 and len(factorsB) > 0:
        for n in factorsA:
            factorsA[n] /= 2
        for n in factorsB:
            factorsB[n] /= 2
        new_squares = prime_factors_to_int(factorsA), prime_factors_to_int(factorsB)
    return new_squares


def find_congruence_of_squares(z_list):
    """
    takes a list of tuples of tuples
    each tuple contains a tuple that contains a number and a dictionary of its prime factors
    these paired tuples should be smooth on some unknown bound, it doesn't really matter here
    it takes these pairs and sees if it can get a congruence of squares from the products of two pairs
    it also checks each pair itself for a conguence of squares, cuz you never know right
    :param z_list: list of tupled tuples. contains pairs of numbers and their prime factors
    :return: a list containing tupled tuples of congruences of squares
    """
    congruent_squares = []
    if len(z_list) > 1:
        for a in range(len(z_list)):
            factors1 = z_list[a]
            z_factors1 = factors1[0][1]
            nz_factors1 = factors1[1][1]

            # checks if the relation itself is a congruence of squares
            squares = get_congruence_of_squares(z_factors1, nz_factors1)
            if squares:
                congruent_squares.append(squares)

            # multiply with each relation ahead in the z_list and check if they are congruent squares
            for b in range(a + 1, len(z_list)):
                factors2 = z_list[b]
                z_combined = multiply_prime_factors(z_factors1, factors2[0][1])
                nz_combined = multiply_prime_factors(nz_factors1, factors2[1][1])

                # checks if the multiplied relations are a congruence of squares
                squares = get_congruence_of_squares(z_combined, nz_combined)
                if squares:
                    congruent_squares.append(squares)
    return congruent_squares


def get_factors_from_congruence_of_squares(n, congruent_squares):
    """
    takes the tupled tuples that fulfill the congruence of squares
    attemps to get two non-trivial factors a and b of n by applying the identity (a + b)(a - b) = 0 (mod n)
    :param n: semiprime n
    :param congruent_squares: list of tupled tuples that are a congruence of squares and a dict of their prime factors
    :return: two non-trivial factors of n, or an empty list
    """
    prime_factors = []
    for pair in congruent_squares:
        # formula for producing the factors
        gcdA = GCD_Stein(math.fabs(pair[1] - pair[0]), n)
        gcdB = GCD_Stein(pair[0] + pair[1], n)
        # if these factors are nontrivial, then we DID IT (theoretically)
        if gcdA != 1 and gcdA != n and gcdB != 1 and gcdB != n:
            prime_factors = [gcdA, gcdB]
            return prime_factors
    return prime_factors


def Legendre_symbol(a, p):
    """
    calculates the Legendre symbol of a number and prime using Euler's criterion
    :param a: number to be tested
    :param p: prime that goes with it
    :return: Legendre symbol (1, 0, or -1)
    """
    power = int((p-1)/2)
    legendre = mod_pow(a, power, p)
    if legendre == p-1:
        legendre = -1
    return legendre


def generate_quadratic_residue_primes(amount, n):
    """
    generates an amount of primes where n is a quadratic residue (mod prime)
    :param amount: amount of primes to be generated
    :param n: number that is tested to be a quadratic residue
    :return: list of primes
    """
    primes = [2]
    # basically just goes through every odd number and:
    #   does a primality test using sympy
    #   does a quadratic residue test using Euler's criterion
    i = 3
    while len(primes) < amount:
        if sympy.isprime(i) and Legendre_symbol(n, i) == 1:
            primes.append(i)
        i += 2
    return primes


def is_perfect_square(n):
    return math.sqrt(n).is_integer()


def is_perfect_square_factors(factor_dict):
    """
    determines if a number is a perfect square from their prime factors
    idk if this is gonna be used
    :param factor_dict:
    :return:
    """
    for key in factor_dict:
        if factor_dict[key] % 2 != 0:
            return False
    return True


def is_smooth(a, n):
    """
    checks if a number 'a' is n-smooth with the help of sympy factoring
    :param a: number to be tested
    :param n: smooth value
    :return: bool in a is n-smooth
    """

    greatest_factor = list(sympy.factorint(a).keys())[-1]
    if greatest_factor <= n:
        return True
    return False


def GCD_Stein(a, b):
    """
    go look in PrimeCalc
    """
    return PrimeCalc.GCD_Stein(a, b)


def prime_factors_to_parity(factor_dict, primes):
    vector = []
    for prime in primes:
        if prime in factor_dict:
            vector.append(factor_dict[prime]%2)
        else:
            vector.append(0)
    return vector


def prime_factors_to_int(factor_dict):
    """
    multiplies a number's prime factors together from their prime factor dict
    :param factor_dict: prime factor dict
    :return: int value of the number
    """
    num = 1
    for key in factor_dict:
        num *= math.pow(key, factor_dict[key])
    return num


def multiply_prime_factors(dict1, dict2):
    """
    multiply two prime factors dicts together
    could be useful for the sieve
    :param dict1: factor dict
    :param dict2: factor dict
    :return: product dict
    """
    new_dict = dict(dict2)
    for key in dict1:
        if key in new_dict:
            new_dict[key] += dict1[key]
        else:
            new_dict[key] = dict1[key]
    return new_dict


def mod_pow(base, exp, mod):
    """
    modular exponentiation, iterative
    took from GeeksforGeeks
    """

    answer = 1
    base %= mod

    while exp > 0:
        if exp & 1:
            answer = answer * base % mod
        exp >>= 1
        base = base * base % mod
    return answer
