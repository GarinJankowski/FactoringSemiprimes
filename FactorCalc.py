import copy
import math

import sympy

import PrimeCalc
import MatrixCalc as Matrix
import ModCalc as Mod

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


def quadratic_sieve(n):
    """
    MUCH faster than the rational sieve
    but just as inconsistent (I need to optimize it)

    This time, the sieve is a bit more particular.
    The factor base only consists of primes that make n a quadratic residue (mod prime) to remove unnecessary primes.
    Checks a bunch of integers starting at 0 by putting them into a polynomial and getting all the y(x) as a list.
    Next it uses a sieve of Eratosthenes to get the numbers from that list that factor uniquely into primes from the base.
    Then it puts all those y(x) into a matrix, each row being the exponents of the prime factors of a y(x) in mod 2.
    Gauss-Jordan elimination is then performed on the matrix augmented with an identity to find its left null space.
    The resulting identity from the elimination tells us which combinations of rows produce a perfect square.
    We then get a congruence of squares from the relation of the square to a product of the x values of the square's factors.
    Finally, we can pop those roots into Euclid's algorithm and hopefully get some non-trivial factors.

    Beware: This thing is pretty unoptimized, and by that I mean it doesn't work most of the time.
    We'll see if we can get some better polynomials going.

    :param n: semiprime to be factored
    :param base_size: size of the factor base
    :return: non-trivial factors of semiprime n
    """
    # nice little loading indicator for ya
    # loading = "."
    # if (base_size-3) % 50 == 0:
    #     loading += "\r"
    # print(loading, sep=' ', end='', flush=True)
    print(n, "\nBits:", n.bit_length())
# STEP 1
    # choose a base_size instead of a bound this time
    # this is because we want a number of primes, not all primes below a maximum
    # this is due to the fact that we want quadratic residue primes,
    # so an amount is probably better than a maximum since idk residue distribution

    # alright so APPARENTLY according to this paper https://www.cs.virginia.edu/crab/QFS_Simple.pdf,
    # the optimum size of the factor base is (e**(ln(n)*ln(ln(n)))**(sqrt(2)/4) I guess
    # and the optimum size of the sieving interval is that number cubed
    # I finally figured out why the sieving is gonna take awhile
    scary_calculation = (math.e**(math.sqrt(math.log(n)*math.log(math.log(n)))))**(math.sqrt(2)/4)
    base_size = int(scary_calculation)
    sieving_interval = int(scary_calculation**3)
    print("B:", base_size, "\nM:", sieving_interval)


# STEP 2
    # for this sieve, you can't just use any old primes for factor base
    # for each prime p, n must be a quadratic residue mod p
    # aka, there must exist an integer x where x^2 = n (mod p)
    # this cuts out any unnecessary primes so we can search on a more relevant set
    factor_base = generate_quadratic_residue_primes(base_size, n)
    # print(factor_base)

# STEP 3
    # go through a certain range of X values to put into a polynomial
    # idk how to choose the right polynomial so I got a default one, some people even use multiple at once
    # list the resulting y values to get a base for your sieve
    # range for the input is another variable that we can mess with, but I have no clue what to actually make it
    # consult the Fermat-Kraitchik method, maybe some clues there
    sqrt_n = math.sqrt(n)
    if not sqrt_n.is_integer():
        sqrt_n = int(sqrt_n+1)
    # the list to hold our values that we're gonna sieve on in STEP 4
    y_list = {}
    # goes every number under the interval, also covering the negative side, to get a range of -interval to interval
    for i in range(sieving_interval):
        # POLYNOMIAL
        # try to make it better and tailor it to the number being factored
        y_list[i] = int(math.pow(i + sqrt_n, 2) - n)
        y_list[-i] = int(math.pow(-i + sqrt_n, 2) - n)
    # print("y list:", y_list)

# STEP 4
    # sieve your y values for numbers where all their prime factors are in the factor base and also unique (power of 1)
    # this gives us a list of numbers where we can try all of their product combinations to find a perfect square

    # make a copy of the y_list so we can divide its values
    y_list_copy = copy.deepcopy(y_list)
    # the z_list is the resulting pairs of (index+sqrt_n) and (y_list[index]) that will give us our congruence factors
    z_list = []
    # this is a matrix of all y_list values that make it into z_list so we can see which values make a perfect square
    # basically is a prime factor exponent matrix (mod 2) of the factors of the values of the z_list,
    # I am explaining this poorly
    parity_matrix = []

    # use the Tonelli-Shanks algorithm on each prime to locate the y_list[indices] which are a multiple of the prime
    # Tonelli-Shanks basically gives us the values which make n a quadratic residue (mod prime),
    # which will work as indices for the y_list (I don't know why but it does)
    # this locates the first 1 or 2 indices, and every index+(x*prime) after these will also be a multiple of the prime
    #   that x variable in the comment just means every integer from 0 onwards
    #   this works like a sieve of Eratosthenes, since we know the next index will also work
    factor_indices_p = [Tonelli_Shanks(n, p) for p in factor_base]
    factor_indices_n = copy.deepcopy(factor_indices_p)
    for i in range(len(factor_indices_p)):
        for j in range(len(factor_indices_p[i])):
            factor_indices_p[i][j] = (factor_indices_p[i][j] - sqrt_n) % factor_base[i]
            factor_indices_n[i][j] = factor_indices_p[i][j] - factor_base[i]
    # print(factor_indices_p)
    # print(factor_indices_n)
    # a copy of the factor base where we can eliminate values when needed so the main while loop eventually ends
    factor_base_copy = list(factor_base)
    # this loop runs until every prime has gone through the y_list_copy, dividing everything they need to
    while len(factor_base_copy) > 0:
        p = 0
        # this next loop goes through each prime, doing a few things:
        #   divides the number at the factor_indices
        #   check if that number = 1, if so add it to the z_list and parity_matrix
        #   adds the current prime to each of its factor_indices, advancing its position in the y_list
        #       this works since we know that the value at y_list[index+prime] is going to be divisible by the prime
        while p < len(factor_base_copy):
            # does the next while loop for both the positive and negative direction of x,
            # since the sieving interval covers the negative side too
            current_indices = factor_indices_p
            negative = 1
            for j in range(2):
                # this last while loop should run a maximum of 2 times,
                # since there is a max of 2 possible results from Tonelli_Shanks()
                i = 0
                while i < len(current_indices[p]):
                    factor_i = current_indices[p][i]
                    y_list_copy[factor_i] /= factor_base_copy[p]
                    # if the divided value finally equals 1, then it is both smooth and unique over our factor base,
                    # so it is a viable value to be added to the z_list/parity_matrix
                    if y_list_copy[factor_i] == 1:
                        # the z_list receives a tuple containing the index+sqrt_n and the original y_list value
                        # this is why we made a y_list_copy to modify, because we still want the original, non-divided value
                        index = factor_i + sqrt_n
                        value = y_list[factor_i]
                        tup = index, value
                        # not sure if we'll ever run into duplicates, always nice to make sure
                        if tup not in z_list:
                            z_list.append(tup)
                            # the matrix receives a new row,
                            # equivalent to the parity of the exponents of the number's prime factors from the factor base
                            # basically, each value in the row corresponds to the prime in the factor base,
                            # where it is 1 if the prime is a factor and 0 if it is not
                            # NOTE: its likely faster to record the prime factors during the loop instead of using sympy
                            #   We can fix this later
                            parity_matrix.append(prime_factors_to_parity(sympy.factorint(value), factor_base))
                        else:
                            print("THIS CONDITIONAL IS USEFUL")
                    # advance this y_list index to the next one
                    current_indices[p][i] += factor_base_copy[p]*negative
                    # if the factor index has gone through all the its y_list values, remove it
                    if abs(current_indices[p][i]) >= sieving_interval:
                        current_indices[p].pop(i)
                    else:
                        i += 1
                current_indices = factor_indices_n
                negative = -1
            # if there are no more indices for this prime, it is complete
            # remove the prime from the factor_base_copy and remove its indices from factor_indices
            if len(factor_indices_p[p]) == 0 and len(factor_indices_n[p]) == 0:
                factor_indices_p.pop(p)
                factor_indices_n.pop(p)
                factor_base_copy.pop(p)
                p -= 1
            p += 1

    # print("y list:", y_list)
    # print("y list:", y_list_copy)
    print("z list:", z_list)
    Matrix.print_matrix(parity_matrix)
    if len(z_list) == 0:
        return None

# STEP 5
    # gets all possible combinations of the matrix's rows that form a zero matrix
    # does this by getting the left null space of the matrix
    # FOR SOME REASON it doesn't always find all the combinations?
    #   if you do n = 15247 and manually set M = 10,000, it misses a valid combo, idk why
    left_null_space = Matrix.left_null_space(parity_matrix, 2)
    parity_matrix = Matrix.mod_matrix(left_null_space[0], 2)
    new_identity_matrix = Matrix.mod_matrix(left_null_space[1], 2)
    valid_combos = []
    # Matrix.print_matrix(parity_matrix)
    # Matrix.print_matrix(new_identity_matrix)

    for r in range(len(parity_matrix)):
        if sum(parity_matrix[r]) == 0:
            valid_combos.append(new_identity_matrix[r])

    print("valid combos:", valid_combos)
    if len(valid_combos) == 0 or len(valid_combos[0]) > len(z_list):
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
        print(combo, square_a, square_b, "\r")
        # print(square_a, square_b)
        factor1 = GCD_Stein(math.fabs(square_a - square_b), n)
        factor2 = 0
        print(factor1)
        if factor1 != 1 and factor1 != n:
            return factor1, n/factor1
    return None


def exp_quadratic_sieve(n, base_size=5):
    """
    experimental
    this is just for messing around with the sieve

    :param n: semiprime to be factored
    :param base_size: size of the factor base
    :return: non-trivial factors of semiprime n
    """

    factor_base = generate_quadratic_residue_primes(base_size, n)
    print("Bound:", base_size)
    print(factor_base)

    sqrt_n = math.sqrt(n)
    if not sqrt_n.is_integer():
        sqrt_n = int(sqrt_n+1)
    # y_list = []
    # for i in range(1000):
    #     # POLYNOMIAL
    #     # try to make it better and tailor it to the number being factored
    #     value = int(math.pow(i + sqrt_n, 2) - n)
    #     y_list.append(value)
    # # print("y list:", y_list)
    #
    # y_list_copy = copy.deepcopy(y_list)
    # z_list = []
    # parity_matrix = []
    # for i in range(len(y_list)):
    #     # maybe also optimize this with Shanks-Tonelli
    #     for prime in factor_base:
    #         if y_list_copy[i] % prime == 0:
    #             y_list_copy[i] /= prime
    #     if y_list_copy[i] == 1:
    #         index = i + sqrt_n
    #         # value = int(math.pow(i + sqrt_n, 2) - n)
    #         value = y_list[i]
    #         z_list.append((index, value))
    #         parity_matrix.append(prime_factors_to_parity(sympy.factorint(value), factor_base))

    i = 0
    z_list = []
    parity_matrix = []
    max_num = 1
    value = 0
    for prime in factor_base:
        max_num *= prime
    print("Max:", max_num)
    # this sieve doesn't stop until it gets the parity matrix into a linear dependency
    # there is also a maximum where the sieve knows it won't find any more values,
    # which is the product of the factor base, so it just stops
    while value <= max_num and len(z_list) <= len(factor_base):
        y = int(math.pow(i + sqrt_n, 2) - n)
        value = y
        for prime in factor_base:
            if y % prime == 0:
                y /= prime
        if y == 1:
            index = i + sqrt_n
            # value = int(math.pow(i + sqrt_n, 2) - n)
            z_list.append((index, value))
            parity_matrix.append(prime_factors_to_parity(sympy.factorint(value), factor_base))
        i += 1

        # loading
        if i % 200000 == 0:
            loading = "."
            if i % 10000000 == 0:
                print(z_list)
                loading += "\r"
            print(loading, sep=' ', end='', flush=True)


    # print("y list:", y_list)
    # print("y list:", y_list_copy)
    print("z list:", z_list)
    Matrix.print_matrix(parity_matrix)
    if len(z_list) == 0:
        return None

    left_null_space = Matrix.left_null_space(parity_matrix, 2)
    parity_matrix = Matrix.mod_matrix(left_null_space[0], 2)
    new_identity_matrix = Matrix.mod_matrix(left_null_space[1], 2)
    valid_combos = []

    for r in range(len(parity_matrix)):
        if sum(parity_matrix[r]) == 0:
            valid_combos.append(new_identity_matrix[r])

    print("valid combos:", valid_combos)
    if len(valid_combos) == 0 or len(valid_combos[0]) > len(z_list):
        return None

    for combo in valid_combos:
        # squares of indices
        square_a = 1
        # squares of polynomial results
        square_b = 1
        for i in range(len(combo)):
            square_a *= z_list[i][0]
            square_b *= z_list[i][1]
        square_b = int(math.sqrt(square_b))
        print(combo, square_a, square_b, "\r")
        # print(square_a, square_b)
        factor1 = GCD_Stein(math.fabs(square_a - square_b), n)
        factor2 = 0
        print(factor1)
        if factor1 != 1 and factor1 != n:
            return factor1, n/factor1
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
def Tonelli_Shanks(n, p):
    """
    calculates the integer r that fulfills n being a quadratic residue mod p in r^2 = n (mod p), if any
    thanks Wikipedia
    :param n: quadratic residue
    :param p: prime modulus of the quadratic residue
    :return: a list containing any r that fulfills r^2 = n (mod p), returns None if n is a non-residue
    """
    # added a case for 2 so I don't have to do that in other places
    if p == 2:
        return [n % 2]

    p1_factors = sympy.factorint(p-1)
    S = p1_factors[2]
    Q = 1
    for key in p1_factors:
        if key != 2:
            Q *= key**p1_factors[key]

    z = -1
    i = 0
    while z == -1 and i < p:
        i += 1
        if Legendre_symbol(i, p) == -1:
            z = i

    M = S
    c = mod_pow(z, Q, p)
    t = mod_pow(n, Q, p)
    R = mod_pow(n, (Q+1)/2, p)

    if t == 0:
        return [0]
    elif t == 1:
        # also returns -R because if R works mod p, then so does -R mod p
        # same for the return at the end of this method
        return [R, p-R]
    i = 1
    b = 0
    while t != 1:
        i = 1
        while mod_pow(t, 2**i, p) != 1 and i < M:
            i += 1
        # added this break to account for quadratic non-residues that will never reach an answer
        if i >= M:
            break
        b = mod_pow(c, 2**(M-i-1), p)
        M = i
        c = mod_pow(b, 2, p)
        t = t*c % p
        R = R*b % p
    if t == 1:
        return [R, p-R]
    return []


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
    if a == 0:
        return 0
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

    base = int(base)
    exp = int(exp)

    while exp > 0:
        if exp & 1:
            answer = answer * base % mod
        exp >>= 1
        base = base * base % mod
    return answer
