import sympy


def mod(number, modulo):
    if modulo != 0:
        number %= modulo
    return number


def mod_pow(base, exp, mod):
    """
    modular exponentiation, iterative
    taken from GeeksforGeeks
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


def extended_euclidean_algorithm(d, mod):
    """
    :return: modular multiplicative inverse of firstkey in mod(totient)
    :author: William Retert
    """
    inverse = 0
    tempinverse = 1
    m = mod
    nd = d
    while nd != 0:
        quotient = m//nd

        newtempinverse = inverse - (quotient*tempinverse)
        (inverse, tempinverse) = (tempinverse, newtempinverse)

        new_nd = m - (quotient*nd)
        (m, nd) = (nd, new_nd)
    if m > 1:
        return "no inverse found"
    if inverse < 0:
        inverse += mod
    return inverse


def GCD_Stein_while(a, b):
    """
    Stein's Algorithm, a slightly different take on Euclid's algorithm
    gives the greatest common divisor of the params
    Using while loop instead of recursion because recursion limit
    :param a: positive integer
    :param b: positive integer
    :return: GCD of v and u
    """

    multiplier = 1
    while True:
        # if a is the same as b, their gcd is themselves
        if a == b:
            return a
        # if one number is 0, their gcd is the other number
        if a == 0:
            return b
        if b == 0:
            return a

        # if a and b are both even, their gcd is the gcd of a/2 and b/2 then multiplied by 2
        if a % 2 == 0:
            if b % 2 == 0:
                multiplier *= 2
                a //= 2
                b //= 2
            # if only one is even, let's say 'a', then their gcd is the same as the gcd of a/2 and b
            else:
                a //= 2
        elif b % 2 == 0:
            b //= 2
        # if neither are even, there is one more thing you can check:
        # their gcd should be the same as the gcd of half their difference and the smallest value
        else:
            if a > b:
                a = (a-b)//2
            else:
                b = (b-a)//2


# ----------------------------------
#      Quadratic Residue section
# ----------------------------------


def Legendre_symbol(a, p):
    """
    calculates the Legendre symbol of a number and prime using Euler's criterion
    :param a: number to be tested
    :param p: prime that goes with it
    :return: Legendre symbol (1, 0, or -1)
    """
    a = a % p
    if a == 0:
        return 0
    power = (p-1)//2
    legendre = mod_pow(a, power, p)
    if legendre == p-1:
        legendre = -1
    return legendre


def generate_quadratic_residue_primes(amount, n, minimum=0, direction=1):
    """
    generates an amount of primes where n is a quadratic residue (mod prime)
    :param amount: amount of primes to be generated
    :param n: number that is tested to be a quadratic residue
    :return: list of primes
    """
    if minimum <= 2:
        primes = [2]
        i = 3
    else:
        primes = []
        i = int(minimum)
        if i % 2 == 0:
            i += direction
    # basically just goes through every odd number and:
    #   does a primality test using sympy
    #   does a quadratic residue test using Euler's criterion
    while len(primes) < amount:
        if sympy.isprime(i) and Legendre_symbol(n, i) == 1:
            primes.append(i)
        i += 2*direction
    return primes


def Tonelli_Shanks(n, p, Dickson=0):
    """
    calculates the integer r that fulfills n being a quadratic residue mod p in r^2 = n (mod p), if any
    thanks Wikipedia

    :param n: quadratic residue
    :param p: prime modulus of the quadratic residue
    :param Dickson: prime power exponent of p, if any
    :return: a list containing any r that fulfills r^2 = n (mod p), returns None if n is a non-residue
    """
    # added a case for 2 so I don't have to do that in other places
    if p == 2:
        squares = [n % 2]
        if Dickson > 0:
            squares = Tonelli_Dickson(squares, n, p, Dickson)
        return squares

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
    R = mod_pow(n, (Q+1)//2, p)

    if t == 0:
        return [0]
    elif t == 1:
        # also returns -R because if R works mod p, then so does -R mod p
        # same for the return at the end of this method
        squares = [R, p-R]
        if Dickson > 0:
            squares = Tonelli_Dickson(squares, n, p, Dickson)
        return squares
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
        squares = [R, p-R]
        if Dickson > 0:
            squares = Tonelli_Dickson(squares, n, p, Dickson)
        return squares
    return []


def Tonelli_Dickson(squares, n, p, d):
    """
    Dickson said that Tonelli gave a formula for solving squares of
    quadratic residues in the mod of a prime power so here it is I guess.
    hope it works

    :param squares: Tonelli-Shanks list of n (mod p)
    :param n: quadratic residue
    :param p: prime root
    :param d: power
    :return: list of roots r where r^2=n(mod p^d)
    """
    new_squares = []
    for r in squares:
        powmod = p ** d
        powmod2 = p ** (d - 1)
        nsq = (mod_pow(r, powmod2, powmod) * mod_pow(n, (powmod - (2 * powmod2) + 1) // 2, powmod)) % powmod
        new_squares.append(nsq)
    return new_squares


# ------------------------------------------
#      Prime Factor Dictionary section
#
#  (deals with prime factors as a dictionary
#  with the primes as keys and their exponents
#                 as values)
# ------------------------------------------


def prime_factors_to_int(factor_dict):
    """
    multiplies a number's prime factors together from their prime factor dict
    :param factor_dict: prime factor dict
    :return: int value of the number
    """
    num = 1
    for key in factor_dict:
        num *= key**factor_dict[key]
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


# -------------------------
#      Matrix section
# -------------------------


def left_null_space(matrix, modulus=0):
    """
    augments the matrix with an identity matrix and performs Gauss-Jordan elimination
    this records the steps of the elimination onto the identity
    returns the resulting matrix and identity
    :param matrix: matrix to be solved
    :return: resulting matrix and identity in a tuple
    """
    rows = len(matrix)
    columns = len(matrix[0])
    identity = generate_identity(rows)
    for i in range(rows):
        matrix[i] += identity[i]

    Gauss_Jordan_elimination(matrix, modulus)

    for i in range(rows):
        whole_row = list(matrix[i])
        matrix[i] = whole_row[0:columns]
        identity[i] = whole_row[columns: len(whole_row)]

    return matrix, identity


def generate_identity(size):
    identity = []
    for r in range(size):
        identity_row = []
        for c in range(size):
            if c == r:
                identity_row.append(1)
            else:
                identity_row.append(0)
        identity.append(identity_row)
    return identity


def Gauss_Jordan_elimination(matrix, modulus=0):
    """
    puts a matrix in reduced row echelon form
    beware of numerical instability, idk what to do about it
    although hopefully it doesn't matter too much since we're only really using this mod 2
    :param matrix: matrix to be solved
    :return matrix: solved matrix
    """
    rows = len(matrix)
    columns = len(matrix[0])
    r = 0
    c = 0
    while r < rows and c < columns:
        # find the row index of the number in column c with the greatest absolute value
        # this finds the location of the 'c'th pivot
        r_pivot = r
        for i in range(r, rows):
            n = abs(matrix[i][c])
            if n >= abs(matrix[r_pivot][c]):
                r_pivot = i
        # if there is no pivot in this column, aka the whole column is zeros, then skip the column
        if matrix[r_pivot][c] != 0:
            pivot = matrix[r_pivot][c]
            # divide the pivot row by the pivot to make it 1
            for c2 in range(columns):
                matrix[r_pivot][c2] = mod(matrix[r_pivot][c2]/pivot, modulus)
            swap_rows(matrix, r, r_pivot)
            # set everything in the column to zero except for the pivot
            # also subtract the pivot row from each row that gets modified
            for r2 in range(rows):
                if r2 != r:
                    quotient = matrix[r2][c]/matrix[r][c]
                    matrix[r2][c] = 0
                    for c2 in range(c+1, columns):
                        matrix[r2][c2] = mod(matrix[r2][c2]-matrix[r][c2]*quotient, modulus)
            r += 1
        c += 1
    return matrix


def swap_rows(matrix, a, b):
    temp_a = list(matrix[a])
    matrix[a] = list(matrix[b])
    matrix[b] = temp_a


def mod_matrix(matrix, modulus):
    """
    puts a matrix in a specified modulus
    :param matrix: matrix to be modded
    :param modulus: modulus
    :return: modded matrix
    """
    for r in range(len(matrix)):
        for c in range(len(matrix[r])):
            matrix[r][c] %= modulus
    return matrix