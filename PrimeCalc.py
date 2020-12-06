import math
import random

import sympy


def isprime_simple(num):
    """
    checks if num is even
    then checks every integer below sqrt(num) to see if it is a divisor of num
    :param num: number to have its primality determined
    :return: bool if num is prime
    """
    if num < 2 or num % 2 == 0:
        return False
    sqrt = math.sqrt(num)
    i = 3
    while i < sqrt:
        if num % i == 0:
            return True
        i += 2
    return False


def isprime_6k(n: int) -> bool:
    """
    Code from  Introduction to Algorithms (Second ed.). MIT Press and McGraw–Hill. pp. 887–896.
    """
    """Primality test using 6k+-1 optimization."""
    if n <= 3:
        return n > 1
    elif n % 2 == 0 or n % 3 == 0:
        return False
    i: int = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def GCD_Stein(a, b):
    """
    Stein's Algorithm, a slightly different take on Euclid's algorithm
    gives the greatest common divisor of the params
    Works with large numbers thanks to floor division
    :param v: positive integer
    :param u: positive integer
    :return: GCD of v and u
    """

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
            return 2 * GCD_Stein(a // 2, b // 2)
        # if only one is even, let's say 'a', then their gcd is the same as the gcd of a/2 and b
        else:
            return GCD_Stein(a // 2, b)
    elif b % 2 == 0:
        return GCD_Stein(a, b // 2)
    # if neither are even, there is one more thing you can check:
    # their gcd should be the same as the gcd of half their difference and the smallest value
    else:
        if a > b:
            return GCD_Stein((a - b)//2, b)
        else:
            return GCD_Stein((b - a)//2, a)

def GCD_Stein_while(a, b):
    """
    Stein's Algorithm, a slightly different take on Euclid's algorithm
    gives the greatest common divisor of the params
    Works with large numbers thanks to floor division
    :param v: positive integer
    :param u: positive integer
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


def isprime_base_case(n):
    """
    checks if a number is:
        less than 2 (return False)
        equal to 2 or 3 (return True)
        divisible by 2 or 3 (return False)
    these bools are represented as ints
    returns 2 if the number doesn't satisfy any of these
    we can cut out a bunch of numbers like this

    :param n: number to be tested
    :return: an int based on what this method detected
    """
    if n <= 3:
        return n > 1
    elif n % 2 == 0 or n % 3 == 0:
        return False
    return None


def isprime_sympy(n):
    """
    documentation link for sympy's ntheory.primetest:
    http://www.caacle.com/sympy-docs-html-1.4/_modules/sympy/ntheory/primetest.html
    :param n: number to be tested
    :return: bool if n is prime
    """
    return sympy.isprime(n)


def isprime_Fermat(n, iterations=5):
    """
    og Fermat probabilistic primality test, gets destroyed by Carmichael numbers
    :param n: number to be tested
    :param iterations: number of times to run the test, higher number means higher accuracy
    :return: bool if n is a probable prime
    """
    base_case = isprime_base_case(n)
    if base_case is not None:
        return base_case

    for i in range(iterations):
        a = random.randint(2, n - 2)
        if mod_pow(a, n - 1, n) == 1:
            return True
    return False


def isprime_Miller_Rabin(n, iterations=5):
    pass


def generate_primes(max_limit=100.0, min_limit=0, amount=0):
    """
    generates an amount of primes starting at the min_limit
    this can be specified by either the maximum prime or the amount of primes generated
    :param min_limit: inclusive minimum number allowed to be tested
    :param max_limit: inclusive maximum number allowed to be tested
    :param amount: amount of primes to be generated
    :return: list of primes
    """
    primes = []
    if min_limit > max_limit:
        max_limit = 0

    if min_limit <= 2:
        primes.append(2)

    i = min_limit
    if i % 2 == 0:
        i += 1

    # basically just goes through every odd number and does a primality test using sympy
    while (i <= max_limit or max_limit == 0) and (amount == 0 or len(primes) < amount):
        if sympy.isprime(i):
            primes.append(i)
        i += 2
    return primes


def sieve_of_Eratosthenes(max_num):
    """
    prime generator using the Sieve of Eratosthenes
    this is a bit slower than the above method

    :param max: maximum number allowed to be generated
    :return: list of primes under the max
    """
    # takes a list of booleans that are True from 2 to max
    nums = [True for i in range(max_num)]
    nums[0] = False
    nums[1] = False
    i = 0
    sqrt = math.sqrt(max_num)
    while i < sqrt:
        # goes up the list, finding the first True number
        # this True number is found to be prime
        if nums[i]:
            # after finding a True number, mark all multiples of this number in the list False
            for k in range(i+1, max_num):
                if nums[k] and k % i == 0:
                    nums[k] = False
        i += 1
        # continue until the end of the list is reached and only the primes are left True

    # return all of the True values
    primes = []
    for i in range(max_num):
        if nums[i]:
            primes.append(i)
    return primes
