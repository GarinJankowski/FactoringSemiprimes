"""
- CS2911 - 0NN
- Fall 2017
- Lab N
- Names:
  -
  -

B-bit RSA
"""

import random
import sys
import math

# To change the number of bits, change BIT_LENGTH, MAX_PRIME, and MIN_PRIME.
# MAX_PRIME and MIN_PRIME should have exactly BIT_LENGTH/2 bits
BITS_PER_HEX_DIGIT = 4 # binary digits per hex digit -- always 4
BIT_LENGTH = 38  # "B" in the lab handout. Length of n in bits
HEX_LENGTH = (BIT_LENGTH+(BITS_PER_HEX_DIGIT-1))//BITS_PER_HEX_DIGIT # Length of n in hexadecimal digits
MAX_PRIME = 0b1111111111111111111  # The maximum value a prime number can have
MIN_PRIME = 0b1100000000000000001  # The minimum value a prime number can have
PUBLIC_EXPONENT = 17  # The default public exponent


def main():
    """ Provide the user with a variety of encryption-related actions """

    # Get chosen operation from the user.
    while True:
        action = input("Select an option from the menu below:\n"
                       "(1-CK) create_keys\n"
                       "(2-CC) compute_checksum\n"
                       "(3-VC) verify_checksum\n"
                       "(4-EM) encrypt_message\n"
                       "(5-DM) decrypt_message\n"
                       "(6-BK) break_key\n "
                       "Please enter the option you want:\n")
        # Execute the chosen operation.
        if action in ['1', 'CK', 'ck', 'create_keys']:
            create_keys_interactive()
        elif action in ['2', 'CC', 'cc', 'compute_checksum']:
            compute_checksum_interactive()
        elif action in ['3', 'VC', 'vc', 'verify_checksum']:
            verify_checksum_interactive()
        elif action in ['4', 'EM', 'em', 'encrypt_message']:
            encrypt_message_interactive()
        elif action in ['5', 'DM', 'dm', 'decrypt_message']:
            decrypt_message_interactive()
        elif action in ['6', 'BK', 'bk', 'break_key']:
            break_key_interactive()
        else:
            print("Unknown action: '{0}'".format(action))
        print("\r")


def create_keys_interactive():
    """
    Create new public keys

    :return: the private key (d, n) for use by other interactive methods
    """

    key_pair = create_keys()
    pub = get_public_key(key_pair)
    priv = get_private_key(key_pair)
    print("Public key: ")
    print(pub)
    print("Private key: ")
    print(priv)
    return priv


def compute_checksum_interactive():
    """
    Compute the checksum for a message, and encrypt it
    """

    priv = create_keys_interactive()

    message = input('Please enter the message to be checksummed: ')

    hash = compute_checksum(message)
    print('Hash:',as_hex(hash))
    cipher = apply_key(priv, hash)
    print('Encrypted Hash:', as_hex(cipher))


def verify_checksum_interactive():
    """
    Verify a message with its checksum, interactively
    """

    pub = enter_public_key_interactive()
    message = input('Please enter the message to be verified: ')
    recomputed_hash = compute_checksum(message)

    string_hash = input('Please enter the encrypted hash (in hexadecimal): ')
    encrypted_hash = int(string_hash, 16)
    decrypted_hash = apply_key(pub, encrypted_hash)
    print('Recomputed hash:', as_hex(format(recomputed_hash)))
    print('Decrypted hash: ', as_hex(format(decrypted_hash)))
    if recomputed_hash == decrypted_hash:
        print('Hashes match -- message is verified')
    else:
        print('Hashes do not match -- has tampering occured?')


def encrypt_message_interactive():
    """
    Encrypt a message
    """

    message = input('Please enter the message to be encrypted: ')
    pub = enter_public_key_interactive()
    encrypted = ''
    for c in message:
        ciphernumber = apply_key(pub,ord(c))
        encrypted += as_hex(ciphernumber)
    print("Encrypted message:", encrypted)


def decrypt_message_interactive(priv = None):
    """
    Decrypt a message
    """

    encrypted = input('Please enter the message to be decrypted: ')
    if priv is None:
        priv = enter_key_interactive('private')
    message = ''
    for i in range(0, len(encrypted), HEX_LENGTH):
        enc_string = encrypted[i:i + HEX_LENGTH]
        enc = int(enc_string, 16)
        dec = apply_key(priv, enc)
        if dec >= 0 and dec < 256:
            message += chr(dec)
        else:
            print('Warning: Could not decode encrypted entity: ' + enc_string)
            print('         decrypted as: ' + str(dec) + ' which is out of range.')
            print('         inserting _ at position of this character')
            message += '_'
    print("Decrypted message:", message)


def break_key_interactive():
    """
    Break key, interactively
    """

    pub = enter_public_key_interactive()
    priv = break_key(pub)
    print("Private key:")
    print(priv)
    decrypt_message_interactive(priv)


def enter_public_key_interactive():
    """
    Prompt user to enter the public modulus.

    :return: the tuple (e,n)
    """

    print('(Using public exponent = ' + str(PUBLIC_EXPONENT) + ')')
    string_modulus = input('Please enter the modulus (decimal): ')
    modulus = int(string_modulus)
    return (PUBLIC_EXPONENT, modulus)


def enter_key_interactive(key_type):
    """
    Prompt user to enter the exponent and modulus of a key

    :param key_type: either the string 'public' or 'private' -- used to prompt the user on how
                     this key is interpretted by the program.
    :return: the tuple (e,n)
    """
    string_exponent = input('Please enter the ' + key_type + ' exponent (decimal): ')
    exponent = int(string_exponent)
    string_modulus = input('Please enter the modulus (decimal): ')
    modulus = int(string_modulus)
    return (exponent, modulus)


def compute_checksum(string):
    """
    Compute simple hash

    Given a string, compute a simple hash as the sum of characters
    in the string.

    (If the sum goes over sixteen bits, the numbers should "wrap around"
    back into a sixteen bit number.  e.g. 0x3E6A7 should "wrap around" to
    0xE6A7)

    This checksum is similar to the internet checksum used in UDP and TCP
    packets, but it is a two's complement sum rather than a one's
    complement sum.

    :param str string: The string to hash
    :return: the checksum as an integer
    """

    total = 0
    for c in string:
        total += ord(c)
    total %= 0x8000  # Guarantees checksum is only 4 hex digits
    # How many bytes is that?
    #
    # Also guarantees that that the checksum will
    # always be less than the modulus.
    return total


# ---------------------------------------
# Do not modify code above this line
# ---------------------------------------


def create_keys():
    """
    Create the public and private keys.
    :author: Garin Jankowski
    :return: the keys as a three-tuple: (e,d,n)
    """
    values = verify_primes(generate_prime(), generate_prime())
    prime1 = values[0]
    prime2 = values[1]
    modulus = values[2]
    totient = (prime1-1)*(prime2-1)

    public_key = PUBLIC_EXPONENT
    private_key = extended_euclidean_algorithm(public_key, totient)
    return public_key, private_key, modulus


def generate_prime():
    """
    generate a prime number that is within the global range
    :return: prime number in between MIN_PRIME and MAX_PRIME
    :author: Garin Jankowski
    """
    random_num = random.randrange(MIN_PRIME, MAX_PRIME)
    while not isprime(random_num):
        random_num += 2
        if random_num > MAX_PRIME:
            random_num = random.randrange(MIN_PRIME, MAX_PRIME)
    return random_num


def verify_primes(prime1, prime2):
    """
    verifies two things about each prime:
        (prime-1) % modulus = 0
        (prime-1) and the PUBLIC_EXPONENT are coprime
    :param int prime1: first prime
    :param int prime2: second prime
    :return: a valid first prime, valid second prime, and their product
    :author: Garin Jankowski
    """
    modulus = prime1*prime2
    while ((prime1 - 1) % modulus) == 0 or math.gcd(prime1-1, PUBLIC_EXPONENT) > 1:
        prime1 = generate_prime()
        modulus = prime1*prime2
    while ((prime2 - 1) % modulus) == 0 or math.gcd(prime2 - 1, PUBLIC_EXPONENT) > 1 or prime2 == prime1:
        prime2 = generate_prime()
        modulus = prime1 * prime2
    return prime1, prime2, modulus


def isprime(n: int) -> bool:
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


def extended_euclidean_algorithm(firstkey, totient):
    """
    :return: modular multiplicative inverse of firstkey in mod(totient)
    :author: William Retert
    """
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


def apply_key(key, m):
    """
    Apply the key, given as a tuple (e,n) or (d,n) to the message.

    This can be used both for encryption and decryption.

    :author: Garin Jankowski
    :param tuple key: (e,n) or (d,n)
    :param int m: the message as a number 1 < m < n (roughly)
    :return: the message with the key applied. For example,
             if given the public key and a message, encrypts the message
             and returns the ciphertext.
    """
    return m**key[0] % key[1]


def break_key(pub):
    """
    Break a key.  Given the public key, find the private key.
    Factorizes the modulus n to find the prime numbers p and q.

    You can follow the steps in the "optional" part of the in-class
    exercise.

    :author: Vishnu Appalaraju
    :param pub: a tuple containing the public key (e,n)
    :return: a tuple containing the private key (d,n)
    """
    # with n you find p and q
    # totient (z) = (p-1) x (q-1)
    # to find d, you do euclid(e, z)
    factors = get_factors(pub[1])
    totient = calculate_totient(factors[0], factors[1])
    d = extended_euclidean_algorithm(pub[0], totient)
    return int(d), pub[1]


def calculate_totient(p, q):
    """
    Calculates the totient.
    :author: Vishnu Appalaraju
    :param: the factors
    :return: the totient
    """
    return (p-1) * (q-1)


def get_factors(e):
    """
    Gets the factors.
    :author: Vishnu Appalaraju
    :param: the totient
    :return: the factors as a tuple
    """
    p = 0
    q = 0
    sqrt_n = math.sqrt(e)
    i = 0
    while p == 0 and q == 0 and i < sqrt_n:
        if isprime(i) and e%i == 0:
            p = i
            q = e/i
        i += 1
    return p, q


# Your code and additional functions go here. (Replace this line.)

# ---------------------------------------
# Do not modify code below this line
# ---------------------------------------


def get_public_key(key_pair):
    """
    Pulls the public key out of the tuple structure created by
    create_keys()

    :param key_pair: (e,d,n)
    :return: (e,n)
    """

    return (key_pair[0], key_pair[2])


def get_private_key(key_pair):
    """
    Pulls the private key out of the tuple structure created by
    create_keys()

    :param key_pair: (e,d,n)
    :return: (d,n)
    """

    return (key_pair[1], key_pair[2])


def as_hex(number):
    """
    Convert integer to a zero-padded hex string with the required number
    of characters to represent n, d, or and encrypted message.

    :param int number: to format
    :return: The formatted string
    """

    return "{:0{digits}x}".format(number,digits=str(HEX_LENGTH))

main()