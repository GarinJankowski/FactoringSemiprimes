import QSHelper as QS
import math
import copy


class QuadraticSieve:
    def __init__(self):
        self.n = 0

        self.sieving_interval = 0
        self.base_size = 0
        self.factor_base = list()
        self.residue_roots = list()

        self.z_list = list()
        self.z_matrix = list()
        self.matrix_solutions = list()

    def factor(self, n, print_things=True):
        self.n = n

        if print_things:
            print(self.n)
            print("\tBits: ", self.n.bit_length(), ", Digits: ", int(math.log10(self.n)+1), sep="")
        self.__setup_values()

        if print_things:
            print("\tBase size: ", self.base_size, "\n\tM: ", self.sieving_interval, sep="")
            print("Generating relations...")
        num_polys = self.__run_polynomials()

        if print_things:
            print("\tPolynomials:", num_polys)
            print("Solving matrix...")
        self.__find_matrix_solutions()
        factors = self.__run_matrix_solutions()

        if print_things:
            print(factors)
        if factors[0]*factors[1] != self.n:
            raise SievingError("Found factors are INCORRECT.")
        return factors

    def __setup_values(self):
        num_digits = int(math.log10(self.n) + 1)
        self.sieving_interval = int((num_digits ** 4) / 50)

        self.base_size = int((math.e ** (math.sqrt(math.log(self.n) * math.log(math.log(self.n))))) ** (math.sqrt(2) / 4))
        self.factor_base = QS.generate_quadratic_residue_primes(self.base_size, self.n)
        self.residue_roots = [sorted(QS.Tonelli_Shanks(self.n, p)) for p in self.factor_base]

        self.z_list = list()
        self.z_matrix = list()
        self.matrix_solutions = list()

    def __run_polynomials(self):
        num_relations_goal = self.base_size+12
        coefficient_gen = self.CoefficientGenerator(self)
        num_polynomials = 0

        while len(self.z_list) < num_relations_goal:
            self.__run_poly(coefficient_gen)
            num_polynomials += 1
        return num_polynomials

    def __run_poly(self, coefficient_gen):
        a_factors, a, b, c, poly_factor_base, poly_residue_roots = coefficient_gen.get_next_poly_values()
        self.Polynomial(self, a_factors, a, b, c, poly_factor_base, poly_residue_roots)

    def __find_matrix_solutions(self):
        left_null_space = QS.left_null_space(copy.deepcopy(self.z_matrix), 2)
        parity_matrix = QS.mod_matrix(left_null_space[0], 2)
        new_identity_matrix = QS.mod_matrix(left_null_space[1], 2)

        for r in range(len(parity_matrix)):
            if sum(parity_matrix[r]) == 0:
                self.matrix_solutions.append(new_identity_matrix[r])

    def __run_matrix_solutions(self):
        for solution in self.matrix_solutions:
            square_a = 1
            square_b = 1

            for i in range(len(solution)):
                if solution[i] == 1:
                    square_a *= self.z_list[i][0]
                    square_b *= self.z_list[i][1]
            root_a = math.isqrt(square_a)
            root_b = math.isqrt(square_b)

            if square_a % self.n != square_b % self.n:
                # EXCEPTION HERE
                raise SievingError("HEY. Calculated congruence of squares are NOT CONGRUENT.")

            if (root_a % self.n) == (root_b % self.n) or (root_a % self.n) == ((-root_b) % self.n):
                continue

            factor1 = QS.GCD_Stein_while(abs(root_a - root_b), self.n)
            if factor1 != 1 and factor1 != self.n:
                return factor1, self.n // factor1
            else:
                # EXCEPTION HERE
                raise SievingError("UH OH. Results from congruence are not factors of n. The math is messed up somewhere.")
        return None

    class CoefficientGenerator:
        def __init__(self, quadratic_sieve):
            self.n = quadratic_sieve.n

            self.sieving_interval = quadratic_sieve.sieving_interval
            self.factor_base = quadratic_sieve.factor_base
            self.factor_base_dict = {self.factor_base[i]: i for i in range(len(self.factor_base))}
            self.residue_roots = quadratic_sieve.residue_roots

            self.q_min = math.sqrt(math.sqrt(2*self.n)/self.sieving_interval)

        def get_next_poly_values(self):
            a_factors, a, b, c = self.__generate_coefficients()
            poly_fb, poly_rr = self.__get_poly_factor_base(a_factors)
            return a_factors, a, b, c, poly_fb, poly_rr

        def __generate_coefficients(self):
            q = QS.generate_quadratic_residue_primes(1, self.n, self.q_min)[0]
            a_factors = {q: 2}

            a = QS.prime_factors_to_int(a_factors)
            b = QS.Tonelli_Shanks(self.n, q, 2)[0]
            c = (b ** 2 - self.n) // a

            self.q_min = q + 1

            return a_factors, a, b, c

        def __get_poly_factor_base(self, a_factors):
            index_list = []
            for key in a_factors:
                if key in self.factor_base_dict:
                    index_list.append(self.factor_base_dict[key])

            index_list = sorted(index_list)
            index_offset = 0
            poly_factor_base = list(self.factor_base)
            poly_residue_roots = list(self.residue_roots)
            for i in index_list:
                i = i - index_offset
                poly_factor_base.pop(i)
                poly_residue_roots.pop(i)
                index_offset += 1
            return poly_factor_base, poly_residue_roots

    class Polynomial:
        def __init__(self, quadratic_sieve, a_factors, a, b, c, poly_factor_base, poly_residue_roots):
            self.n = quadratic_sieve.n
            self.sieving_interval = quadratic_sieve.sieving_interval
            self.complete_factor_base = quadratic_sieve.factor_base
            self.z_list = quadratic_sieve.z_list
            self.z_matrix = quadratic_sieve.z_matrix

            self.a_factors = a_factors
            self.poly_factor_base = poly_factor_base
            self.poly_residue_roots = poly_residue_roots

            self.a = a
            self.b = b
            self.c = c

            self.y_list = dict()
            self.y_copy = dict()
            self.y_prime_factors = dict()

            self.factor_indices = list()

            self.partial_relations = dict()
            self.partial_prime_factors = dict()
            self.partial_index = 0

            self.__run_polynomial()

        def __run_polynomial(self):
            self.__calculate_indices()
            self.__generate_y_list()
            self.__sieve()

        def __calculate_indices(self):
            self.factor_indices = copy.deepcopy(self.poly_residue_roots)

            for i in range(len(self.poly_factor_base)):
                p = self.poly_factor_base[i]
                a_inverse = QS.extended_euclidean_algorithm(self.a, p)
                for k in range(len(self.factor_indices[i])):
                    t = self.factor_indices[i][k]
                    self.factor_indices[i][k] = (((self.b + t) * (-a_inverse)) % p) - ((self.sieving_interval//p)*p)

        def __Q(self, x):
            return int(self.a*x**2 + 2*self.b*x + self.c)

        def __generate_y_list(self):
            for x in range(self.sieving_interval+1):
                self.y_list[x] = self.__Q(x)
                self.y_list[-x] = self.__Q(-x)
            self.y_copy = copy.deepcopy(self.y_list)

            for x in self.y_list:
                self.y_prime_factors[self.y_list[x]] = {}

        def __sieve(self):
            current_indices = self.factor_indices
            for p in range(len(current_indices)):
                prime = self.poly_factor_base[p]
                for x in current_indices[p]:
                    while x < self.sieving_interval:
                        self.__division_check(x, prime)
                        x += prime
            for x in self.y_list:
                if self.y_copy[x] == 1:
                    self.__add_to_z_list(x)
                else:
                    self.__large_prime_check(x)

        def __division_check(self, x, p):
            true_y = self.y_list[x]
            if self.y_copy[x] < 0:
                self.y_copy[x] *= -1
                self.y_prime_factors[true_y][-1] = 1
            if self.y_copy[x] % p == 0:
                self.y_prime_factors[true_y][p] = 0
                while self.y_copy[x] % p == 0:
                    self.y_copy[x] //= p
                    self.y_prime_factors[true_y][p] += 1
            else:
                # EXCEPTION HERE
                raise SievingError("WRONG. Sieving indices are messed up. Q(x) not factorable by prime.")

        def __large_prime_check(self, x):
            if self.y_copy[x] != 1 and self.poly_factor_base[-1] < self.y_copy[x] <= \
                    self.poly_factor_base[-1] ** 2:
                large_prime = self.y_copy[x]
                if large_prime not in self.partial_relations:
                    partial_y = self.y_list[x]

                    self.partial_relations[large_prime] = x, partial_y
                    self.partial_prime_factors[partial_y] = self.y_prime_factors[partial_y]
                else:
                    partial_x, partial_y = self.partial_relations[large_prime]
                    new_partial_x = x
                    new_partial_y = self.y_list[x]

                    new_x = ((self.a * partial_x + self.b) ** 2) * ((self.a * new_partial_x + self.b) ** 2)
                    new_y = partial_y * new_partial_y * self.a ** 2

                    new_factors = QS.multiply_prime_factors(self.partial_prime_factors[partial_y],
                                                            self.y_prime_factors[self.y_list[x]])
                    new_factors = QS.multiply_prime_factors(new_factors, self.a_factors)
                    new_factors = QS.multiply_prime_factors(new_factors, self.a_factors)

                    tup = new_x, new_y

                    self.partial_relations.pop(large_prime)

                    if tup not in self.z_list:
                        self.z_list.append(tup)
                        self.z_matrix.append(self.__prime_factors_to_parity(new_factors, self.complete_factor_base))

        # not using this right now, probably should
        def __convert_to_relation(self, x, y):
            new_x = (self.a*x + self.b)**2
            new_y = self.a*y
            new_y_factors = QS.multiply_prime_factors(self.y_prime_factors[y], self.a_factors)

            return new_x, new_y, new_y_factors

        def __add_to_z_list(self, x):
            index = (self.a*x + self.b)**2
            value = self.y_list[x]*self.a
            tup = index, value
            if tup not in self.z_list:
                self.z_list.append(tup)
                prime_factors = QS.multiply_prime_factors(self.y_prime_factors[self.y_list[x]], self.a_factors)
                self.z_matrix.append(self.__prime_factors_to_parity(prime_factors, self.complete_factor_base))

        @staticmethod
        def __prime_factors_to_parity(factor_dict, primes):
            vector = []
            primes = [-1] + primes
            for prime in primes:
                if prime in factor_dict:
                    vector.append(factor_dict[prime] % 2)
                else:
                    vector.append(0)
            return vector


class SievingError(Exception):
    def __init__(self, message):
        super().__init__(message)
