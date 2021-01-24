import QSHelper as QS
import math
import copy


class QuadraticSieve:
    n = 0

    sieving_interval = 0
    base_size = 0
    factor_base = list()
    residue_roots = list()

    z_list = list()
    z_matrix = list()
    matrix_solutions = list()

    def factor(self, n):
        self.n = n
        print(self.n)
        print("\tBits: ", self.n.bit_length(), ", Digits: ", int(math.log10(self.n)+1), sep="")
        self.__setup_values()
        print("\tBase size: ", self.base_size, "\n\tM: ", self.sieving_interval, sep="")
        print("Generating relations...")
        self.__run_polynomials()
        print("Solving matrix...")
        self.__find_matrix_solutions()
        factors = self.__run_matrix_solutions()
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

    def __generate_new_a_factors(self, q_min):
        q = QS.generate_quadratic_residue_primes(1, self.n, q_min)[0]
        a_factors = {q: 2}
        return a_factors

    def __run_polynomials(self):
        q_min = self.factor_base[-1]+1
        num_relations_goal = self.base_size+10

        num_polynomials = 0

        while len(self.z_list) < num_relations_goal:
            a_factors = self.__generate_new_a_factors(q_min)
            polynomial = self.QSPolynomial(self.n, a_factors, self.sieving_interval, self.factor_base, self.residue_roots, self.z_list, self.z_matrix)
            new_q_min = polynomial.run_polynomial()

            q_min = new_q_min + 1
            num_polynomials += 1
        print("\tPolynomials:", num_polynomials)

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
                print("HEY:", square_a, square_b)

            if (root_a % self.n) == (root_b % self.n) or (root_a % self.n) == ((-root_b) % self.n):
                continue

            factor1 = QS.GCD_Stein_while(abs(root_a - root_b), self.n)
            if factor1 != 1 and factor1 != self.n:
                return factor1, self.n // factor1
            else:
                # EXCEPTION HERE
                pass
        return None

    class QSPolynomial:
        n = 0

        sieving_interval = 0
        complete_factor_base = list()
        poly_factor_base = list()
        residue_roots = list()

        a_factors = dict()
        a = 0
        b = 0
        c = 0

        y_list = dict()
        y_copy = dict()
        y_prime_factors = dict()

        factor_indices_p = list()
        factor_indices_n = list()

        partial_relations = dict()
        partial_prime_factors = dict()
        partial_index = 0

        z_list = list()
        z_matrix = list()

        def __init__(self, n, a_factors, sieving_interval, factor_base, residue_roots, z_list, z_matrix):
            self.n = n

            self.a_factors = a_factors
            self.sieving_interval = sieving_interval
            self.complete_factor_base = factor_base
            self.residue_roots = residue_roots

            self.a = 0
            self.b = 0
            self.c = 0

            self.y_list = dict()
            self.y_copy = dict()
            self.y_prime_factors = dict()

            self.factor_indices_p = list()
            self.factor_indices_n = list()

            self.partial_relations = dict()
            self.partial_prime_factors = dict()
            self.partial_index = 0

            self.z_list = z_list
            self.z_matrix = z_matrix

        def run_polynomial(self):
            self.__setup_coefficients_and_FB()
            self.__calculate_indices()
            self.__generate_y_list()
            self.__sieve()

            return next(iter(self.a_factors.keys()))

        def __setup_coefficients_and_FB(self):
            self.a = QS.prime_factors_to_int(self.a_factors)
            q = next(iter(self.a_factors.keys()))
            self.b = QS.Tonelli_Shanks(self.n, q, 2)[0]
            self.c = (self.b ** 2 - self.n) // self.a

            # gonna need to change this eventually
            self.poly_factor_base = self.complete_factor_base
            self.residue_roots = self.residue_roots

        def __calculate_indices(self):
            self.factor_indices_p = copy.deepcopy(self.residue_roots)
            self.factor_indices_n = copy.deepcopy(self.residue_roots)

            for i in range(len(self.poly_factor_base)):
                p = self.poly_factor_base[i]
                a_inverse = QS.extended_euclidean_algorithm(self.a, p)
                for k in range(len(self.factor_indices_p[i])):
                    t = self.factor_indices_p[i][k]
                    self.factor_indices_p[i][k] = (self.b + t) * -a_inverse % p
                    self.factor_indices_n[i][k] = self.factor_indices_p[i][k] - p

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
            factor_base_indices = [i for i in range(len(self.poly_factor_base))]

            while len(factor_base_indices) > 0:
                p = 0
                # go through the current indices
                while p < len(factor_base_indices):
                    current_indices = self.factor_indices_p
                    negative = 1
                    f_index = factor_base_indices[p]
                    # do negative and positive direction
                    for j in range(2):
                        i = 0
                        # do division up to the current index, increment current index by prime
                        while i < len(current_indices[f_index]):
                            if abs(current_indices[f_index][i]) >= self.sieving_interval:
                                i += 1
                                continue

                            factor_i = current_indices[f_index][i]
                            self.__division_check(factor_i, self.poly_factor_base[f_index])

                            if self.y_copy[factor_i] == 1:
                                self.__add_to_z_list(factor_i)

                            current_indices[f_index][i] += self.poly_factor_base[f_index] * negative
                            i += 1
                        current_indices = self.factor_indices_n
                        negative = -1
                    cross_out_prime = True
                    for f in range(len(self.factor_indices_p[f_index])):
                        if self.factor_indices_p[f_index][f] < self.sieving_interval or abs(
                                self.factor_indices_n[f_index][f]) < self.sieving_interval:
                            cross_out_prime = False
                    if cross_out_prime:
                        factor_base_indices.pop(p)
                    else:
                        p += 1
                self.__large_prime_procedure()

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
                # print(factor_base[f_index], factor_i, y_list[factor_i], y_list_copy[factor_i], n)
                print("WRONG WRONG WRONG")
                pass

        def __large_prime_procedure(self):
            while self.partial_index < self.factor_indices_p[0][0] and self.partial_index < self.sieving_interval:
                negative = 1
                for i in range(2):
                    factor_i = self.partial_index * negative
                    if self.y_copy[factor_i] != 1 and self.complete_factor_base[-1] < self.y_copy[factor_i] <= self.complete_factor_base[-1] ** 2:
                        large_prime = self.y_copy[factor_i]
                        if large_prime not in self.partial_relations:
                            partial_y = self.y_list[factor_i]

                            self.partial_relations[large_prime] = factor_i, partial_y
                            self.partial_prime_factors[partial_y] = self.y_prime_factors[partial_y]
                        else:
                            partial_x, partial_y = self.partial_relations[large_prime]

                            new_partial_x = factor_i
                            new_partial_y = self.y_list[factor_i]

                            new_x = ((self.a * partial_x + self.b) ** 2) * ((self.a * new_partial_x + self.b) ** 2)
                            new_y = partial_y * new_partial_y * self.a**2
                            new_factors = QS.multiply_prime_factors(self.partial_prime_factors[partial_y],
                                                                    self.y_prime_factors[self.y_list[factor_i]])
                            new_factors = QS.multiply_prime_factors(new_factors, self.a_factors)
                            new_factors = QS.multiply_prime_factors(new_factors, self.a_factors)

                            tup = new_x, new_y

                            self.partial_relations.pop(large_prime)

                            if tup not in self.z_list:
                                self.z_list.append(tup)
                                self.z_matrix.append(self.prime_factors_to_parity(new_factors, self.complete_factor_base))
                    negative = -1
                self.partial_index += 1

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
                self.z_matrix.append(self.prime_factors_to_parity(prime_factors, self.complete_factor_base))

        @staticmethod
        def prime_factors_to_parity(factor_dict, primes):
            vector = []
            primes = [-1] + primes
            for prime in primes:
                if prime in factor_dict:
                    vector.append(factor_dict[prime] % 2)
                else:
                    vector.append(0)
            return vector
