

def add_rows(ra, rb):
    for i in range(len(ra)):
        rb[i] += ra[i]
    return rb


def multiply_row(r, a):
    for i in range(len(r)):
        r[i] *= a
    return r


def swap_rows(matrix, a, b):
    temp_a = list(matrix[a])
    matrix[a] = list(matrix[b])
    matrix[b] = temp_a
    return matrix


def add_columns(matrix, c1, c2):
    for row in matrix:
        row[c2] += row[c1]
    return matrix


def multiply_column(matrix, c, a):
    for row in matrix:
        row[c] *= a
    return matrix


def swap_columns(matrix, c1, c2):
    for row in matrix:
        temp_c = row[c1]
        row[c1] = row[c2]
        row[c2] = temp_c
    return matrix


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


def Gauss_Jordan_elimination(matrix):
    """
    puts a matrix in reduced row echelon form
    beware of numerical instability, idk what to do about it
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
                matrix[r_pivot][c2] /= pivot
            swap_rows(matrix, r, r_pivot)
            # set everything in the column to zero except for the pivot
            # also subtract the pivot row from each row that gets modified
            for r2 in range(rows):
                if r2 != r:
                    quotient = matrix[r2][c]/matrix[r][c]
                    matrix[r2][c] = 0
                    for c2 in range(c+1, columns):
                        matrix[r2][c2] -= matrix[r][c2]*quotient
            r += 1
        c += 1
    return matrix


def sideways_Gauss_Jordan_elimination(matrix):
    """
    DOES NOT WORK CORRECTLY
    puts a matrix in reduced column echelon form
    :param matrix: matrix to be solved
    :return matrix: solved matrix
    """
    rows = len(matrix)
    columns = len(matrix[0])
    r = rows-1
    c = 0
    while r >= 0 and c < columns:
        # find the row index of the number in column c with the greatest absolute value
        # this finds the location of the 'c'th pivot
        c_pivot = c
        for i in range(c, columns):
            n = abs(matrix[r][i])
            if n >= abs(matrix[r][c_pivot]):
                c_pivot = i
        # if there is no pivot in this column, aka the whole column is zeros, then skip the column
        if matrix[r][c_pivot] != 0:
            pivot = matrix[r][c_pivot]
            # divide the pivot row by the pivot to make it 1
            for r2 in range(rows-1, -1, -1):
                matrix[r2][c_pivot] /= pivot
            swap_columns(matrix, c, c_pivot)
            # set everything in the column to zero except for the pivot
            # also subtract the pivot row from each row that gets modified
            for c2 in range(columns):
                if c2 != c:
                    quotient = matrix[r][c2]/matrix[r][c]
                    matrix[r][c2] = 0
                    for r2 in range(r+1, rows):
                        matrix[r2][c2] -= matrix[r2][c]*quotient
            c += 1
        r += 1
    return matrix


def left_null_space(matrix):
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

    print_matrix(matrix)
    Gauss_Jordan_elimination(matrix)
    print()
    print_matrix(matrix)
    print()

    for i in range(rows):
        whole_row = list(matrix[i])
        matrix[i] = whole_row[0:columns]
        identity[i] = whole_row[columns: len(whole_row)]

    return matrix, identity


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


def round_matrix(matrix):
    """
    rounds everything in the matrix to integers
    note that matrices can have fractional/floating point answers so maybe this isn't really useful at all
    :param matrix: matrix to be rounded
    :return matrix: rounded matrix
    """
    rows = len(matrix)
    columns = len(matrix[0])
    for r in range(rows):
        for c in range(columns):
            if matrix[r][c] < 0:
                matrix[r][c] = int(matrix[r][c] - 0.5)
            else:
                matrix[r][c] = int(matrix[r][c] + 0.5)
    return matrix


def print_matrix(matrix):
    matrix_str = ""
    for row in matrix:
        matrix_str += str(row) + "\n"
    print(matrix_str)