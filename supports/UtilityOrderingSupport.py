import numpy as np
from scipy.linalg import cholesky
from statistics import mean


def CorrelateUtils(matrix_u, matrix_q, epsilon, delta):
    """
    reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4794424/#APP1
    :param matrix_u: (numpy array) independent QoL weight matrix (nxs). n: row, PSA; s: col, health states
    :param matrix_q: (numpy array) preference order matrix (sxs). s: # of health states
    :param epsilon: % violation threshold
    :param delta: step for correlation coefficient reduction (rho)
    """
    # Step 1: initializations
    n = matrix_u.shape[0]  # number of PSA samples
    s = matrix_u.shape[1]  # number of health states

    # reference matrix
    matrix_r = np.random.normal(loc=0, scale=1, size=(n, s))
    # desired correlation matrix
    matrix_c = np.zeros(shape=(s, s))
    # violation matrix
    Viol = np.zeros(shape=(s, s))

    # Step 2: find minimum correlation coefficients for pairwise health states
    # [j, k] is the selected pair for state comparison
    for j in range(1, s):
        for k in range(0, j):
            rho = 1  # bi-variate correlation
            matrix_x = matrix_u[:, [j, k]]      # select column of U (independent QoL weight matrix)
            matrix_y = matrix_r[:, [j, k]]      # select column of R (reference matrix)
            viol = 0
            while (viol < epsilon) and (rho >= 0):  # if conditions met, continues
                rho = rho - delta  # reduce correlation
                matrix_x_star = induceRankCorrelation(matrix_x, matrix_y, rho)
                # compute correlations between the col vectors (calculate the % of violation)
                viol_list = (matrix_q[j, k] * matrix_x_star[:, 0]) < (matrix_q[j, k] * matrix_x_star[:, 1])
                viol_list = viol_list.tolist()
                viol = mean(viol_list)
            matrix_c[j, k] = rho + delta    # record the desired correlation into desired correlation matrix

    matrix_c += matrix_c.T
    for j in range(1, s):
        matrix_c[j, j] = 1      # fill the diagonal ones

    # Step 3: Eigenvector and Eigenvalues correlation of C
    B_orig, V = np.linalg.eig(matrix_c)
    # set eigenvalues<=0 to a very small positive number
    B = []
    for value in B_orig:
        if value <= 0:
            value = 0.0001
        B.append(value)
    # reconstruct C (C* = V \times B \times V-1)
    matrix_c_star = np.dot(np.dot(V, np.diag(B)), np.linalg.inv(V))
    # similar to above, induce the correlation
    matrix_u_star = induceRankCorrelation(matrix_u, matrix_r, matrix_c_star)

    return matrix_u_star


def induceRankCorrelation(matrix_x, matrix_y, sigma):
    # if Sigma is a single value, convert to 2x2 matrix
    if isinstance(sigma, float):
        sigma = np.asarray([[1, sigma], [sigma, 1]])

    n = matrix_x.shape[0]   # num of rows
    s = matrix_x.shape[1]   # num of columns

    # Initial matrices (not needed?)
    matrix_x_sorted = np.zeros(shape=(n, s))
    matrix_y_rank = np.zeros(shape=(n, s))
    matrix_x_star = np.zeros(shape=(n, s))

    # compute the upper triangular matrix: Sigma = PP'
    matrix_p = cholesky(sigma)  # np.linalg.cholesky gives lower triangular matrix
    # sort the values in the reference factors by multiplying by matrix_p
    matrix_y_star = np.dot(matrix_y, matrix_p)
    # cor(Ystar) in R
    # print('cor() function in R:')
    # print(np.corrcoef(matrix_y_star))

    # sort X, then reorder X_sorted based on the rank order in Y*
    matrix_x_sorted = np.sort(matrix_x, axis=0)
    matrix_y_rank = np.argsort(matrix_y_star, axis=0)
    matrix_x_star = np.array(list(map(lambda x, y: y[x], matrix_y_rank.T, matrix_x_sorted.T)))

    return matrix_x_star.T


def calculate_violation_percentage(matrix):
    i = 0
    n_row = matrix.shape[0]
    for n in range(0, n_row):
        if matrix[n, 0] >= matrix[n, 1] >= matrix[n, 2]:
            i += 1

    print(1 - i / n_row)
