import math
import numpy as np
import scipy.linalg as linalg
import importScripts.exceptions as exceptions


def adaptiveWheeler(m, n, time):

    nMoms = 2*n

    sigma = np.zeros(n*(n + 1))

    # Assigning the bottom row of the 2D sigma matrix
    sigma[0:nMoms] = m

    # Defining array for diagonal elements of the Jacobi matrix
    a = np.zeros(n)

    # Assigning the first element of array "a"
    atmp = sigma[1]/sigma[0]
    a[0] = atmp

    # Assigning the second row of the sigma matrix
    for i in range(nMoms - 2):
        sigma[i + nMoms] = sigma[i + 2] - atmp*sigma[i + 1]

    # Defining array for off-diagonal elements of the Jacobi matrix
    b = np.zeros(n-1)

    # Assigning the first element of the array "b"
    btmp = sigma[nMoms] / sigma[0]
    b[0] = btmp

    # Defining zeta array for realizability check
    zeta = np.zeros(nMoms - 1)

    zeta[0] = atmp

    zetatmp = btmp / atmp
    zeta[1] = zetatmp

    if zetatmp > 0:
        # Assigning the second element of the array "a"
        atmp = sigma[nMoms + 1] / sigma[nMoms] - atmp
        a[1] = atmp

        zeta[2] = atmp - zetatmp

    # Loop over the next rows of the sigma matrix
    for j in range(2, n):

        if (zetatmp <= 0):
            break

        nCol = nMoms - 2*j

        # Auxiliary indexes for mapping the sigma matrix from 2D to 1D
        fb = nMoms*j + j*(1 - j)
        fb_1 = fb - nCol - 2
        fb_2 = fb_1 - nCol - 4

        # Loop over the columns of each row and assign the elements
        for i in range(nCol):
            sigma[i + fb] = sigma[i + 2 + fb_1] - atmp * sigma[i + 1 + fb_1] \
                - btmp * sigma[i + 2 + fb_2]

        atmp = sigma[fb + 1]/sigma[fb] - sigma[fb_1 + 1]/sigma[fb_1]
        btmp = sigma[fb] / sigma[fb_1]

        a[j] = atmp
        b[j-1] = btmp

        zetatmp = btmp / zeta[2*j-2]
        zeta[2*j-1] = zetatmp
        zeta[2*j] = atmp - zeta[2*j-1]

    # Reduce number of nodes (n) if the moments are unrealizable
    _, rN = nodeReduction(zeta, n, time)

    # Reduced "a" and "b" arrays (in case of unrealizable moment set)
    aR = np.zeros(rN)
    bR = np.zeros(rN-1)

    for i in range(rN-1):
        aR[i] = a[i]
        bR[i] = -1.0 * (b[i]**0.5)

    aR[rN-1] = a[rN-1]

    D, V = linalg.eigh_tridiagonal(aR, bR)

    w = np.zeros(n)
    x = np.zeros(n)

    for i in range(0, rN):
        w[i] = (V[0, i]**2)*m[0]
        if w[i] < 0.0:
            raise exceptions.RealizabilityErr(m)
        elif w[i] > 1e-15:
            x[i] = D[i]
            if x[i] < 0.0:
                raise exceptions.RealizabilityErr(m, x[i], negNode=True)
        else:
            x[i] = 1e-20

    for i in range(rN, n):
        w[i] = 0.0
        x[i] = 1e-20  # This size is not important since the weight is zero

    return w, x


def wheeler(m, n):

    sigma = np.zeros([2*n + 1, 2*n + 1])

    for i in range(1, 2*n + 1):
        sigma[1, i] = m[i-1]

    a = np.zeros(n)
    b = np.zeros(n)

    a[0] = m[1]/m[0]

    b[0] = 0

    for i in range(2, n+1):
        for j in range(i, 2*n - i + 2):
            sigma[i, j] = sigma[i - 1, j + 1] - a[i - 2]*sigma[i - 1, j] \
                - b[i - 2]*sigma[i - 2, j]

        a[i - 1] = sigma[i, i + 1]/sigma[i, i] \
            - sigma[i - 1, i]/sigma[i - 1, i - 1]
        b[i - 1] = sigma[i, i]/sigma[i - 1, i - 1]

    negSqrt_b = np.zeros(n-1)

    for i in range(0, n-1):
        b_i = b[i + 1]

        if b_i < 0:
            raise exceptions.RealizabilityErr(m)

        negSqrt_b[i] = -math.sqrt(b_i)

    D, V = linalg.eigh_tridiagonal(a, negSqrt_b)

    w = np.zeros(n)
    x = np.zeros(n)

    for i in range(0, n):
        w[i] = (V[0, i]**2)*m[0]
        if w[i] < 0:
            raise exceptions.RealizabilityErr(m)

        x[i] = D[i]
        if x[i] < 0:
            raise exceptions.RealizabilityErr(m, x[i], negNode=True)

    return w, x


def PD(m, n):

    fIndex = 2*n + 1

    P = np.zeros(int(fIndex*(fIndex + 1)/2))

    P[0] = 1.0

    for i in range(fIndex - 1):
        P[i + fIndex] = ((-1.0)**i) * m[i]

    zeta = np.zeros(fIndex - 1)

    for j in range(2, fIndex):
        fb = int(j*(2*fIndex + 1 - j) / 2)
        fb_1 = fb - fIndex + j - 1
        fb_2 = fb_1 - fIndex + j - 2

        for i in range(fIndex - j):
            P[i + fb] = P[fb_1] * P[i + 1 + fb_2] - P[fb_2] * P[i + 1 + fb_1]

            product = P[fb_1] * P[fb_2]

            if (product > 0):
                zeta[j - 1] = P[fb] / product

    a = np.zeros(n)
    b = np.zeros(n-1)

    for i in range(n - 1):
        a[i] = zeta[2*i + 1] + zeta[2*i]
        b[i] = -1.0 * np.sqrt(zeta[2*i + 2] * zeta[2*i + 1])

    a[n-1] = zeta[fIndex - 2] + zeta[fIndex - 3]


    D, V = linalg.eigh_tridiagonal(a, b)

    w = np.zeros(n)
    x = np.zeros(n)

    for i in range(0, n):
        w[i] = (V[0, i]**2)*m[0]
        if w[i] < 0:
            raise exceptions.RealizabilityErr(m)

        x[i] = D[i]
        if x[i] < 0:
            raise exceptions.RealizabilityErr(m, x[i], negNode=True)

    return w, x


def nodeReduction(zeta, n, time):

    isRealizable = True

    rN = 1

    if (zeta[0] <= 0.0):
        # print("Warning:\nNumber of nodes is reduced to " +
        #       "{0} at time {1}\n".format(rN, time))
        return False, rN

    for i in range(2, 2*n-1, 2):
        if (zeta[i] <= 0 or zeta[i-1] <= 0):
            isRealizable = False
            # print("Warning:\nNumber of nodes is reduced to " +
            #       "{0} at time {1}\n".format(rN, time))
            break
        rN += 1

    return isRealizable, rN
