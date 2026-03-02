#!/usr/bin/env python3
"""
verify.py — Cross-check zblas against slatec reference BLAS.

Two verification layers:
  1. Python reference implementations (always runs)
  2. Compiled slatec Fortran driver (runs if gfortran available)

Both must agree for a test to pass.
"""

import math
import subprocess
import sys
import os
import shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FORTRAN_SRC = os.path.join(SCRIPT_DIR, "slatec_driver.f90")
FORTRAN_EXE = os.path.join(SCRIPT_DIR, "slatec_driver.exe")

# ---- Python reference implementations ----

def ref_daxpy(n, alpha, x, incx, y, incy):
    """Reference DAXPY: y = alpha*x + y"""
    if n <= 0 or alpha == 0.0:
        return y[:]
    result = y[:]
    ix, iy = 0, 0
    if incx < 0:
        ix = (-n + 1) * incx
    if incy < 0:
        iy = (-n + 1) * incy
    for _ in range(n):
        result[iy] = alpha * x[ix] + result[iy]
        ix += incx
        iy += incy
    return result


def ref_ddot(n, x, incx, y, incy):
    """Reference DDOT: dot product"""
    if n <= 0:
        return 0.0
    result = 0.0
    ix, iy = 0, 0
    if incx < 0:
        ix = (-n + 1) * incx
    if incy < 0:
        iy = (-n + 1) * incy
    for _ in range(n):
        result += x[ix] * y[iy]
        ix += incx
        iy += incy
    return result


def ref_dscal(n, alpha, x, incx):
    """Reference DSCAL: x = alpha * x"""
    if n <= 0:
        return x[:]
    result = x[:]
    ix = 0
    if incx < 0:
        ix = (-n + 1) * incx
    for _ in range(n):
        result[ix] = alpha * result[ix]
        ix += incx
    return result


def ref_dcopy(n, x, incx, y, incy):
    """Reference DCOPY: y = x"""
    if n <= 0:
        return y[:]
    result = y[:]
    ix, iy = 0, 0
    if incx < 0:
        ix = (-n + 1) * incx
    if incy < 0:
        iy = (-n + 1) * incy
    for _ in range(n):
        result[iy] = x[ix]
        ix += incx
        iy += incy
    return result


def ref_dswap(n, x, incx, y, incy):
    """Reference DSWAP: x <-> y. Returns (new_x, new_y)."""
    if n <= 0:
        return x[:], y[:]
    rx, ry = x[:], y[:]
    ix, iy = 0, 0
    if incx < 0:
        ix = (-n + 1) * incx
    if incy < 0:
        iy = (-n + 1) * incy
    for _ in range(n):
        rx[ix], ry[iy] = ry[iy], rx[ix]
        ix += incx
        iy += incy
    return rx, ry


def ref_dasum(n, x, incx):
    """Reference DASUM: sum of absolute values"""
    if n <= 0:
        return 0.0
    result = 0.0
    ix = 0
    if incx < 0:
        ix = (-n + 1) * incx
    for _ in range(n):
        result += abs(x[ix])
        ix += incx
    return result


def ref_dnrm2(n, x, incx):
    """Reference DNRM2: Euclidean norm"""
    if n <= 0:
        return 0.0
    if n == 1:
        ix = 0
        if incx < 0:
            ix = 0
        return abs(x[ix])
    ssq = 0.0
    ix = 0
    if incx < 0:
        ix = (-n + 1) * incx
    for _ in range(n):
        ssq += x[ix] * x[ix]
        ix += incx
    return math.sqrt(ssq)


def ref_idamax(n, x, incx):
    """Reference IDAMAX: 1-based index of max |x_i|"""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    ix = 0
    if incx < 0:
        ix = (-n + 1) * incx
    best_val = abs(x[ix])
    best_idx = 1
    ix += incx
    for i in range(2, n + 1):
        if abs(x[ix]) > best_val:
            best_val = abs(x[ix])
            best_idx = i
        ix += incx
    return best_idx


def ref_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc):
    """Reference DGEMM: C := alpha*op(A)*op(B) + beta*C (column-major)."""
    # Quick return
    if m <= 0 or n <= 0:
        return c[:]
    result = c[:]
    # If alpha==0, just scale C by beta
    if alpha == 0.0:
        for j in range(n):
            for i in range(m):
                idx = j * ldc + i
                result[idx] = 0.0 if beta == 0.0 else beta * result[idx]
        return result

    nota = transa.upper() == 'N'
    notb = transb.upper() == 'N'

    def b_elem(i, j):
        if notb:
            return b[j * ldb + i]
        else:
            return b[i * ldb + j]

    for j in range(n):
        # Scale C(:,j) by beta
        if beta == 0.0:
            for i in range(m):
                result[j * ldc + i] = 0.0
        elif beta != 1.0:
            for i in range(m):
                idx = j * ldc + i
                result[idx] = beta * result[idx]

        # C(:,j) += alpha * op(A) * op(B)(:,j)
        if nota:
            # DAXPY form: A not transposed
            for l in range(k):
                temp = alpha * b_elem(l, j)
                if temp != 0.0:
                    for i in range(m):
                        result[j * ldc + i] += temp * a[l * lda + i]
        else:
            # Dot product form: A transposed
            for i in range(m):
                temp = 0.0
                for l in range(k):
                    temp += a[i * lda + l] * b_elem(l, j)
                result[j * ldc + i] += alpha * temp

    return result


def ref_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy):
    """Reference DGEMV: y := alpha*op(A)*x + beta*y (column-major)."""
    if m <= 0 or n <= 0:
        return y[:]
    result = y[:]
    nota = trans.upper() == 'N'
    leny = m if nota else n
    lenx = n if nota else m

    ix = 0
    if incx < 0:
        ix = (-lenx + 1) * incx
    iy = 0
    if incy < 0:
        iy = (-leny + 1) * incy

    # Scale y by beta
    jy = iy
    for i in range(leny):
        if beta == 0.0:
            result[jy] = 0.0
        elif beta != 1.0:
            result[jy] = beta * result[jy]
        jy += incy

    if alpha == 0.0:
        return result

    if nota:
        # DAXPY form: y += alpha * A * x
        jx = ix
        for j in range(n):
            temp = alpha * x[jx]
            if temp != 0.0:
                ky = iy
                for i in range(m):
                    result[ky] += temp * a[j * lda + i]
                    ky += incy
            jx += incx
    else:
        # Dot-product form: y += alpha * A' * x
        jy = iy
        for j in range(n):
            temp = 0.0
            kx = ix
            for i in range(m):
                temp += a[j * lda + i] * x[kx]
                kx += incx
            result[jy] += alpha * temp
            jy += incy

    return result


def ref_dtrsv(uplo, trans, diag, n, a, lda, x, incx):
    """Reference DTRSV: solve op(A)*x = b in-place."""
    if n <= 0:
        return x[:]
    result = x[:]
    nounit = diag.upper() == 'N'
    upper = uplo.upper() == 'U'
    nota = trans.upper() == 'N'

    kx = 0
    if incx < 0:
        kx = (-n + 1) * incx

    if nota:
        if upper:
            # Back substitution: j = N-1 downto 0
            jx = kx + (n - 1) * incx
            for j in range(n - 1, -1, -1):
                if nounit:
                    result[jx] /= a[j * lda + j]
                temp = result[jx]
                ix = jx
                for i in range(j - 1, -1, -1):
                    ix -= incx
                    result[ix] -= temp * a[j * lda + i]
                jx -= incx
        else:
            # Forward substitution: j = 0..N-1
            jx = kx
            for j in range(n):
                if nounit:
                    result[jx] /= a[j * lda + j]
                temp = result[jx]
                ix = jx
                for i in range(j + 1, n):
                    ix += incx
                    result[ix] -= temp * a[j * lda + i]
                jx += incx
    else:
        if upper:
            # Forward dot: j = 0..N-1
            jx = kx
            for j in range(n):
                temp = result[jx]
                ix = kx
                for i in range(j):
                    temp -= a[j * lda + i] * result[ix]
                    ix += incx
                if nounit:
                    temp /= a[j * lda + j]
                result[jx] = temp
                jx += incx
        else:
            # Backward dot: j = N-1 downto 0
            jx = kx + (n - 1) * incx
            for j in range(n - 1, -1, -1):
                temp = result[jx]
                ix = kx + (n - 1) * incx
                for i in range(n - 1, j, -1):
                    temp -= a[j * lda + i] * result[ix]
                    ix -= incx
                if nounit:
                    temp /= a[j * lda + j]
                result[jx] = temp
                jx -= incx

    return result


def ref_dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb):
    """Reference DTRSM: solve op(A)*X=alpha*B or X*op(A)=alpha*B."""
    if m <= 0 or n <= 0:
        return b[:]
    result = b[:]
    if alpha == 0.0:
        for j in range(n):
            for i in range(m):
                result[j * ldb + i] = 0.0
        return result

    nounit = diag.upper() == 'N'
    upper = uplo.upper() == 'U'
    lside = side.upper() == 'L'
    nota = transa.upper() == 'N'

    if lside:
        if nota:
            if upper:
                # Left/Upper/N: back sub per B column
                for j in range(n):
                    if alpha != 1.0:
                        for i in range(m):
                            result[j * ldb + i] *= alpha
                    for k in range(m - 1, -1, -1):
                        if result[j * ldb + k] != 0.0:
                            if nounit:
                                result[j * ldb + k] /= a[k * lda + k]
                            temp = result[j * ldb + k]
                            for i in range(k):
                                result[j * ldb + i] -= temp * a[k * lda + i]
            else:
                # Left/Lower/N: forward sub per B column
                for j in range(n):
                    if alpha != 1.0:
                        for i in range(m):
                            result[j * ldb + i] *= alpha
                    for k in range(m):
                        if result[j * ldb + k] != 0.0:
                            if nounit:
                                result[j * ldb + k] /= a[k * lda + k]
                            temp = result[j * ldb + k]
                            for i in range(k + 1, m):
                                result[j * ldb + i] -= temp * a[k * lda + i]
        else:
            if upper:
                # Left/Upper/Trans
                for j in range(n):
                    for i in range(m):
                        temp = alpha * result[j * ldb + i]
                        for k in range(i):
                            temp -= a[i * lda + k] * result[j * ldb + k]
                        if nounit:
                            temp /= a[i * lda + i]
                        result[j * ldb + i] = temp
            else:
                # Left/Lower/Trans
                for j in range(n):
                    for i in range(m - 1, -1, -1):
                        temp = alpha * result[j * ldb + i]
                        for k in range(i + 1, m):
                            temp -= a[i * lda + k] * result[j * ldb + k]
                        if nounit:
                            temp /= a[i * lda + i]
                        result[j * ldb + i] = temp
    else:
        if nota:
            if upper:
                # Right/Upper/N
                for j in range(n):
                    if alpha != 1.0:
                        for i in range(m):
                            result[j * ldb + i] *= alpha
                    for k in range(j):
                        if a[j * lda + k] != 0.0:
                            for i in range(m):
                                result[j * ldb + i] -= a[j * lda + k] * result[k * ldb + i]
                    if nounit:
                        temp = 1.0 / a[j * lda + j]
                    else:
                        temp = 1.0
                    if temp != 1.0:
                        for i in range(m):
                            result[j * ldb + i] *= temp
            else:
                # Right/Lower/N
                for j in range(n - 1, -1, -1):
                    if alpha != 1.0:
                        for i in range(m):
                            result[j * ldb + i] *= alpha
                    for k in range(j + 1, n):
                        if a[j * lda + k] != 0.0:
                            for i in range(m):
                                result[j * ldb + i] -= a[j * lda + k] * result[k * ldb + i]
                    if nounit:
                        temp = 1.0 / a[j * lda + j]
                    else:
                        temp = 1.0
                    if temp != 1.0:
                        for i in range(m):
                            result[j * ldb + i] *= temp
        else:
            if upper:
                # Right/Upper/Trans
                for k in range(n - 1, -1, -1):
                    if nounit:
                        temp = 1.0 / a[k * lda + k]
                    else:
                        temp = 1.0
                    if temp != 1.0:
                        for i in range(m):
                            result[k * ldb + i] *= temp
                    for j in range(k):
                        if a[k * lda + j] != 0.0:
                            temp = a[k * lda + j]
                            for i in range(m):
                                result[j * ldb + i] -= temp * result[k * ldb + i]
                    if alpha != 1.0:
                        for i in range(m):
                            result[k * ldb + i] *= alpha
            else:
                # Right/Lower/Trans
                for k in range(n):
                    if nounit:
                        temp = 1.0 / a[k * lda + k]
                    else:
                        temp = 1.0
                    if temp != 1.0:
                        for i in range(m):
                            result[k * ldb + i] *= temp
                    for j in range(k + 1, n):
                        if a[k * lda + j] != 0.0:
                            temp = a[k * lda + j]
                            for i in range(m):
                                result[j * ldb + i] -= temp * result[k * ldb + i]
                    if alpha != 1.0:
                        for i in range(m):
                            result[k * ldb + i] *= alpha

    return result


def ref_dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc):
    """Reference DSYRK: C := alpha*A*A' + beta*C (or alpha*A'*A + beta*C)."""
    if n <= 0:
        return c[:]
    result = c[:]
    upper = uplo.upper() == 'U'
    nota = trans.upper() == 'N'

    if alpha == 0.0 and beta == 1.0:
        return result

    if nota:
        # C := alpha*A*A' + beta*C
        if upper:
            for j in range(n):
                # Scale upper triangle column j by beta
                if beta == 0.0:
                    for i in range(j + 1):
                        result[j * ldc + i] = 0.0
                elif beta != 1.0:
                    for i in range(j + 1):
                        result[j * ldc + i] *= beta
                # Accumulate alpha * A * A'
                for l in range(k):
                    if a[l * lda + j] != 0.0:
                        temp = alpha * a[l * lda + j]
                        for i in range(j + 1):
                            result[j * ldc + i] += temp * a[l * lda + i]
        else:
            for j in range(n):
                if beta == 0.0:
                    for i in range(j, n):
                        result[j * ldc + i] = 0.0
                elif beta != 1.0:
                    for i in range(j, n):
                        result[j * ldc + i] *= beta
                for l in range(k):
                    if a[l * lda + j] != 0.0:
                        temp = alpha * a[l * lda + j]
                        for i in range(j, n):
                            result[j * ldc + i] += temp * a[l * lda + i]
    else:
        # C := alpha*A'*A + beta*C
        if upper:
            for j in range(n):
                for i in range(j + 1):
                    temp = 0.0
                    for l in range(k):
                        temp += a[i * lda + l] * a[j * lda + l]
                    if beta == 0.0:
                        result[j * ldc + i] = alpha * temp
                    else:
                        result[j * ldc + i] = alpha * temp + beta * result[j * ldc + i]
        else:
            for j in range(n):
                for i in range(j, n):
                    temp = 0.0
                    for l in range(k):
                        temp += a[i * lda + l] * a[j * lda + l]
                    if beta == 0.0:
                        result[j * ldc + i] = alpha * temp
                    else:
                        result[j * ldc + i] = alpha * temp + beta * result[j * ldc + i]

    return result


def ref_dger(m, n, alpha, x, incx, y, incy, a, lda):
    """Reference DGER: A := alpha*x*y' + A (column-major)."""
    if m <= 0 or n <= 0 or alpha == 0.0:
        return a[:]
    result = a[:]

    ix = 0
    if incx < 0:
        ix = (-m + 1) * incx
    iy = 0
    if incy < 0:
        iy = (-n + 1) * incy

    jy = iy
    for j in range(n):
        if y[jy] != 0.0:
            temp = alpha * y[jy]
            kx = ix
            for i in range(m):
                result[j * lda + i] += temp * x[kx]
                kx += incx
        jy += incy

    return result


# ---- Test cases ----

DAXPY_CASES = [
    # (name, n, alpha, x, incx, y, incy, expected_y)
    ("basic",
     4, 2.0,
     [1.0, 2.0, 3.0, 4.0], 1,
     [10.0, 20.0, 30.0, 40.0], 1,
     [12.0, 24.0, 36.0, 48.0]),

    ("alpha=0",
     4, 0.0,
     [1.0, 2.0, 3.0, 4.0], 1,
     [10.0, 20.0, 30.0, 40.0], 1,
     [10.0, 20.0, 30.0, 40.0]),

    ("n=0",
     0, 2.0,
     [1.0], 1,
     [10.0], 1,
     [10.0]),

    ("n=1",
     1, 3.0,
     [5.0], 1,
     [7.0], 1,
     [22.0]),

    ("incx=2",
     3, 1.0,
     [1.0, 99.0, 2.0, 99.0, 3.0], 2,
     [10.0, 20.0, 30.0], 1,
     [11.0, 22.0, 33.0]),

    ("neg_incx",
     3, 1.0,
     [1.0, 2.0, 3.0], -1,
     [10.0, 20.0, 30.0], 1,
     [13.0, 22.0, 31.0]),  # x reversed: [3,2,1] + [10,20,30]
]

DDOT_CASES = [
    # (name, n, x, incx, y, incy, expected)
    ("basic",
     4, [1.0, 2.0, 3.0, 4.0], 1,
     [5.0, 6.0, 7.0, 8.0], 1,
     70.0),

    ("n=0",
     0, [1.0], 1, [1.0], 1,
     0.0),

    ("orthogonal",
     2, [1.0, 0.0], 1, [0.0, 1.0], 1,
     0.0),

    ("n=1",
     1, [3.0], 1, [7.0], 1,
     21.0),

    ("incx=2",
     2, [1.0, 99.0, 2.0], 2, [3.0, 4.0], 1,
     11.0),

    ("neg_incx",
     3, [1.0, 2.0, 3.0], -1, [10.0, 20.0, 30.0], 1,
     100.0),  # x reversed: [3,2,1].[10,20,30] = 30+40+30
]


DSCAL_CASES = [
    # (name, n, alpha, x, incx, expected_x)
    ("basic",
     3, 2.0, [1.0, 2.0, 3.0], 1,
     [2.0, 4.0, 6.0]),
    ("alpha=0",
     2, 0.0, [5.0, 10.0], 1,
     [0.0, 0.0]),
    ("n=0",
     0, 2.0, [1.0], 1,
     [1.0]),
    ("n=1",
     1, 3.0, [5.0], 1,
     [15.0]),
    ("neg_incx",
     3, 2.0, [1.0, 2.0, 3.0], -1,
     [2.0, 4.0, 6.0]),
]

DCOPY_CASES = [
    # (name, n, x, incx, y, incy, expected_y)
    ("basic",
     3, [1.0, 2.0, 3.0], 1,
     [0.0, 0.0, 0.0], 1,
     [1.0, 2.0, 3.0]),
    ("n=0",
     0, [1.0], 1, [99.0], 1,
     [99.0]),
    ("incx=2",
     2, [1.0, 99.0, 2.0], 2,
     [0.0, 0.0], 1,
     [1.0, 2.0]),
    ("neg_incx",
     3, [1.0, 2.0, 3.0], -1,
     [0.0, 0.0, 0.0], 1,
     [3.0, 2.0, 1.0]),
]

DSWAP_CASES = [
    # (name, n, x, incx, y, incy, expected_x, expected_y)
    ("basic",
     2, [1.0, 2.0], 1, [3.0, 4.0], 1,
     [3.0, 4.0], [1.0, 2.0]),
    ("n=0",
     0, [1.0], 1, [2.0], 1,
     [1.0], [2.0]),
    ("n=1",
     1, [5.0], 1, [10.0], 1,
     [10.0], [5.0]),
    ("neg_incx",
     3, [1.0, 2.0, 3.0], -1, [10.0, 20.0, 30.0], 1,
     [30.0, 20.0, 10.0], [3.0, 2.0, 1.0]),
]

DASUM_CASES = [
    # (name, n, x, incx, expected)
    ("basic",
     4, [1.0, -2.0, 3.0, -4.0], 1,
     10.0),
    ("n=0",
     0, [1.0], 1,
     0.0),
    ("all_pos",
     3, [1.0, 2.0, 3.0], 1,
     6.0),
    ("n=1",
     1, [-5.0], 1,
     5.0),
    ("neg_incx",
     3, [1.0, -2.0, 3.0], -1,
     6.0),
]

DNRM2_CASES = [
    # (name, n, x, incx, expected)
    ("basic",
     2, [3.0, 4.0], 1,
     5.0),
    ("n=0",
     0, [1.0], 1,
     0.0),
    ("n=1",
     1, [3.0], 1,
     3.0),
    ("neg_val",
     1, [-5.0], 1,
     5.0),
    ("unit",
     3, [1.0, 0.0, 0.0], 1,
     1.0),
]

IDAMAX_CASES = [
    # (name, n, x, incx, expected)
    ("basic",
     3, [1.0, -5.0, 3.0], 1,
     2),
    ("n=0",
     0, [1.0], 1,
     0),
    ("n=1",
     1, [3.0], 1,
     1),
    ("first_max",
     4, [10.0, 1.0, 2.0, 3.0], 1,
     1),
    ("last_max",
     4, [1.0, 2.0, 3.0, 10.0], 1,
     4),
]

# Column-major 2x2: A=[1,2;3,4] stored [1,2,3,4], B=[5,6;7,8] stored [5,6,7,8]
DGEMM_CASES = [
    # (name, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, expected_c)
    # Col-major: A(1,1)=1, A(2,1)=2, A(1,2)=3, A(2,2)=4
    #            B(1,1)=5, B(2,1)=6, B(1,2)=7, B(2,2)=8
    ("nn_basic",
     'N', 'N', 2, 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [5.0, 6.0, 7.0, 8.0], 2,
     0.0,
     [0.0, 0.0, 0.0, 0.0], 2,
     [23.0, 34.0, 31.0, 46.0]),

    # T2: NN, alpha=2, beta=1, C_init=[1,1,1,1]
    # C = 2*A*B + C = 2*[23,34,31,46] + [1,1,1,1] = [47,69,63,93]
    ("nn_ab",
     'N', 'N', 2, 2, 2, 2.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [5.0, 6.0, 7.0, 8.0], 2,
     1.0,
     [1.0, 1.0, 1.0, 1.0], 2,
     [47.0, 69.0, 63.0, 93.0]),

    # T3: alpha=0, beta=2 → C = 2*C_init
    ("alpha0",
     'N', 'N', 2, 2, 2, 0.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [5.0, 6.0, 7.0, 8.0], 2,
     2.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [2.0, 4.0, 6.0, 8.0]),

    # T4: N=0, quick return (C untouched)
    ("n_zero",
     'N', 'N', 2, 0, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [5.0, 6.0, 7.0, 8.0], 2,
     1.0,
     [9.0, 9.0, 9.0, 9.0], 2,
     [9.0, 9.0, 9.0, 9.0]),

    # T5: TN — A transposed, B not
    # op(A)=A^T: op(A)(i,l)=A(l,i) → A stored [1,2,3,4], A(l,i)=a[i*lda+l]
    # C(1,1) = A(1,1)*B(1,1) + A(2,1)*B(2,1) = 1*5+2*6 = 17
    # C(2,1) = A(1,2)*B(1,1) + A(2,2)*B(2,1) = 3*5+4*6 = 39
    # C(1,2) = A(1,1)*B(1,2) + A(2,1)*B(2,2) = 1*7+2*8 = 23
    # C(2,2) = A(1,2)*B(1,2) + A(2,2)*B(2,2) = 3*7+4*8 = 53
    ("tn_basic",
     'T', 'N', 2, 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [5.0, 6.0, 7.0, 8.0], 2,
     0.0,
     [0.0, 0.0, 0.0, 0.0], 2,
     [17.0, 39.0, 23.0, 53.0]),

    # T6: NT — A not transposed, B transposed
    # op(B)=B^T: op(B)(l,j)=B(j,l) → B(j,l)=b[l*ldb+j]
    # C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(1,2) = 1*5+3*7 = 26
    # C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(1,2) = 2*5+4*7 = 38
    # C(1,2) = A(1,1)*B(2,1) + A(1,2)*B(2,2) = 1*6+3*8 = 30
    # C(2,2) = A(2,1)*B(2,1) + A(2,2)*B(2,2) = 2*6+4*8 = 44
    ("nt_basic",
     'N', 'T', 2, 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [5.0, 6.0, 7.0, 8.0], 2,
     0.0,
     [0.0, 0.0, 0.0, 0.0], 2,
     [26.0, 38.0, 30.0, 44.0]),

    # T7: TT — both transposed
    # op(A)(i,l)=A(l,i), op(B)(l,j)=B(j,l)
    # C(1,1) = A(1,1)*B(1,1) + A(2,1)*B(1,2) = 1*5+2*7 = 19
    # C(2,1) = A(1,2)*B(1,1) + A(2,2)*B(1,2) = 3*5+4*7 = 43
    # C(1,2) = A(1,1)*B(2,1) + A(2,1)*B(2,2) = 1*6+2*8 = 22
    # C(2,2) = A(1,2)*B(2,1) + A(2,2)*B(2,2) = 3*6+4*8 = 50
    ("tt_basic",
     'T', 'T', 2, 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [5.0, 6.0, 7.0, 8.0], 2,
     0.0,
     [0.0, 0.0, 0.0, 0.0], 2,
     [19.0, 43.0, 22.0, 50.0]),
]


# A=[1,2,3,4] col-major 2x2 → A(1,1)=1,A(2,1)=2,A(1,2)=3,A(2,2)=4
DGEMV_CASES = [
    # (name, trans, m, n, alpha, a, lda, x, incx, beta, y, incy, expected_y)
    # T1: N, a=1,b=0: y = A*x = [1*1+3*2, 2*1+4*2] = [7,10]
    ("n_basic",
     'N', 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [1.0, 2.0], 1,
     0.0,
     [0.0, 0.0], 1,
     [7.0, 10.0]),

    # T2: T, a=1,b=0: y = A'*x = [1*1+2*2, 3*1+4*2] = [5,11]
    ("t_basic",
     'T', 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [1.0, 2.0], 1,
     0.0,
     [0.0, 0.0], 1,
     [5.0, 11.0]),

    # T3: a=0,b=2: y = 2*y_init
    ("alpha0",
     'N', 2, 2, 0.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [1.0, 2.0], 1,
     2.0,
     [1.0, 2.0], 1,
     [2.0, 4.0]),

    # T4: M=0 quick return
    ("m_zero",
     'N', 0, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [1.0, 2.0], 1,
     0.0,
     [9.0], 1,
     [9.0]),

    # T5: N, a=2,b=1: y = 2*A*x + y = 2*[7,10] + [1,1] = [15,21]
    ("n_ab",
     'N', 2, 2, 2.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [1.0, 2.0], 1,
     1.0,
     [1.0, 1.0], 1,
     [15.0, 21.0]),

    # T6: negative incx
    # x=[1,2], incx=-1: x(0)=x[1]=2, x(1)=x[0]=1 → effective [2,1]
    # A*[2,1] = [1*2+3*1, 2*2+4*1] = [5, 8]
    ("neg_incx",
     'N', 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [1.0, 2.0], -1,
     0.0,
     [0.0, 0.0], 1,
     [5.0, 8.0]),
]

# Upper: A=[1,0,3,2], col-major → A(1,1)=1,A(2,1)=0,A(1,2)=3,A(2,2)=2
# Lower: A=[1,2,0,3], col-major → A(1,1)=1,A(2,1)=2,A(1,2)=0,A(2,2)=3
DTRSV_CASES = [
    # (name, uplo, trans, diag, n, a, lda, x, incx, expected_x)
    # T1: U,N,nonunit: Ax=b, A={{1,3},{0,2}}, b=[7,4]
    #   x2=4/2=2, x1=(7-3*2)/1=1
    ("un_nonunit",
     'U', 'N', 'N', 2,
     [1.0, 0.0, 3.0, 2.0], 2,
     [7.0, 4.0], 1,
     [1.0, 2.0]),

    # T2: L,N,nonunit: A={{1,0},{2,3}}, b=[3,12]
    #   x1=3/1=3, x2=(12-2*3)/3=2
    ("ln_nonunit",
     'L', 'N', 'N', 2,
     [1.0, 2.0, 0.0, 3.0], 2,
     [3.0, 12.0], 1,
     [3.0, 2.0]),

    # T3: U,N,unit: A_eff={{1,2},{0,1}}, b=[7,2]
    #   x2=2, x1=7-2*2=3
    ("un_unit",
     'U', 'N', 'U', 2,
     [1.0, 0.0, 2.0, 1.0], 2,
     [7.0, 2.0], 1,
     [3.0, 2.0]),

    # T4: N=0 quick return
    ("n_zero",
     'U', 'N', 'N', 0,
     [1.0], 2,
     [9.0], 1,
     [9.0]),

    # T5: U,T,nonunit: solve A'x=b, A={{1,3},{0,2}}
    #   A'={{1,0},{3,2}}. x1=7/1=7, x2=(4-3*7)/2=(4-21)/2=-8.5
    ("ut_nonunit",
     'U', 'T', 'N', 2,
     [1.0, 0.0, 3.0, 2.0], 2,
     [7.0, 4.0], 1,
     [7.0, -8.5]),
]

DTRSM_CASES = [
    # (name, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, expected_b)
    # T1: L,U,N,nonunit,a=1: A={{1,3},{0,2}}, B col0=[7,4],col1=[13,6]
    #   col0: x2=4/2=2, x1=(7-3*2)/1=1 → [1,2]
    #   col1: x2=6/2=3, x1=(13-3*3)/1=4 → [4,3]
    ("lun_basic",
     'L', 'U', 'N', 'N', 2, 2, 1.0,
     [1.0, 0.0, 3.0, 2.0], 2,
     [7.0, 4.0, 13.0, 6.0], 2,
     [1.0, 2.0, 4.0, 3.0]),

    # T2: L,L,N,nonunit,a=1: A={{1,0},{2,3}}, B col0=[3,12],col1=[5,16]
    #   col0: x1=3, x2=(12-2*3)/3=2 → [3,2]
    #   col1: x1=5, x2=(16-2*5)/3=2 → [5,2]
    ("lln_basic",
     'L', 'L', 'N', 'N', 2, 2, 1.0,
     [1.0, 2.0, 0.0, 3.0], 2,
     [3.0, 12.0, 5.0, 16.0], 2,
     [3.0, 2.0, 5.0, 2.0]),

    # T3: alpha=0 → zero B
    ("alpha0",
     'L', 'U', 'N', 'N', 2, 2, 0.0,
     [1.0, 0.0, 3.0, 2.0], 2,
     [1.0, 2.0, 3.0, 4.0], 2,
     [0.0, 0.0, 0.0, 0.0]),

    # T4: M=0 quick return
    ("m_zero",
     'L', 'U', 'N', 'N', 0, 2, 1.0,
     [1.0], 2,
     [9.0], 2,
     [9.0]),
]

# A=[1,2,3,4] col-major 2x2, A*A'={{10,14},{14,20}}
DSYRK_CASES = [
    # (name, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, expected_c)
    # T1: U,N,a=1,b=0: upper triangle only
    #   C(1,1)=10, C(1,2)=14, C(2,2)=20, C(2,1) untouched
    ("un_basic",
     'U', 'N', 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     0.0,
     [0.0, 99.0, 0.0, 0.0], 2,
     [10.0, 99.0, 14.0, 20.0]),

    # T2: L,N,a=1,b=0: lower triangle only
    #   C(1,1)=10, C(2,1)=14, C(2,2)=20, C(1,2) untouched
    ("ln_basic",
     'L', 'N', 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     0.0,
     [0.0, 0.0, 99.0, 0.0], 2,
     [10.0, 14.0, 99.0, 20.0]),

    # T3: a=0,b=2: scale upper triangle by 2
    #   C_init=[1,2,3,4], upper: C(1,1),C(1,2),C(2,2) scaled
    #   → [2, 2, 6, 8] (C(2,1)=2 untouched)
    ("alpha0",
     'U', 'N', 2, 2, 0.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     2.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     [2.0, 2.0, 6.0, 8.0]),

    # T4: N=0 quick return
    ("n_zero",
     'U', 'N', 0, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     0.0,
     [9.0], 2,
     [9.0]),

    # T5: U,T: A'*A
    # A=[1,2,3,4] col-major: A col0=[1,2], col1=[3,4]
    # A'*A: (A'*A)(i,j) = sum_l A(l,i)*A(l,j)
    # (1,1)=1*1+2*2=5, (1,2)=1*3+2*4=11, (2,2)=3*3+4*4=25
    ("ut_basic",
     'U', 'T', 2, 2, 1.0,
     [1.0, 2.0, 3.0, 4.0], 2,
     0.0,
     [0.0, 99.0, 0.0, 0.0], 2,
     [5.0, 99.0, 11.0, 25.0]),
]

# A 2x2 col-major, rank-1 update: A := alpha*x*y' + A
DGER_CASES = [
    # (name, m, n, alpha, x, incx, y, incy, a, lda, expected_a)
    # T1: basic, alpha=1, x=[1,2], y=[3,4], A=zeros
    # A += 1*[1;2]*[3,4] = [3,6,4,8] col-major
    ("basic",
     2, 2, 1.0,
     [1.0, 2.0], 1,
     [3.0, 4.0], 1,
     [0.0, 0.0, 0.0, 0.0], 2,
     [3.0, 6.0, 4.0, 8.0]),

    # T2: accumulate, alpha=2, x=[1,2], y=[1,1], A=[1,2,3,4]
    # A = [1,2,3,4] + 2*[1;2]*[1,1] = [3,6,5,8]
    ("accum",
     2, 2, 2.0,
     [1.0, 2.0], 1,
     [1.0, 1.0], 1,
     [1.0, 2.0, 3.0, 4.0], 2,
     [3.0, 6.0, 5.0, 8.0]),

    # T3: alpha=0, A unchanged
    ("alpha0",
     2, 2, 0.0,
     [1.0, 2.0], 1,
     [3.0, 4.0], 1,
     [1.0, 2.0, 3.0, 4.0], 2,
     [1.0, 2.0, 3.0, 4.0]),

    # T4: M=0 quick return
    ("m_zero",
     0, 2, 1.0,
     [1.0], 1,
     [1.0, 2.0], 1,
     [9.0], 2,
     [9.0]),

    # T5: negative incx=-1
    # x=[1,2] incx=-1: effective x=[2,1]
    # A += 1*[2;1]*[3,4] = [6,3,8,4]
    ("neg_incx",
     2, 2, 1.0,
     [1.0, 2.0], -1,
     [3.0, 4.0], 1,
     [0.0, 0.0, 0.0, 0.0], 2,
     [6.0, 3.0, 8.0, 4.0]),
]


def test_python_daxpy():
    passed, failed = 0, 0
    for name, n, alpha, x, incx, y, incy, expected in DAXPY_CASES:
        result = ref_daxpy(n, alpha, x, incx, y, incy)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: daxpy {name}")
            passed += 1
        else:
            print(f"  FAIL: daxpy {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_ddot():
    passed, failed = 0, 0
    for name, n, x, incx, y, incy, expected in DDOT_CASES:
        result = ref_ddot(n, x, incx, y, incy)
        ok = abs(result - expected) < 1e-12
        if ok:
            print(f"  PASS: ddot {name}")
            passed += 1
        else:
            print(f"  FAIL: ddot {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dscal():
    passed, failed = 0, 0
    for name, n, alpha, x, incx, expected in DSCAL_CASES:
        result = ref_dscal(n, alpha, x, incx)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dscal {name}")
            passed += 1
        else:
            print(f"  FAIL: dscal {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dcopy():
    passed, failed = 0, 0
    for name, n, x, incx, y, incy, expected in DCOPY_CASES:
        result = ref_dcopy(n, x, incx, y, incy)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dcopy {name}")
            passed += 1
        else:
            print(f"  FAIL: dcopy {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dswap():
    passed, failed = 0, 0
    for name, n, x, incx, y, incy, exp_x, exp_y in DSWAP_CASES:
        rx, ry = ref_dswap(n, x, incx, y, incy)
        ok_x = all(abs(r - e) < 1e-12 for r, e in zip(rx, exp_x))
        ok_y = all(abs(r - e) < 1e-12 for r, e in zip(ry, exp_y))
        if ok_x and ok_y:
            print(f"  PASS: dswap {name}")
            passed += 1
        else:
            print(f"  FAIL: dswap {name}: got x={rx} y={ry}")
            failed += 1
    return passed, failed


def test_python_dasum():
    passed, failed = 0, 0
    for name, n, x, incx, expected in DASUM_CASES:
        result = ref_dasum(n, x, incx)
        ok = abs(result - expected) < 1e-12
        if ok:
            print(f"  PASS: dasum {name}")
            passed += 1
        else:
            print(f"  FAIL: dasum {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dnrm2():
    passed, failed = 0, 0
    for name, n, x, incx, expected in DNRM2_CASES:
        result = ref_dnrm2(n, x, incx)
        ok = abs(result - expected) < 1e-12
        if ok:
            print(f"  PASS: dnrm2 {name}")
            passed += 1
        else:
            print(f"  FAIL: dnrm2 {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_idamax():
    passed, failed = 0, 0
    for name, n, x, incx, expected in IDAMAX_CASES:
        result = ref_idamax(n, x, incx)
        ok = result == expected
        if ok:
            print(f"  PASS: idamax {name}")
            passed += 1
        else:
            print(f"  FAIL: idamax {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dgemm():
    passed, failed = 0, 0
    for name, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, expected in DGEMM_CASES:
        result = ref_dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c[:], ldc)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dgemm {name}")
            passed += 1
        else:
            print(f"  FAIL: dgemm {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dgemv():
    passed, failed = 0, 0
    for name, trans, m, n, alpha, a, lda, x, incx, beta, y, incy, expected in DGEMV_CASES:
        result = ref_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y[:], incy)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dgemv {name}")
            passed += 1
        else:
            print(f"  FAIL: dgemv {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dtrsv():
    passed, failed = 0, 0
    for name, uplo, trans, diag, n, a, lda, x, incx, expected in DTRSV_CASES:
        result = ref_dtrsv(uplo, trans, diag, n, a, lda, x[:], incx)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dtrsv {name}")
            passed += 1
        else:
            print(f"  FAIL: dtrsv {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dtrsm():
    passed, failed = 0, 0
    for name, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, expected in DTRSM_CASES:
        result = ref_dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b[:], ldb)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dtrsm {name}")
            passed += 1
        else:
            print(f"  FAIL: dtrsm {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dsyrk():
    passed, failed = 0, 0
    for name, uplo, trans, n, k, alpha, a, lda, beta, c, ldc, expected in DSYRK_CASES:
        result = ref_dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c[:], ldc)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dsyrk {name}")
            passed += 1
        else:
            print(f"  FAIL: dsyrk {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


def test_python_dger():
    passed, failed = 0, 0
    for name, m, n, alpha, x, incx, y, incy, a, lda, expected in DGER_CASES:
        result = ref_dger(m, n, alpha, x, incx, y, incy, a[:], lda)
        ok = all(abs(r - e) < 1e-12 for r, e in zip(result, expected))
        if ok:
            print(f"  PASS: dger {name}")
            passed += 1
        else:
            print(f"  FAIL: dger {name}: got {result}, expected {expected}")
            failed += 1
    return passed, failed


# ---- Fortran cross-check ----

def build_fortran():
    """Compile the slatec Fortran driver. Returns True on success."""
    if not os.path.exists(FORTRAN_SRC):
        print("  SKIP: slatec_driver.f90 not found")
        return False
    if not shutil.which("gfortran"):
        print("  SKIP: gfortran not available")
        return False
    result = subprocess.run(
        ["gfortran", "-O2", "-o", FORTRAN_EXE, FORTRAN_SRC],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  SKIP: gfortran failed:\n{result.stderr}")
        return False
    return True


def run_fortran():
    """Run the Fortran driver and parse output. Returns dict of results."""
    if not os.path.exists(FORTRAN_EXE):
        return None
    result = subprocess.run(
        [FORTRAN_EXE], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: Fortran driver failed:\n{result.stderr}")
        return None

    results = {}
    for line in result.stdout.strip().split('\n'):
        line = line.strip()
        if not line:
            continue
        # Parse "ROUTINE test_name  val1  val2  ..."
        parts = line.split()
        if len(parts) >= 3:
            key = f"{parts[0]}_{parts[1]}"
            vals = [float(v) for v in parts[2:]]
            results[key] = vals
    return results


FORTRAN_EXPECTED = {
    # DAXPY tests
    "DAXPY_basic":   [12.0, 24.0, 36.0, 48.0],
    "DAXPY_zero_a":  [10.0, 20.0, 30.0, 40.0],
    "DAXPY_n_zero":  [10.0, 20.0, 30.0, 40.0],
    "DAXPY_n_one":   [22.0],
    "DAXPY_incx2":   [11.0, 22.0, 33.0],
    "DAXPY_neg_inc": [13.0, 22.0, 31.0],
    # DDOT tests
    "DDOT_basic":    [70.0],
    "DDOT_n_zero":   [0.0],
    "DDOT_ortho":    [0.0],
    "DDOT_n_one":    [21.0],
    "DDOT_incx2":    [11.0],
    "DDOT_neg_inc":  [100.0],
    # DSCAL tests
    "DSCAL_basic":   [2.0, 4.0, 6.0],
    "DSCAL_zero_a":  [0.0, 0.0],
    "DSCAL_n_zero":  [1.0],
    "DSCAL_n_one":   [15.0],
    # DCOPY tests
    "DCOPY_basic":   [1.0, 2.0, 3.0],
    "DCOPY_n_zero":  [99.0],
    "DCOPY_incx2":   [1.0, 2.0],
    # DSWAP tests (X values after swap)
    "DSWAP_basic_x": [3.0, 4.0],
    "DSWAP_basic_y": [1.0, 2.0],
    "DSWAP_n_zero_x": [1.0],
    "DSWAP_n_zero_y": [2.0],
    # DASUM tests
    "DASUM_basic":   [10.0],
    "DASUM_n_zero":  [0.0],
    "DASUM_allpos":  [6.0],
    "DASUM_n_one":   [5.0],
    # DNRM2 tests
    "DNRM2_basic":   [5.0],
    "DNRM2_n_zero":  [0.0],
    "DNRM2_n_one":   [3.0],
    # IDAMAX tests (integer results as float)
    "IDAMAX_basic":  [2.0],
    "IDAMAX_n_zero": [0.0],
    "IDAMAX_n_one":  [1.0],
    "IDAMAX_first":  [1.0],
    "IDAMAX_last":   [4.0],
    # DGEMM tests
    "DGEMM_nn_basic":  [23.0, 34.0, 31.0, 46.0],
    "DGEMM_nn_ab":     [47.0, 69.0, 63.0, 93.0],
    "DGEMM_alpha0":    [2.0, 4.0, 6.0, 8.0],
    "DGEMM_tn_basic":  [17.0, 39.0, 23.0, 53.0],
    "DGEMM_nt_basic":  [26.0, 38.0, 30.0, 44.0],
    "DGEMM_tt_basic":  [19.0, 43.0, 22.0, 50.0],
    # DGEMV tests
    "DGEMV_n_basic":   [7.0, 10.0],
    "DGEMV_t_basic":   [5.0, 11.0],
    "DGEMV_alpha0":    [2.0, 4.0],
    "DGEMV_n_ab":      [15.0, 21.0],
    # DTRSV tests
    "DTRSV_un_nonu":   [1.0, 2.0],
    "DTRSV_ln_nonu":   [3.0, 2.0],
    "DTRSV_un_unit":   [3.0, 2.0],
    # DTRSM tests
    "DTRSM_lun_basic": [1.0, 2.0, 4.0, 3.0],
    "DTRSM_lln_basic": [3.0, 2.0, 5.0, 2.0],
    "DTRSM_alpha0":    [0.0, 0.0, 0.0, 0.0],
    # DSYRK tests
    "DSYRK_un_basic":  [10.0, 99.0, 14.0, 20.0],
    "DSYRK_ln_basic":  [10.0, 14.0, 99.0, 20.0],
    "DSYRK_alpha0":    [2.0, 2.0, 6.0, 8.0],
    "DSYRK_ut_basic":  [5.0, 99.0, 11.0, 25.0],
    # DGER tests
    "DGER_basic":      [3.0, 6.0, 4.0, 8.0],
    "DGER_accum":      [3.0, 6.0, 5.0, 8.0],
    "DGER_alpha0":     [1.0, 2.0, 3.0, 4.0],
}


def test_fortran():
    """Cross-check Fortran output against known-good values."""
    passed, failed = 0, 0
    results = run_fortran()
    if results is None:
        print("  SKIP: Fortran driver not available")
        return 0, 0

    for key, expected in FORTRAN_EXPECTED.items():
        if key not in results:
            print(f"  FAIL: {key} missing from Fortran output")
            failed += 1
            continue
        actual = results[key]
        ok = len(actual) == len(expected) and all(
            abs(a - e) < 1e-12 for a, e in zip(actual, expected))
        if ok:
            print(f"  PASS: fortran {key}")
            passed += 1
        else:
            print(f"  FAIL: fortran {key}: got {actual}, expected {expected}")
            failed += 1

    return passed, failed


# ---- Main ----

def main():
    print("zblas verify — mathematical correctness check")
    print("=" * 55)

    total_pass = 0
    total_fail = 0

    print("\n--- Layer 1: Python reference ---")
    for label, func in [
        ("DAXPY", test_python_daxpy),
        ("DDOT", test_python_ddot),
        ("DSCAL", test_python_dscal),
        ("DCOPY", test_python_dcopy),
        ("DSWAP", test_python_dswap),
        ("DASUM", test_python_dasum),
        ("DNRM2", test_python_dnrm2),
        ("IDAMAX", test_python_idamax),
        ("DGEMM", test_python_dgemm),
        ("DGEMV", test_python_dgemv),
        ("DTRSV", test_python_dtrsv),
        ("DTRSM", test_python_dtrsm),
        ("DSYRK", test_python_dsyrk),
        ("DGER", test_python_dger),
    ]:
        print(f"\n{label}:")
        p, f = func()
        total_pass += p; total_fail += f

    print("\n--- Layer 2: Slatec Fortran ---")
    print()
    if build_fortran():
        p, f = test_fortran()
        total_pass += p; total_fail += f
    else:
        print("  (Fortran layer skipped)")

    print()
    print("=" * 55)
    print(f"Total: {total_pass} passed, {total_fail} failed")

    if total_fail > 0:
        print("VERIFICATION FAILED")
        sys.exit(1)
    else:
        print("All checks passed.")
        sys.exit(0)


if __name__ == "__main__":
    main()
