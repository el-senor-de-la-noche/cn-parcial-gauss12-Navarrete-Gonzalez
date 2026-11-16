#Integrantes: Elias Gonzalez, Daniel Navarrete
#Ramo: Computacion Numerica
#Profesor: Juan Cabezas

"""Verificación independiente: solución exacta, residuo en precisión extendida,
refinamiento iterativo, condicionamiento y control con biblioteca.

La precisión extendida se hace con fractions.Fraction (aritmética racional
exacta): A y b se leen de las cadenas decimales del CSV, sin pérdida de cifras.
"""

import math
from fractions import Fraction
from decimal import Decimal

import mpmath

from src.gauss import resolver

EPS = 2.220446049250313e-16   # epsilon de máquina, doble precisión
mpmath.mp.dps = 50


def _frac(cadena):
    return Fraction(Decimal(cadena))


def a_fracciones(A_str, b_str):
    A = [[_frac(v) for v in fila] for fila in A_str]
    b = [_frac(v) for v in b_str]
    return A, b


# --------------------------------------------------------------------------
# Solución exacta (referencia para contar dígitos correctos)
# Eliminación propia con Fraction: no usa rutinas de biblioteca para resolver.
# --------------------------------------------------------------------------
def solucion_exacta(A_str, b_str):
    A, b = a_fracciones(A_str, b_str)
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0)
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= m * M[k][j]
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


# --------------------------------------------------------------------------
# Errores y dígitos
# --------------------------------------------------------------------------
def error_relativo_inf(x, x_exacto):
    """||x - x*||_inf / ||x*||_inf, calculado en precisión extendida."""
    num = max(abs(Fraction(xi) - xe) for xi, xe in zip(x, x_exacto))
    den = max(abs(xe) for xe in x_exacto)
    return float(num / den)


def digitos_correctos(xi, xe, maximo=15):
    """Dígitos significativos correctos de xi respecto de xe (entero)."""
    err = abs(Fraction(xi) - xe)
    if err == 0:
        return maximo
    rel = err / abs(xe)
    d = -mpmath.log10(mpmath.mpf(rel.numerator) / mpmath.mpf(rel.denominator))
    return max(0, min(maximo, int(mpmath.floor(d))))


# --------------------------------------------------------------------------
# Residuo en precisión extendida
# --------------------------------------------------------------------------
def residuo_exacto(A_str, b_str, x):
    """r = b - A x con A, b exactos y x (float) tomado tal cual es."""
    A, b = a_fracciones(A_str, b_str)
    xf = [Fraction(v) for v in x]
    return [b[i] - sum(A[i][j] * xf[j] for j in range(len(A))) for i in range(len(A))]


def norma_inf(v):
    return max(abs(vi) for vi in v)


def fraccion_a_mp(f):
    return mpmath.mpf(f.numerator) / mpmath.mpf(f.denominator)


def formatear(f, cifras=6):
    return mpmath.nstr(fraccion_a_mp(f), cifras)


# --------------------------------------------------------------------------
# Refinamiento iterativo (residuo exacto, corrección con L y U reutilizados)
# --------------------------------------------------------------------------
def refinar(res, A_str, b_str, x0, x_exacto, pasadas=3):
    x = x0[:]
    historia = []
    for k in range(pasadas + 1):
        r = residuo_exacto(A_str, b_str, x)
        historia.append({
            "pasada": k,
            "x": x[:],
            "norma_r": norma_inf(r),
            "error_rel": error_relativo_inf(x, x_exacto),
            "min_digitos": min(digitos_correctos(xi, xe)
                               for xi, xe in zip(x, x_exacto)),
        })
        if k == pasadas:
            break
        d = resolver(res, [float(ri) for ri in r])
        x = [xi + di for xi, di in zip(x, d)]
    return historia


# --------------------------------------------------------------------------
# Condicionamiento
# --------------------------------------------------------------------------
def condicionamiento(A):
    import numpy as np
    kappa = float(np.linalg.cond(np.array(A, dtype=float), 2))
    return {"kappa2": kappa, "cota": kappa * EPS,
            "digitos_cota": -math.log10(kappa * EPS)}


# --------------------------------------------------------------------------
# Control con biblioteca (NO sustituye la implementación propia)
# --------------------------------------------------------------------------
def control_numpy(A, b, x_exacto):
    import numpy as np
    x = np.linalg.solve(np.array(A, dtype=float), np.array(b, dtype=float))
    x = [float(v) for v in x]
    return {
        "x": x,
        "error_rel": error_relativo_inf(x, x_exacto),
        "min_digitos": min(digitos_correctos(xi, xe) for xi, xe in zip(x, x_exacto)),
    }
