# Integrantes: Elias Gonzalez, Daniel Navarrete
# Ramo: Computacion Numerica
# Profesor: Juan Cabezas

# ARCHIVO: VERIFICACIÓN DE RESULTADOS
# Calcula una solución de referencia exacta, mide los errores y los residuos,
# aplica refinamiento iterativo y evalúa el condicionamiento del sistema.
# Fraction permite operar con fracciones exactas a partir de los valores
# decimales originales del CSV; así se evita redondear los datos al leerlos.
# La solución de NumPy se usa exclusivamente como control independiente.

import math
from fractions import Fraction
from decimal import Decimal

import mpmath

from src.gauss import resolver

# Se usa el épsilon de máquina de float64 para estimar la limitación
# de precisión del método en función del condicionamiento.
EPS = 2.220446049250313e-16   # epsilon de máquina, doble precisión
mpmath.mp.dps = 50


# Convierte el texto decimal primero a Decimal y luego a Fraction;
# así se conserva el valor escrito, sin pasar por float.
def _frac(cadena):
    return Fraction(Decimal(cadena))


# Crea copias exactas de la matriz y el vector a partir de los CSV.
def a_fracciones(A_str, b_str):
    A = [[_frac(v) for v in fila] for fila in A_str]
    b = [_frac(v) for v in b_str]
    return A, b


#
# SOLUCIÓN EXACTA: referencia para evaluar las cifras correctas.
# Se usa eliminación propia con Fraction, sin resolver con librerías.
#
# Calcula la solución de referencia con eliminación gaussiana hecha
# con Fraction, sin delegar la resolución a una biblioteca numérica.
def solucion_exacta(A_str, b_str):
    A, b = a_fracciones(A_str, b_str)
    n = len(A)
    # Forma la matriz aumentada para operar sobre A y b simultáneamente.
    M = [A[i][:] + [b[i]] for i in range(n)]
    # Si hace falta, coloca arriba una fila con pivote distinto de cero;
    # luego elimina los valores que están debajo de ese pivote.
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0)
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= m * M[k][j]
    # Obtiene las incógnitas en orden inverso con aritmética racional.
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / M[i][i]
    return x


#
# MEDIDAS DE ERROR Y DÍGITOS CORRECTOS
#
def error_relativo_inf(x, x_exacto):
    # Calcula el error relativo de la solución usando la norma infinito:
    # máximo error absoluto de las componentes dividido por el máximo valor
    # absoluto de la solución exacta. Fraction evita redondeos intermedios.
    # El numerador es el mayor error por componente; el denominador
    # es el valor máximo en módulo de la solución exacta.
    num = max(abs(Fraction(xi) - xe) for xi, xe in zip(x, x_exacto))
    den = max(abs(xe) for xe in x_exacto)
    return float(num / den)


def digitos_correctos(xi, xe, maximo=15):
    # Estima cuántos dígitos significativos de xi coinciden con el valor
    # exacto xe, utilizando -log10 del error relativo y un máximo de 15.
    err = abs(Fraction(xi) - xe)
    # Si no existe diferencia, se informa el máximo de dígitos permitido.
    if err == 0:
        return maximo
    # Los dígitos se obtienen del orden de magnitud del error relativo.
    rel = err / abs(xe)
    d = -mpmath.log10(mpmath.mpf(rel.numerator) / mpmath.mpf(rel.denominator))
    return max(0, min(maximo, int(mpmath.floor(d))))


#
# RESIDUO CALCULADO CON FRACCIONES EXACTAS
#
def residuo_exacto(A_str, b_str, x):
    # Evalúa r=b-A*x con A y b exactos y los valores float de x tal como
    # están almacenados, sin volver a redondear los datos de entrada.
    A, b = a_fracciones(A_str, b_str)
    xf = [Fraction(v) for v in x]
    return [b[i] - sum(A[i][j] * xf[j] for j in range(len(A))) for i in range(len(A))]


# Calcula la norma infinito: el mayor valor absoluto del vector.
def norma_inf(v):
    return max(abs(vi) for vi in v)


# Convierte una fracción exacta al formato de alta precisión de mpmath.
def fraccion_a_mp(f):
    return mpmath.mpf(f.numerator) / mpmath.mpf(f.denominator)


# Presenta las fracciones con una cantidad controlada de cifras.
def formatear(f, cifras=6):
    return mpmath.nstr(fraccion_a_mp(f), cifras)


#
# REFINAMIENTO ITERATIVO: utiliza L y U sin repetir la factorización.
#
# Mejora la solución inicial calculando el residuo exacto en cada pasada.
# Para obtener la corrección se reutilizan L y U; no se repite Gauss.
# Conserva el historial completo, incluida la solución inicial (pasada 0).
def refinar(res, A_str, b_str, x0, x_exacto, pasadas=3):
    x = x0[:]
    historia = []
    # Guarda para cada pasada el residuo, error relativo y dígitos válidos.
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
        # Se convierte el residuo a float para resolver la corrección con L y U.
        d = resolver(res, [float(ri) for ri in r])
        # Se suma la corrección a cada componente de la solución actual.
        x = [xi + di for xi, di in zip(x, d)]
    return historia


#
# CONDICIONAMIENTO: mide la sensibilidad del sistema a errores.
#
# Calcula el número de condición en norma 2 mediante NumPy.
# El producto kappa*EPS indica una cota aproximada de sensibilidad y
# -log10 de ese valor da una estimación de dígitos garantizables.
def condicionamiento(A):
    import numpy as np
    kappa = float(np.linalg.cond(np.array(A, dtype=float), 2))
    return {"kappa2": kappa, "cota": kappa * EPS,
            "digitos_cota": -math.log10(kappa * EPS)}


#
# CONTROL NUMPY: comparación independiente, no reemplaza el método propio.
#
# Resuelve con numpy.linalg.solve únicamente como verificación externa;
# luego mide su error respecto de la misma solución exacta.
def control_numpy(A, b, x_exacto):
    import numpy as np
    x = np.linalg.solve(np.array(A, dtype=float), np.array(b, dtype=float))
    x = [float(v) for v in x]
    return {
        "x": x,
        "error_rel": error_relativo_inf(x, x_exacto),
        "min_digitos": min(digitos_correctos(xi, xe) for xi, xe in zip(x, x_exacto)),
    }
