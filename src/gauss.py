#Integrantes: Elias Gonzalez, Daniel Navarrete
#Ramo: Computacion Numerica
#Profesor: Juan Cabezas


"""Eliminación gaussiana con y sin pivoteo, con contador de flops.

- Una sola factorización PA = LU (L guarda los multiplicadores, cada m(i,k)
  se calcula UNA vez y se reutiliza en toda la fila).
- El barrido del paso k va de la columna k en adelante (no toca lo ya anulado).
- Nunca se calcula la inversa de A.
- Todo en doble precisión (float de Python, IEEE 754, 53 bits de mantisa).
"""

import csv
import math
from dataclasses import dataclass, field


# --------------------------------------------------------------------------
# Lectura de datos (se conservan las cadenas originales para la verificación
# exacta: nada se trunca ni se redondea al leer).
# --------------------------------------------------------------------------
def leer_matriz_csv(ruta):
    with open(ruta, newline="", encoding="utf-8") as f:
        filas = [[v.strip() for v in fila] for fila in csv.reader(f) if fila]
    return filas


def leer_vector_csv(ruta):
    with open(ruta, newline="", encoding="utf-8") as f:
        return [fila[0].strip() for fila in csv.reader(f) if fila]


def a_float(matriz_o_vector):
    if matriz_o_vector and isinstance(matriz_o_vector[0], list):
        return [[float(v) for v in fila] for fila in matriz_o_vector]
    return [float(v) for v in matriz_o_vector]


# --------------------------------------------------------------------------
# Resultado de una corrida
# --------------------------------------------------------------------------
@dataclass
class Resultado:
    n: int
    pivoteo: str
    U: list
    L: list
    perm: list
    c: list                      # b transformado (= L^-1 P b)
    x: list
    trazas: list = field(default_factory=list)   # (paso, A, b) tras cada paso
    max_multiplicador: float = 0.0
    factor_crecimiento: float = 0.0
    max_elemento_inicial: float = 0.0
    max_elemento_durante: float = 0.0
    bits_perdidos_max: float = 0.0
    intercambios: int = 0


def _copia(M):
    return [fila[:] for fila in M]


def ceros_bajo_diagonal(U):
    n = len(U)
    return sum(1 for i in range(n) for j in range(i) if U[i][j] == 0.0)


# --------------------------------------------------------------------------
# Eliminación
# --------------------------------------------------------------------------
def eliminar(A, b, pivoteo, contador):
    """Triangulariza [A|b]. pivoteo: 'ninguno' o 'parcial_escalado'."""
    if pivoteo not in ("ninguno", "parcial_escalado"):
        raise ValueError("pivoteo debe ser 'ninguno' o 'parcial_escalado'")
    n = len(A)
    A = _copia(A)
    b = b[:]
    L = [[0.0] * n for _ in range(n)]
    perm = list(range(n))

    max_ini = max(abs(v) for fila in A for v in fila)
    max_dur = max_ini
    max_mult = 0.0
    bits_max = 0.0
    intercambios = 0
    trazas = []

    # Escalas de fila (se calculan una vez, sobre la matriz original).
    escala = None
    if pivoteo == "parcial_escalado":
        escala = []
        for i in range(n):
            s = abs(A[i][0])
            for j in range(1, n):
                contador.pivote_escala()
                if abs(A[i][j]) > s:
                    s = abs(A[i][j])
            escala.append(s)

    for k in range(n - 1):
        # --- búsqueda de pivote (sobrecosto de estabilidad) ---------------
        if pivoteo == "parcial_escalado":
            p = k
            contador.pivote_razon()
            mejor = abs(A[k][k]) / escala[k]
            for i in range(k + 1, n):
                contador.pivote_razon()
                razon = abs(A[i][k]) / escala[i]
                contador.pivote_busqueda()
                if razon > mejor:
                    mejor = razon
                    p = i
            if p != k:
                A[k], A[p] = A[p], A[k]
                b[k], b[p] = b[p], b[k]
                L[k], L[p] = L[p], L[k]
                escala[k], escala[p] = escala[p], escala[k]
                perm[k], perm[p] = perm[p], perm[k]
                contador.pivote_intercambio()
                intercambios += 1

        if A[k][k] == 0.0:
            raise ZeroDivisionError(f"pivote nulo en el paso {k + 1}")

        # --- eliminación del paso k ---------------------------------------
        contador.fase_eliminacion(k)
        pivote = A[k][k]
        for i in range(k + 1, n):
            # multiplicador: una sola división por fila, reutilizado después
            contador.div()
            m = A[i][k] / pivote
            L[i][k] = m
            if abs(m) > max_mult:
                max_mult = abs(m)

            # barrido desde la columna k (la columna k queda exactamente 0)
            for j in range(k, n):
                contador.mul()
                prod = m * A[k][j]
                contador.sumres()
                nuevo = A[i][j] - prod
                if j > k and prod != 0.0 and A[i][j] != 0.0:
                    bits = math.log2(abs(prod) / abs(A[i][j]))
                    if bits > bits_max:
                        bits_max = bits
                A[i][j] = nuevo
                if j > k and abs(nuevo) > max_dur:
                    max_dur = abs(nuevo)
            A[i][k] = 0.0

            # actualización del lado derecho (con el mismo m)
            contador.mul()
            prod_b = m * b[k]
            contador.sumres()
            b[i] = b[i] - prod_b

        trazas.append((k + 1, _copia(A), b[:]))

    # --- sustitución regresiva ---------------------------------------------
    contador.fase_regresiva()
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = b[i]
        for j in range(i + 1, n):
            contador.mul()
            contador.sumres()
            s -= A[i][j] * x[j]
        contador.div()
        x[i] = s / A[i][i]

    return Resultado(
        n=n, pivoteo=pivoteo, U=A, L=L, perm=perm, c=b, x=x, trazas=trazas,
        max_multiplicador=max_mult,
        factor_crecimiento=max_dur / max_ini,
        max_elemento_inicial=max_ini,
        max_elemento_durante=max_dur,
        bits_perdidos_max=min(bits_max, 53.0),
        intercambios=intercambios,
    )


# --------------------------------------------------------------------------
# Resolver otro lado derecho reutilizando L y U (sin refactorizar)
# --------------------------------------------------------------------------
def resolver(res, rhs, contador=None):
    n = res.n
    pb = [rhs[res.perm[i]] for i in range(n)]
    if contador is not None:
        contador.fase_progresiva()
    y = [0.0] * n
    for i in range(n):
        s = pb[i]
        for j in range(i):
            if contador is not None:
                contador.mul()
                contador.sumres()
            s -= res.L[i][j] * y[j]
        y[i] = s
    if contador is not None:
        contador.fase_regresiva()
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        for j in range(i + 1, n):
            if contador is not None:
                contador.mul()
                contador.sumres()
            s -= res.U[i][j] * x[j]
        if contador is not None:
            contador.div()
        x[i] = s / res.U[i][i]
    return x
