# Integrantes: Elias Gonzalez, Daniel Navarrete
# Ramo: Computacion Numerica
# Profesor: Juan Cabezas


# ARCHIVO: ELIMINACIÓN GAUSSIANA
# Resuelve sistemas lineales mediante eliminación gaussiana, tanto sin
# pivoteo como con pivoteo parcial escalado, y registra los flops efectuados.
# Al eliminar, cada multiplicador se calcula una sola vez por fila y se
# reutiliza. La matriz L guarda dichos multiplicadores y U queda triangular.
# El recorrido empieza en la columna del pivote, sin rehacer columnas previas.
# No se calcula la matriz inversa; los cálculos usan float de doble precisión
# (IEEE 754, con 53 bits de precisión en la mantisa).

import csv
import math
from dataclasses import dataclass, field


#
# LECTURA DE DATOS: se conservan los decimales originales para la verificación
# exacta, sin truncarlos ni redondearlos al abrir los archivos.
#
# Lee una matriz desde un CSV y conserva cada valor como texto.
# Esto permite verificar después los resultados con aritmética exacta.
def leer_matriz_csv(ruta):
    with open(ruta, newline="", encoding="utf-8") as f:
        filas = [[v.strip() for v in fila] for fila in csv.reader(f) if fila]
    return filas


# Lee la primera columna del CSV para construir el vector b, también
# sin convertir ni redondear los decimales durante la lectura.
def leer_vector_csv(ruta):
    with open(ruta, newline="", encoding="utf-8") as f:
        return [fila[0].strip() for fila in csv.reader(f) if fila]


# Convierte los textos a float para efectuar el cálculo numérico.
# Acepta tanto una matriz (lista de filas) como un vector.
def a_float(matriz_o_vector):
    if matriz_o_vector and isinstance(matriz_o_vector[0], list):
        return [[float(v) for v in fila] for fila in matriz_o_vector]
    return [float(v) for v in matriz_o_vector]


#
# ESTRUCTURA DEL RESULTADO DE CADA CORRIDA
#
# Reúne el resultado del algoritmo en una sola estructura para facilitar
# la verificación, las trazas y el reporte de indicadores de estabilidad.
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


# Copia las filas para evitar modificar la matriz original al eliminar.
def _copia(M):
    return [fila[:] for fila in M]


# Cuenta cuántas posiciones bajo la diagonal son exactamente cero;
# permite comprobar si la eliminación logró triangularizar la matriz.
def ceros_bajo_diagonal(U):
    n = len(U)
    return sum(1 for i in range(n) for j in range(i) if U[i][j] == 0.0)


#
# ELIMINACIÓN GAUSSIANA
#
def eliminar(A, b, pivoteo, contador):
    # Transforma el sistema A*x=b hasta obtener una matriz triangular superior.
    # El parámetro pivoteo indica si se conservan las filas o si se selecciona
    # un pivote proporcionalmente grande según la escala de cada fila.
    # Se guardan los multiplicadores en L, los cambios de filas y la solución.
    if pivoteo not in ("ninguno", "parcial_escalado"):
        raise ValueError("pivoteo debe ser 'ninguno' o 'parcial_escalado'")
    # Se trabaja sobre copias de A y b. L almacena los multiplicadores y
    # perm recuerda el orden original de las filas después de intercambiarlas.
    n = len(A)
    A = _copia(A)
    b = b[:]
    L = [[0.0] * n for _ in range(n)]
    perm = list(range(n))

    # Inicializa los indicadores de estabilidad que se van observando:
    # crecimiento de elementos, multiplicador máximo y cambios de filas.
    max_ini = max(abs(v) for fila in A for v in fila)
    max_dur = max_ini
    max_mult = 0.0
    bits_max = 0.0
    intercambios = 0
    trazas = []

    # Cada escala guarda el mayor valor absoluto de su fila original.
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

    # Cada vuelta fija el pivote de la columna k y elimina los términos
    # que quedan debajo de su diagonal.
    for k in range(n - 1):
        # BÚSQUEDA DEL PIVOTE: se comparan razones escaladas por fila.
        if pivoteo == "parcial_escalado":
            # Se parte de la fila actual y se compara |a_ik| / escala_i.
            # La fila con la mayor razón pasa a ser el nuevo pivote.
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
            # Si se eligió otra fila, se intercambian A, b, L, escalas y perm
            # para mantener consistente la factorización y el orden de variables.
            if p != k:
                A[k], A[p] = A[p], A[k]
                b[k], b[p] = b[p], b[k]
                L[k], L[p] = L[p], L[k]
                escala[k], escala[p] = escala[p], escala[k]
                perm[k], perm[p] = perm[p], perm[k]
                contador.pivote_intercambio()
                intercambios += 1

        # Un pivote nulo impide dividir: se informa el paso problemático.
        if A[k][k] == 0.0:
            raise ZeroDivisionError(f"pivote nulo en el paso {k + 1}")

        # ELIMINACIÓN: se anulan los valores bajo el pivote actual.
        # El contador separa el costo de este paso de los demás.
        contador.fase_eliminacion(k)
        pivote = A[k][k]
        for i in range(k + 1, n):
            # Se calcula un multiplicador por fila y se reutiliza en todo el barrido.
            contador.div()
            m = A[i][k] / pivote
            # Guarda el factor que se usará para anular el elemento A[i][k].
            L[i][k] = m
            if abs(m) > max_mult:
                max_mult = abs(m)

            # Se actualizan las columnas desde k; los valores anteriores no se modifican.
            for j in range(k, n):
                contador.mul()
                prod = m * A[k][j]
                contador.sumres()
                nuevo = A[i][j] - prod
                # Se observa la relación entre el término restado y el anterior
                # para estimar posibles pérdidas de precisión por cancelación.
                if j > k and prod != 0.0 and A[i][j] != 0.0:
                    bits = math.log2(abs(prod) / abs(A[i][j]))
                    if bits > bits_max:
                        bits_max = bits
                A[i][j] = nuevo
                if j > k and abs(nuevo) > max_dur:
                    max_dur = abs(nuevo)
            # El elemento eliminado se fija explícitamente en cero.
            A[i][k] = 0.0

            # El mismo multiplicador se aplica al término independiente b[i].
            contador.mul()
            prod_b = m * b[k]
            contador.sumres()
            b[i] = b[i] - prod_b

        # Conserva una copia del estado tras el paso para el archivo de trazas.
        trazas.append((k + 1, _copia(A), b[:]))

    # SUSTITUCIÓN REGRESIVA: despeja x desde la última ecuación.
    contador.fase_regresiva()
    # Reserva el vector de soluciones para la sustitución regresiva.
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = b[i]
        for j in range(i + 1, n):
            contador.mul()
            contador.sumres()
            s -= A[i][j] * x[j]
        contador.div()
        x[i] = s / A[i][i]

    # Entrega U, L, la permutación, la solución y los indicadores
    # acumulados para que otros módulos puedan analizarlos.
    return Resultado(
        n=n, pivoteo=pivoteo, U=A, L=L, perm=perm, c=b, x=x, trazas=trazas,
        max_multiplicador=max_mult,
        factor_crecimiento=max_dur / max_ini,
        max_elemento_inicial=max_ini,
        max_elemento_durante=max_dur,
        bits_perdidos_max=min(bits_max, 53.0),
        intercambios=intercambios,
    )


#
# RESOLUCIÓN CON FACTORES L Y U EXISTENTES
#
# Resuelve un nuevo sistema con el mismo A sin volver a factorizarlo.
# Primero aplica la permutación al lado derecho, luego resuelve L*y=P*b
# (sustitución progresiva) y finalmente U*x=y (regresiva).
def resolver(res, rhs, contador=None):
    n = res.n
    # El nuevo lado derecho debe respetar los intercambios guardados.
    pb = [rhs[res.perm[i]] for i in range(n)]
    if contador is not None:
        contador.fase_progresiva()
    # Resuelve la matriz triangular inferior L de arriba hacia abajo.
    # La diagonal de L es implícitamente unitaria.
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
    # Reserva el vector de soluciones para la sustitución regresiva.
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
