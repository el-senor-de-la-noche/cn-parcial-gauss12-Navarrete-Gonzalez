# Integrantes: Elias Gonzalez, Daniel Navarrete
# Ramo: Computacion Numerica
# Profesor: Juan Cabezas

# ARCHIVO: EJECUCIÓN PRINCIPAL DEL PROYECTO
# Coordina la lectura de datos, las dos versiones de Gauss, la verificación
# numérica, el refinamiento de la solución y la escritura de los resultados.
# Para ejecutarlo desde la raíz del proyecto: python src/main.py
# Los archivos generados se guardan en la carpeta resultados/.

import math
import platform
import sys
from pathlib import Path

import mpmath
import numpy

from src.conteo import (Contador, cota_referencia, dentro_de_tolerancia,
                    reporte_markdown)
from src.gauss import (a_float, ceros_bajo_diagonal, eliminar, leer_matriz_csv,
                   leer_vector_csv)
import src.verificacion as ver

# Define las rutas del proyecto a partir de la ubicación de este archivo,
# de modo que la ejecución no dependa de la carpeta actual.
RAIZ = Path(__file__).resolve().parent.parent
DATOS = RAIZ / "datos"
SALIDA = RAIZ / "resultados"
PASADAS_REFINAMIENTO = 3


# Da formato científico a cada fila de la matriz aumentada [A | b],
# alineando los valores para que los pasos se puedan revisar visualmente.
def formato_matriz(A, b):
    lineas = []
    for fila, bi in zip(A, b):
        celdas = " ".join(f"{v:>17.9e}" for v in fila)
        lineas.append(f"{celdas} | {bi:>17.9e}")
    return "\n".join(lineas)


# Escribe la matriz inicial y todos los estados de eliminación de cada
# corrida, junto con las permutaciones y la cantidad de ceros obtenidos.
def escribir_traza(ruta, A0, b0, corridas):
    n = len(A0)
    out = []
    for nombre, res in corridas:
        out.append("=" * 78)
        out.append(f"CORRIDA: {nombre}")
        out.append("=" * 78)
        out.append("Matriz aumentada [A | b] inicial:\n")
        out.append(formato_matriz(A0, b0))
        # Recorre los estados guardados después de cada columna eliminada.
        for paso, A, b in res.trazas:
            # Cuenta los ceros bajo la diagonal de las columnas ya procesadas.
            nulos = sum(1 for i in range(n) for j in range(min(i, paso))
                        if A[i][j] == 0.0)
            out.append("")
            out.append(f"--- Tras el paso {paso} de {n - 1} "
                       f"(ceros bajo la diagonal en columnas 1..{paso}: "
                       f"{nulos}) ---\n")
            out.append(formato_matriz(A, b))
        out.append("")
        out.append(f"Permutación de filas final (índices originales, base 1): "
                   f"{[p + 1 for p in res.perm]}")
        out.append(f"Elementos nulos bajo la diagonal: "
                   f"{ceros_bajo_diagonal(res.U)} de {n * (n - 1) // 2}")
        out.append("")
    ruta.write_text("\n".join(out), encoding="utf-8")


# Punto de entrada: este bloque se ejecuta solo al iniciar el archivo
# directamente. Organiza las etapas del proyecto y guarda los resultados.
if __name__ == "__main__":
    SALIDA.mkdir(exist_ok=True)

    # LECTURA DE DATOS: conserva textos exactos y crea copias float para Gauss.
    A_str = leer_matriz_csv(DATOS / "matriz_a.csv")
    b_str = leer_vector_csv(DATOS / "vector_b.csv")
    A, b = a_float(A_str), a_float(b_str)
    # El trabajo está configurado para un sistema cuadrado de 12 ecuaciones.
    n = len(A)
    assert n == 12 and all(len(f) == n for f in A) and len(b) == n

    # REFERENCIA EXACTA: permite evaluar error y dígitos de las soluciones.
    x_exacto = ver.solucion_exacta(A_str, b_str)

    # T1 / T2: ejecuta Gauss con y sin pivoteo para comparar ambos métodos.
    c_sin, c_con = Contador(n), Contador(n)
    r_sin = eliminar(A, b, "ninguno", c_sin)
    r_con = eliminar(A, b, "parcial_escalado", c_con)

    # Evalúa el error relativo y los dígitos correctos de cada corrida
    # para comprobar el efecto numérico del pivoteo.
    err_sin = ver.error_relativo_inf(r_sin.x, x_exacto)