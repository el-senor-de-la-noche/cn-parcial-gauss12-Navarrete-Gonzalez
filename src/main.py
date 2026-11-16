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

# Agrega la raíz del proyecto al path de búsqueda de módulos.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

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
    err_con = ver.error_relativo_inf(r_con.x, x_exacto)
    dig_sin = min(ver.digitos_correctos(v, e) for v, e in zip(r_sin.x, x_exacto))
    dig_con = min(ver.digitos_correctos(v, e) for v, e in zip(r_con.x, x_exacto))

    # T4: refina tres veces la solución obtenida con pivoteo.
    historia = ver.refinar(r_con, A_str, b_str, r_con.x, x_exacto,
                           PASADAS_REFINAMIENTO)
    # La última entrada del historial contiene la solución refinada.
    x_final = historia[-1]["x"]

    # T3: exporta la solución sin afirmar más precisión de la medida.
    filas_sol = []
    # Para cada incógnita se calcula cuántos dígitos se pueden afirmar,
    # y se conserva solo esa precisión al escribir solucion.csv.
    for i, (xi, xe) in enumerate(zip(x_final, x_exacto), start=1):
        d = ver.digitos_correctos(xi, xe)
        filas_sol.append((i, xi, d))
    with open(SALIDA / "solucion.csv", "w", encoding="utf-8", newline="") as f:
        f.write("i,x,digitos_correctos_afirmados\n")
        for i, xi, d in filas_sol:
            f.write(f"{i},{xi:.{max(d, 1) - 1}e},{d}\n")

    # TRAZAS Y CONTEO: registra los pasos y las operaciones utilizadas.
    escribir_traza(SALIDA / "traza_eliminacion.txt", A, b,
                   [("pivoteo parcial escalado (corrida principal)", r_con),
                    ("sin pivoteo", r_sin)])

    # Se prepara el informe de operaciones y se le agregan indicadores
    # para comparar el comportamiento de ambos métodos.
    md = reporte_markdown(c_sin, c_con)
    md += "\n## Indicadores de estabilidad (T2)\n\n"
    md += ("| Indicador | Sin pivoteo | Pivoteo parcial escalado |\n"
           "| --- | ---: | ---: |\n")
    md += (f"| Multiplicador máximo | {r_sin.max_multiplicador:.4e} | "
           f"{r_con.max_multiplicador:.4e} |\n")
    md += (f"| Factor de crecimiento | {r_sin.factor_crecimiento:.4e} | "
           f"{r_con.factor_crecimiento:.4e} |\n")
    md += (f"| log2(multiplicador máximo) | "
           f"{math.log2(r_sin.max_multiplicador):.1f} | "
           f"{math.log2(r_con.max_multiplicador):.1f} |\n")
    md += (f"| Bits de mantisa desplazados (máx., de 53) | "
           f"{r_sin.bits_perdidos_max:.1f} | {r_con.bits_perdidos_max:.1f} |\n")
    md += (f"| Intercambios de filas | {r_sin.intercambios} | "
           f"{r_con.intercambios} |\n")
    md += (f"| Error relativo final (norma inf, vs. exacto) | {err_sin:.3e} | "
           f"{err_con:.3e} |\n")
    md += (f"| Dígitos correctos (mínimo entre las 12 incógnitas) | "
           f"{dig_sin} | {dig_con} |\n")
    (SALIDA / "conteo_flops.md").write_text(md, encoding="utf-8")

    # VERIFICACIÓN: analiza residuos, condicionamiento y referencia NumPy.
    cond = ver.condicionamiento(A)
    ctrl = ver.control_numpy(A, b, x_exacto)
    r0 = ver.residuo_exacto(A_str, b_str, r_sin.x)

    # Se prepara un informe de residuos, historial de refinamiento,
    # condicionamiento y comparación contra NumPy.
    L = []
    L.append("RESIDUO r = b - A x (precisión extendida: aritmética racional "
             "exacta con fractions)\n")
    L.append(f"Residuo de la corrida SIN pivoteo:  ||r||_inf = "
             f"{ver.formatear(ver.norma_inf(r0), 6)}")
    L.append("")
    L.append("Refinamiento iterativo sobre la corrida con pivoteo parcial "
             "escalado")
    L.append("(residuo exacto, corrección d = U^-1 L^-1 P r con L y U "
             "reutilizados):\n")
    L.append(f"{'pasada':>6} {'||r||_inf':>14} {'error relativo':>16} "
             f"{'dígitos correctos (mín.)':>26}")
    for h in historia:
        L.append(f"{h['pasada']:>6} {ver.formatear(h['norma_r'], 6):>14} "
                 f"{h['error_rel']:>16.3e} {h['min_digitos']:>26}")
    # Calcula el factor de mejora entre el error inicial y el final;
    # el mínimo del denominador evita una división entre cero.
    mejora = historia[0]["error_rel"] / max(historia[-1]["error_rel"], 1e-300)
    L.append(f"\nMejora del error relativo con el refinamiento: "
             f"x{mejora:.3g} (la norma del residuo puede no bajar: queda "
             f"limitada por el redondeo de x a doble precisión)")
    L.append("")
    # Recalcula el residuo de la solución final usando los datos exactos.
    r_fin = ver.residuo_exacto(A_str, b_str, x_final)
    L.append("Vector residuo final r (solución reportada):")
    for i, ri in enumerate(r_fin, start=1):
        L.append(f"  r[{i:2d}] = {ver.formatear(ri, 8)}")
    L.append(f"\n||r||_inf (final) = {ver.formatear(ver.norma_inf(r_fin), 8)}")
    L.append("")
    L.append("CONDICIONAMIENTO")
    L.append(f"  kappa_2(A)                 = {cond['kappa2']:.4e}")
    L.append(f"  cota kappa*eps (eps=2.22e-16) = {cond['cota']:.3e}")
    L.append(f"  dígitos que garantiza la cota = {cond['digitos_cota']:.2f}")
    L.append("")
    L.append("CONTROL CON BIBLIOTECA (numpy.linalg.solve, solo como control)")
    L.append(f"  error relativo = {ctrl['error_rel']:.3e}   "
             f"dígitos correctos (mín.) = {ctrl['min_digitos']}")
    (SALIDA / "residuo.txt").write_text("\n".join(L) + "\n", encoding="utf-8")

    # RESUMEN: muestra los indicadores principales al terminar.
    ref = cota_referencia(n)
    print(f"Python {platform.python_version()} | numpy {numpy.__version__} | "
          f"mpmath {mpmath.__version__}")
    print(f"Triangularización: {ceros_bajo_diagonal(r_con.U)}/66 ceros "
          f"bajo la diagonal (pivoteo) y {ceros_bajo_diagonal(r_sin.U)}/66 "
          f"(sin pivoteo)")
    print(f"Flops base: {c_con.total_base()} (cota {ref}, "
          f"{'dentro' if dentro_de_tolerancia(c_con.total_base(), ref) else 'fuera'} "
          f"de ±10 %) | sobrecosto de pivoteo: "
          f"{c_con.total_sobrecosto_pivoteo()} operaciones")
    print(f"Sin pivoteo : mult. máx {r_sin.max_multiplicador:.3e}, "
          f"crecimiento {r_sin.factor_crecimiento:.3e}, "
          f"error {err_sin:.3e}, dígitos {dig_sin}")
    print(f"Con pivoteo : mult. máx {r_con.max_multiplicador:.3e}, "
          f"crecimiento {r_con.factor_crecimiento:.3e}, "
          f"error {err_con:.3e}, dígitos {dig_con}")
    print(f"kappa_2(A) = {cond['kappa2']:.3e} | cota kappa*eps = "
          f"{cond['cota']:.2e}")
    print(f"||r||_inf final = {ver.formatear(ver.norma_inf(r_fin), 6)}")
    print(f"Resultados escritos en: {SALIDA}")

