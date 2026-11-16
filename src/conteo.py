#Integrantes: Elias Gonzalez, Daniel Navarrete
#Ramo: Computacion Numerica
#Profesor: Juan Cabezas


"""Instrumentación y reporte del costo operacional.

Convención de conteo: 1 flop = 1 operación aritmética elemental
(multiplicación, división, suma o resta). Las comparaciones, los valores
absolutos y los intercambios de filas NO se cuentan como flops: se reportan
aparte como sobrecosto de estabilidad (pivoteo).

El contador se incrementa dentro de los bucles de gauss.py, operación por
operación: los totales salen del código, no de una fórmula escrita a mano.
"""


class Contador:
    def __init__(self, n):
        self.n = n
        # Eliminación: un par [mul_div, sum_res] por cada uno de los n-1 pasos.
        self.elim = [[0, 0] for _ in range(n - 1)]
        self.regresiva = [0, 0]    # sustitución regresiva
        self.progresiva = [0, 0]   # sustitución progresiva (solo en refinamiento)
        # Sobrecosto de estabilidad (no son flops).
        self.pivote = {
            "comparaciones_escala": 0,   # cálculo de s_i = max|a_ij| por fila
            "divisiones_razon": 0,       # razones |a_ik| / s_i
            "comparaciones_busqueda": 0, # búsqueda del mayor |a_ik| / s_i
            "intercambios": 0,           # intercambios de filas realizados
        }
        self._destino = None

    # --- selección de fase -------------------------------------------------
    def fase_eliminacion(self, paso):
        self._destino = self.elim[paso]

    def fase_regresiva(self):
        self._destino = self.regresiva

    def fase_progresiva(self):
        self._destino = self.progresiva

    # --- operaciones aritméticas (flops) -----------------------------------
    def mul(self, k=1):
        self._destino[0] += k

    def div(self, k=1):
        self._destino[0] += k

    def sumres(self, k=1):
        self._destino[1] += k

    # --- sobrecosto de pivoteo ---------------------------------------------
    def pivote_escala(self, k=1):
        self.pivote["comparaciones_escala"] += k

    def pivote_razon(self, k=1):
        self.pivote["divisiones_razon"] += k

    def pivote_busqueda(self, k=1):
        self.pivote["comparaciones_busqueda"] += k

    def pivote_intercambio(self, k=1):
        self.pivote["intercambios"] += k

    # --- totales -----------------------------------------------------------
    def total_eliminacion(self):
        return sum(m + s for m, s in self.elim)

    def total_regresiva(self):
        return sum(self.regresiva)

    def total_progresiva(self):
        return sum(self.progresiva)

    def total_base(self):
        """Flops de eliminación (incluye actualización de b) + regresiva."""
        return self.total_eliminacion() + self.total_regresiva()

    def total_sobrecosto_pivoteo(self):
        p = self.pivote
        return (p["comparaciones_escala"] + p["divisiones_razon"]
                + p["comparaciones_busqueda"])


def cota_referencia(n):
    """Mínimo teórico de la eliminación densa (para comparar, no para contar).

    Paso k (m = n-k filas por debajo): m divisiones (multiplicadores),
    2*m*(m+1) para actualizar las columnas k..n de A (barrido desde la
    columna k, como exige el enunciado) y 2*m para b. Más la sustitución
    regresiva: n(n-1) + n.
    """
    elim = sum(m + 2 * m * (m + 1) + 2 * m for m in range(1, n))
    regresiva = n * (n - 1) + n
    return elim + regresiva


def dentro_de_tolerancia(total, referencia, tolerancia=0.10):
    return abs(total - referencia) <= tolerancia * referencia


def reporte_markdown(c_sin, c_con):
    """Tabla de operaciones por paso (11 pasos) para ambas corridas."""
    n = c_con.n
    ref = cota_referencia(n)
    L = []
    L.append("# Conteo de operaciones (flops)\n")
    L.append("Convención: 1 flop = 1 operación aritmética elemental "
             "(multiplicación, división, suma o resta). Comparaciones, valores "
             "absolutos e intercambios de filas se reportan aparte como "
             "sobrecosto de estabilidad. El barrido del paso k va de la columna "
             "k en adelante; el multiplicador m(i,k) se calcula una sola vez por "
             "fila y se reutiliza. Los totales provienen del contador "
             "instrumentado en `src/gauss.py`.\n")

    L.append("## Eliminación por paso (corrida con pivoteo parcial escalado)\n")
    L.append("| Paso k | Mult/Div | Sumas/Restas | Total |")
    L.append("| --- | ---: | ---: | ---: |")
    for k, (m, s) in enumerate(c_con.elim, start=1):
        L.append(f"| {k} | {m} | {s} | {m + s} |")
    em = sum(m for m, _ in c_con.elim)
    es = sum(s for _, s in c_con.elim)
    L.append(f"| **Eliminación** | **{em}** | **{es}** | **{em + es}** |")
    L.append(f"| Sustitución regresiva | {c_con.regresiva[0]} | "
             f"{c_con.regresiva[1]} | {sum(c_con.regresiva)} |")
    L.append(f"| **Total base** | **{em + c_con.regresiva[0]}** | "
             f"**{es + c_con.regresiva[1]}** | **{c_con.total_base()}** |\n")

    L.append("## Comparación con la cota\n")
    L.append("| Corrida | Flops contados | Cota de referencia | Desviación | "
             "¿Dentro de ±10 %? |")
    L.append("| --- | ---: | ---: | ---: | --- |")
    for nombre, c in (("Sin pivoteo", c_sin), ("Pivoteo parcial escalado", c_con)):
        t = c.total_base()
        dev = 100.0 * (t - ref) / ref
        ok = "sí" if dentro_de_tolerancia(t, ref) else "no"
        L.append(f"| {nombre} | {t} | {ref} | {dev:+.2f} % | {ok} |")
    L.append("")

    p = c_con.pivote
    L.append("## Sobrecosto de estabilidad (pivoteo parcial escalado)\n")
    L.append("Se reporta aparte; no está incluido en el total base.\n")
    L.append("| Concepto | Cantidad |")
    L.append("| --- | ---: |")
    L.append(f"| Comparaciones para calcular las escalas s_i | "
             f"{p['comparaciones_escala']} |")
    L.append(f"| Divisiones |a_ik| / s_i | {p['divisiones_razon']} |")
    L.append(f"| Comparaciones de la búsqueda de pivote | "
             f"{p['comparaciones_busqueda']} |")
    L.append(f"| Intercambios de filas realizados | {p['intercambios']} |")
    L.append(f"| **Sobrecosto total (operaciones de pivoteo)** | "
             f"**{c_con.total_sobrecosto_pivoteo()}** |\n")

    L.append("## Costo del refinamiento iterativo (por pasada)\n")
    L.append("Reutiliza L y U (no se vuelve a factorizar): una sustitución "
             "progresiva y una regresiva por pasada.\n")
    nprog = n * (n - 1)
    nreg = n * (n - 1) + n
    L.append(f"- Sustitución progresiva: {nprog} flops")
    L.append(f"- Sustitución regresiva: {nreg} flops")
    L.append(f"- Actualización x + d: {n} flops")
    L.append(f"- **Total por pasada: {nprog + nreg + n} flops** "
             "(el residuo en precisión extendida se calcula aparte, con "
             "aritmética exacta).\n")
    return "\n".join(L)
