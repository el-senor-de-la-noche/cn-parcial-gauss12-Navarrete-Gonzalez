# Integrantes: Elias Gonzalez, Daniel Navarrete
# Ramo: Computacion Numerica
# Profesor: Juan Cabezas


# ARCHIVO: CONTEO DE OPERACIONES
# Este módulo registra el costo de la eliminación gaussiana y genera un
# informe para comparar las operaciones realizadas con una referencia teórica
# Un flop representa una suma, resta, multiplicación o división
# Las comparaciones y los intercambios de filas se miden por separado
# como sobrecosto asociado a la selección del pivote
# Los flops se cuentan cuando gauss.py ejecuta cada operación, no mediante
# una estimación posterior


# Contador de operaciones separado por fases del algoritmo.
# Cada grupo guarda [multiplicaciones/divisiones, sumas/restas]
class Contador:
    # Prepara un registro por cada paso de eliminación (n-1 pasos), más
    # registros independientes para las sustituciones del sistema
    def __init__(self, n):
        self.n = n
        # En cada paso se guardan dos conteos: multiplicaciones/divisiones y sumas/restas.
        self.elim = [[0, 0] for _ in range(n - 1)]
        self.regresiva = [0, 0]    # sustitución regresiva
        self.progresiva = [0, 0]   # sustitución progresiva (solo en refinamiento)
        # Métricas del pivoteo, reportadas al margen de los flops principales.
        self.pivote = {
            "comparaciones_escala": 0,   # cálculo de s_i = max|a_ij| por fila
            "divisiones_razon": 0,       # razones |a_ik| / s_i
            "comparaciones_busqueda": 0, # búsqueda del mayor |a_ik| / s_i
            "intercambios": 0,           # intercambios de filas realizados
        }
        self._destino = None

    # selección de fase
    # Envía las próximas operaciones al paso de eliminación indicado
    def fase_eliminacion(self, paso):
        self._destino = self.elim[paso]

    # Activa el registro usado al despejar las incógnitas de abajo arriba
    def fase_regresiva(self):
        self._destino = self.regresiva

    # Activa el registro usado para resolver con la matriz L
    def fase_progresiva(self):
        self._destino = self.progresiva

    # operaciones aritméticas (flops)
    # Suma k multiplicaciones al contador de la fase activa
    def mul(self, k=1):
        self._destino[0] += k

    # Suma k divisiones al mismo grupo de multiplicaciones/divisiones
    def div(self, k=1):
        self._destino[0] += k

    # Suma k adiciones o sustracciones a la fase activa
    def sumres(self, k=1):
        self._destino[1] += k

    # sobrecosto de pivoteo
    # Cuenta comparaciones para obtener la escala máxima de cada fila
    def pivote_escala(self, k=1):
        self.pivote["comparaciones_escala"] += k

    # Cuenta las divisiones necesarias para comparar pivotes normalizados
    def pivote_razon(self, k=1):
        self.pivote["divisiones_razon"] += k

    # Cuenta las comparaciones realizadas al buscar el mejor pivote
    def pivote_busqueda(self, k=1):
        self.pivote["comparaciones_busqueda"] += k

    # Registra los intercambios efectivos de filas
    def pivote_intercambio(self, k=1):
        self.pivote["intercambios"] += k

    # totales
    # Suma las operaciones de todos los pasos de eliminación
    def total_eliminacion(self):
        return sum(m + s for m, s in self.elim)

    # Devuelve el total de la sustitución regresiva
    def total_regresiva(self):
        return sum(self.regresiva)

    # Devuelve el total registrado en la sustitución progresiva
    def total_progresiva(self):
        return sum(self.progresiva)

    def total_base(self):
        # Totaliza el trabajo principal de Gauss: eliminación de A y actualización
        # de b, más las operaciones de la sustitución regresiva
        return self.total_eliminacion() + self.total_regresiva()

    # Agrupa comparaciones y divisiones del pivoteo; los intercambios
    # de filas se informan aparte y no se suman a este total
    def total_sobrecosto_pivoteo(self):
        p = self.pivote
        return (p["comparaciones_escala"] + p["divisiones_razon"]
                + p["comparaciones_busqueda"])


def cota_referencia(n):
    # Calcula una cota teórica de referencia para una matriz densa n x n
    # En cada etapa quedan m filas por eliminar: se cuentan m divisiones para
    # los multiplicadores, 2*m*(m+1) operaciones sobre A y 2*m sobre b
    # Después se agregan n*(n-1)+n operaciones de sustitución regresiva
    # Esta fórmula sirve para comparar; no reemplaza el conteo instrumentado
    elim = sum(m + 2 * m * (m + 1) + 2 * m for m in range(1, n))
    regresiva = n * (n - 1) + n
    return elim + regresiva


# Verifica si el conteo real difiere como máximo un 10 % de la cota
# teórica, salvo que se proporcione otra tolerancia
def dentro_de_tolerancia(total, referencia, tolerancia=0.10):
    return abs(total - referencia) <= tolerancia * referencia


def reporte_markdown(c_sin, c_con):
    # Construye el informe Markdown con el desglose por etapa, la comparación
    # entre corridas, el costo del pivoteo y el refinamiento iterativo
    # La corrida con pivoteo determina el tamaño del sistema y la cota
    # empleada para contrastar ambos conteos
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
    # Se agrega una fila de la tabla por cada paso de eliminación, separando
    # multiplicaciones/divisiones de sumas/restas
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
    # Contrasta ambos métodos contra la misma referencia; una desviación
    # positiva significa que se contaron más operaciones que la cota
    for nombre, c in (("Sin pivoteo", c_sin), ("Pivoteo parcial escalado", c_con)):
        t = c.total_base()
        dev = 100.0 * (t - ref) / ref
        ok = "sí" if dentro_de_tolerancia(t, ref) else "no"
        L.append(f"| {nombre} | {t} | {ref} | {dev:+.2f} % | {ok} |")
    L.append("")

    # Las estadísticas de estabilidad se muestran por separado para no
    # mezclarlas con el costo aritmético principal
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
    # Se estima una pasada de corrección con L y U ya disponibles:
    # sustitución progresiva + regresiva + suma de la corrección a x
    nprog = n * (n - 1)
    nreg = n * (n - 1) + n
    L.append(f"- Sustitución progresiva: {nprog} flops")
    L.append(f"- Sustitución regresiva: {nreg} flops")
    L.append(f"- Actualización x + d: {n} flops")
    L.append(f"- **Total por pasada: {nprog + nreg + n} flops** "
             "(el residuo en precisión extendida se calcula aparte, con "
             "aritmética exacta).\n")
    return "\n".join(L)
