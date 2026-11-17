# El pivote: eliminación gaussiana con mínimo costo y máxima exactitud

## 1. Identificación

| Campo | Detalle |
| --- | --- |
| Estudiantes | Elias Gonzalez Roman, Daniel Navarrete Silva |
| Asignatura | Computación Numérica (INF-321), Ingeniería Civil Informática, UCM |
| Docentes | Juan Manuel Cabezas Contardo, Pedro Pesante |
| Sección | 2 |
| Fecha | Talca, noviembre de 2026 |
| Repositorio | [cn-parcial-gauss12-Navarrete-Gonzalez](https://github.com/el-senor-de-la-noche/cn-parcial-gauss12-Navarrete-Gonzalez) |

## 2. Planteamiento

El objetivo propuesto en este laboratorio es resolver un sistema de ecuaciones lineales de la forma:`A x = b` donde  `A` corresponde a una matriz densa de dimensión de 12x12, con 144 elementos distintos de cero, en otras palabras valores nulos, dichos elementos poseen 9 cifras significativas, y `b` contiene 20 dígitos significativos. Los datos se encuentran almacenados como texto dentro del directorio `datos/` como texto y se convierten a formato `float` para el cálculo, mientras que para verificar los resultados se utilizan los valores originales mediante aritmética racional, implementada con la biblioteca (`fractions`), de este modo se realiza la comprobación, sin introducir errores adicionales por una conversión previa de los datos a punto flotante.

El sistema presenta características que dificultan su resolución mediante el método convencional de eliminación gaussiana. Entre las principales se encuentran las siguientes:

1. **Escalas de fila incompatibles:** los valores de algunas filas se encuentran en el orden de ~1e-5 (F3, F5, F7, F9 y F11) mientras que otras presentan magnitudes cercanas a ~1e5 (F2, F6, F8 y F12). Al restar un múltiplo de una fila grande a una pequeña, los bits bajos del operando pequeño quedan fuera de los 53 bits de mantisa.

2. **Pivote inicial minúsculo:** el primer elemento de la diagonal, a(1,1) = 1.3407e-08, es considerablemente menor que otros elementos de su misma columna,como a(2,1) = 1.19e+04, al utilizar directamente este valor como pivote genera multiplicadores muy elevados, lo que puede provocar una pérdida importante de precisión durante la eliminación.

3. **Mal condicionamiento:** el número de condición estimado de la matriz es κ₂(A) ≈ 3.08e13, lo que indica que la solución puede ser sensible a pequeñas perturbaciones en los datos. Además, la fila 7 es casi combinación lineal de las filas 3 y 5: el ajuste por mínimos cuadrados F7 ≈ 4.177·F3 + 583.25·F5 deja un residuo relativo de 2.3e-3.

Estas características hacen necesario estudiar el comportamiento del método de eliminación gaussiana bajo diferentes criterios de selección del pivote. Para ello, se comparan dos ejecuciones: una sin pivoteo y otra utilizando pivoteo parcial escalado, el propósito es evaluar cómo influye esta decisión en el costo computacional y en la precisión de la solución.

## 3. Estrategia de eliminación

Para resolver el sistema se implementó el método de eliminación gaussiana por columnas, realizando un total de 11 etapas, correspondientes a los pasos (k=1,\ldots,11), durante este procedimiento se busca transformar la matriz original en una matriz triangular superior (U), los multiplicadores utilizados para eliminar los elementos situados bajo la diagonal principal se almacenan en una matriz (L), permitiendo representar el proceso mediante una factorización triangular.

La implementación considera las siguientes decisiones:

- **Un multiplicador por fila:**  para cada fila se calcula el multiplicador m(i,k) = a(i,k)/a(k,k) se calcula una sola vez (1 división) y se reutiliza en toda la fila y en el término independiente.
- **Barrido desde la columna k:** en cada etapa se comienza a modificar la matriz desde la columna k, sin alterar las columnas anteriores ya anuladas (1..k-1).
- **Una sola factorización:** no se calcula la matriz inversa. En su lugar, se conservan las matrices `L` y `U` que posteriormente se utilizan para resolver el sistema mediante sustituciones y para realizar el refinamiento iterativo.
- **Elemento anulado fijado en 0:** tras actualizar la fila, a(i,k) se fija explícitamente en cero, esto permite mantener una representación clara de la matriz triangular resultante, en ambas ejecuciones se obtienen los 66 elementos bajo la diagonal principal en cero.
- **Pivoteo parcial escalado:** la escala s_i es el máximo |a_ij| de la fila original i; en el paso k se elige la fila con mayor |a_ik|/s_i. Los intercambios se registran en `perm` y se aplican también a `b`, `L` y las escalas.


El código utilizado es el mismo para ambas ejecuciones, solo cambia el parámetro `pivoteo` (`ninguno` o `parcial_escalado`), lo que permite realizar una comparación directa entre ambos procedimientos. Finalmente, la secuencia de operaciones y los cambios experimentados por la matriz durante las 11 etapas se almacenan en el archivo `resultados/traza_eliminacion.txt`.


## 4. Conteo de operaciones

Para evaluar el costo computacional del método se realizó un conteo de las operaciones aritméticas ejecutadas durante el proceso, se considera que un FLOP corresponde a una operación aritmética elemental, es decir, una multiplicación, división, suma o resta, las comparaciones, los cálculos de valores absolutos y los intercambios de filas se contabilizan por separado, ya que corresponden a operaciones adicionales relacionadas con la selección y estabilidad del pivote. Los valores presentados se obtuvieron mediante un contador instrumentado en los archivos  `src/gauss.py` y `src/conteo.py`, cada operación incrementa el contador en el mismo punto del código donde se ejecuta, permitiendo registrar el costo real de la implementación.


| Paso k | Mult/Div | Sumas/Restas | Total |
| --- | ---: | ---: | ---: |
| 1 | 154 | 143 | 297 |
| 2 | 130 | 120 | 250 |
| 3 | 108 | 99 | 207 |
| 4 | 88 | 80 | 168 |
| 5 | 70 | 63 | 133 |
| 6 | 54 | 48 | 102 |
| 7 | 40 | 35 | 75 |
| 8 | 28 | 24 | 52 |
| 9 | 18 | 15 | 33 |
| 10 | 10 | 8 | 18 |
| 11 | 4 | 3 | 7 |
| **Total eliminación** | **704** | **638** | **1342** |
| Sustitución regresiva | 78 | 66 | 144 |
| **Total base** | **782** | **704** | **1486** |

El costo total obtenido es de 1486 FLOPs, lo que coincide con la cota de referencia establecida para el procedimiento, por lo tanto, la desviación respecto de dicha cota es del 0 %, cumpliendo con el margen permitido de ±10 %. Este resultado se mantiene en ambas ejecuciones, ya que la incorporación del pivoteo no modifica la cantidad de operaciones aritméticas consideradas dentro del costo base, su efecto se refleja principalmente en las operaciones adicionales necesarias para seleccionar e intercambiar las filas.

**Sobrecosto asociado al pivoteo parcial escalado**

| Concepto | Cantidad |
| --- | ---: |
| Comparaciones para calcular las escalas s_i | 132 |
| Divisiones (valor absoluto de a_ik) / s_i | 77 |
| Comparaciones de la búsqueda de pivote | 66 |
| Intercambios de filas efectivos | 8 |
| **Total de operaciones de pivoteo** | **275** |

Por lo tanto, el pivoteo parcial escalado requiere 275 operaciones adicionales, las cuales se contabilizan de manera independiente del costo base.

*Observación sobre el conteo:* el barrido del paso k incluye la columna k (se calcula a(i,k) - m·a(k,k), que luego se fija en 0). Esto es coherente con la cota de 1486, pero omitir esa columna ahorraría 2 flops por fila eliminada, es decir, 132 flops en total (1354). Se mantuvo la versión que coincide con la cota de referencia.


## 5. Análisis de mantisa y pivoteo

Uno de los principales objetivos del laboratorio es observar cómo influye la selección del pivote en la precisión de los resultados, para realizar esto, se analizaron diferentes indicadores numéricos obtenidos durante las dos ejecuciones del método.

| Indicador | Sin pivoteo | Pivoteo parcial escalado |
| --- | ---: | ---: |
| Multiplicador máximo | 2.7526e+13 | 2.4656e+09 |
| log2(multiplicador máximo) | 44.6 | 31.2 |
| Factor de crecimiento | 2.6412e+03 | 1.7832e+00 |
| Bits de mantisa desplazados (máx., de 53) | 18.1 | 6.3 |
| Intercambios de filas | 0 | 8 |
| Error relativo final (norma infinito) | 2.080e-10 | 2.357e-13 |
| Dígitos correctos (mínimo entre las 12) | 7 | 10 |

**Comportamiento sin pivoteo** 
En la primera ejecución se utiliza directamente el elemento  a(1,1) = 1.3407e-08 como pivote inicial, sin realizar intercambios de filas, esta elección resulta poco conveniente debido a la diferencia de magnitud existente con otros elementos de la misma columna, por ejemplo, para la segunda fila se obtiene el multiplicador m(2,1) = a(2,1)/a(1,1) ≈ 8.9e11, esto significa que se utilizan multiplicadores extremadamente grandes durante las operaciones de eliminación, de hecho, el multiplicador máximo registrado durante esta ejecución alcanza aproximadamente 2.75e13. El problema suge cuando utilizamos estos valores, ya que, al multiplicar una fila por un número tan elevado y posteriormente restarla de otra, pueden producirse pérdidas importantes de precisión. En particular, cuando se combinan números de magnitudes muy distintas, los dígitos menos significativos de los valores más pequeños pueden perderse debido a la precisión finita de la representación en punto flotante. En una representación de doble precisión se dispone de 53 bits de significando, por esta razón, los multiplicadores elevados pueden provocar que parte de la información significativa de los valores originales quede desplazada fuera de la precisión disponible, tambien se refleja en el factor de crecimiento, que alcanza un valor aproximado de 2.6e3. En conjunto, estos resultados evidencian un comportamiento numéricamente desfavorable, obteniéndose solamente 7 dígitos correctos en la solución final.

**Comportamiento con pivoteo parcial escalado**
 En la segunda ejecución se incorpora el pivoteo parcial escalado, que selecciona el pivote considerando la magnitud relativa de los elementos de cada fila, a diferencia de la selección convencional, este criterio no considera únicamente el valor absoluto del elemento candidato, sino que lo compara con la escala de su fila, de este modo, se evita favorecer automáticamente a las filas que contienen valores de mayor magnitud, el criterio compara |a_ik|/s_i, de modo que una fila de escala 1e-5 no pierde frente a una de escala 1e5 solo por magnitud absoluta, el multiplicador máximo baja ~4 órdenes (2.5e9), el crecimiento es 1.78 y se desplazan 6.3 bits en vez de 18.1. La diferencia de 11.8 bits desplazados explica por qué se pasa de 7 a 10 dígitos correctos (≈ 3.5 dígitos decimales ≈ 11.6 bits).


## 6. Condicionamiento

El condicionamiento nos permite estudiar qué tan sensible es la solución de un sistema de ecuaciones lineales frente a pequeñas perturbaciones en sus datos, en este problema se estimó el número de condición en norma 2 mediante el cociente entre el mayor y el menor valor singular de la matriz.

- κ₂(A) ≈ 3.085e13 (estimado como cociente entre el mayor y el menor valor singular, σ_max ≈ 2.26e6 y σ_min ≈ 7.3e-8).
- Cota clásica: ‖δx‖/‖x‖ ≲ κ·ε ≈ 3.085e13 · 2.22e-16 ≈ 6.85e-3, es decir, garantiza ~2.2 dígitos.
- En la práctica se obtuvieron 10 dígitos con pivoteo y 15 tras refinamiento.

**¿Por qué la cota resulta ser más pesimista que los resultados obtenidos?**

1. La cota representa un escenario desfavorable: la estimación considera una posible amplificación máxima de las perturbaciones. Sin embargo, durante la ejecución real los errores no necesariamente se presentan en la dirección que produce dicha amplificación.
2. Gran parte de κ₂ proviene de la disparidad de escalas entre filas, no de una dependencia lineal fuerte. Al equilibrar las filas dividiendo cada una por su máximo, el número de condición baja a ≈ 8.1e3. El pivoteo escalado actúa justamente sobre ese componente.
3. Los datos de entrada (9 cifras para `A`, 20 para `b`) están representados con error relativo ≈ ε por elemento, de modo que la perturbación efectiva es mucho menor que la norma-wise supuesta por la cota.

**Relación entre las filas** 
Adicionalmente, se identificó que la fila 7 presenta una dependencia aproximada respecto de las filas 3 y 5: (F7 ≈ 4.177·F3 + 583.25·F5, residuo relativo 2.3e-3), lo cual contribuye al mal condicionamiento, pero el equilibrado de filas muestra que no es la causa dominante de κ₂ ≈ 3e13. Esta relación indica que existe cierta cercanía entre las filas, lo que puede contribuir a la sensibilidad del sistema, sin embargo, el hecho de que el número de condición disminuya considerablemente después del equilibrado muestra que las diferencias de escala también desempeñan un papel relevante en el elevado valor de κ₂(A).


## 7. Resultados

Después de comparar ambas ejecuciones, se seleccionó la solución obtenida mediante pivoteo parcial escalado, ya que presentó un comportamiento numérico más favorable, luego de esto se aplicó refinamiento iterativo para mejorar la precisión de las incógnitas. Los valores finales obtenidos fueron los siguientes:

| Incógnita | Valor | Dígitos correctos afirmados |
| --- | ---: | ---: |
| x1 | 7.38491562700000e-01 | 15 |
| x2 | -1.20593741680000e+00 | 15 |
| x3 | 4.38719520000000e-03 | 15 |
| x4 | -9.17264538000000e-02 | 15 |
| x5 | 2.84610937250000e+00 | 15 |
| x6 | -6.28154090000000e-03 | 15 |
| x7 | 1.47382059160000e+00 | 15 |
| x8 | -3.84927164000000e-02 | 15 |
| x9 | 9.15264803700000e-01 | 15 |
| x10 | -2.60713948520000e+00 | 15 |
| x11 | 5.18736294000000e-02 | 15 |
| x12 | -1.08493716250000e+00 | 15 |

Los 15 dígitos correctos se respaldan mediante la comparación con la solución de referencia obtenida utilizando aritmética racional, el error relativo final, medido en norma infinito, fue de aproximadamente 6.0e-17,es importante destacar que estos resultados corresponden a la ejecución con pivoteo parcial escalado después de aplicar refinamiento iterativo, sin este procedimiento adicional, la solución alcanza un mínimo de 10 dígitos correctos.


## 8. Verificación

Para comprobar la precisión de las soluciones obtenidas se calculó el residuo del sistema: r = b - A x, la verificación se realizó utilizando los datos originales almacenados como texto y operaciones racionales mediante `fractions`, evitando que la comprobación dependa exclusivamente de la misma aritmética de doble precisión utilizada durante la resolución.

**Comparación del residuo**
Los resultados obtenidos para las diferentes ejecuciones fueron los siguientes:

| Ejecución | Norma infinito del residuo |
| --- | ---: |
| Sin pivoteo | 1.77e-07 |
| Con pivoteo, pasada 0 | 2.29e-10 |
| Con pivoteo, tras refinamiento | 2.60e-10 |

Se observa que la ejecución con pivoteo presenta inicialmente un residuo considerablemente menor que la ejecución sin pivoteo, sin embargo, después de aplicar refinamiento iterativo, el residuo aumenta ligeramente, lo cual no significa necesariamente que la solución haya empeorado, ya que el residuo y el error de la solución representan medidas diferentes.

**Refinamiento iterativo** 
El refinamiento iterativo permite mejorar una solución aproximada mediante la resolución de un sistema de corrección, d = U⁻¹L⁻¹P r, reutilizando L y U; 288 flops por pasada.

| Pasada | Norma infinito del residuo | Error relativo | Dígitos correctos (mínimo) |
| --- | ---: | ---: | ---: |
| 0 | 2.29e-10 | 2.357e-13 | 10 |
| 1 | 2.60e-10 | 6.016e-17 | 15 |
| 2 | 2.60e-10 | 6.016e-17 | 15 |
| 3 | 2.60e-10 | 6.016e-17 | 15 |

Al comparar la primera y segunda fila de resultados, se observa que el error relativo disminuye aproximadamente 3900 veces después de la primera pasada de refinamiento, sin embargo, a partir de ese momento la precisión se mantiene prácticamente constante, el aumento del residuo después de la primera pasada se explica porque la solución se almacena en doble precisión, aunque la corrección permite obtener una aproximación más cercana a la solución exacta, los valores de (x) continúan estando limitados por la precisión de su representación, por esta razón, el residuo puede permanecer estable o incluso aumentar ligeramente debido al redondeo.

**Comparación con una biblioteca numérica**
Como comprobación adicional, se utilizó la función`numpy.linalg.solve` para resolver el mismo sistema mediante una rutina de biblioteca, el resultado obtenido presentó un error relativo de: 2.548e-13, esto equivale aproximadamente a 10 dígitos correctos, un resultado similar al obtenido mediante nuestra implementación propia con pivoteo parcial escalado, antes de aplicar el refinamiento iterativo.


## 9. Conclusiones

A partir de los resultados obtenidos, se puede concluir que el pivoteo parcial escalado permite mejorar considerablemente la estabilidad numérica del método de eliminación gaussiana, sin generar un aumento significativo en el costo computacional, en cuanto al costo de ejecución, se mantuvo el total de 1486 flops, coincidiendo con la cota de referencia y obteniendo una desviación del 0 %. A este costo se agregaron 275 operaciones correspondientes al pivoteo, lo que representa aproximadamente un 18 % adicional. Sin embargo, este incremento permitió mejorar la precisión de la solución, pasando de 7 a 10 dígitos correctos respecto de la ejecución sin pivoteo, por otro lado, el refinamiento iterativo permitió alcanzar 15 dígitos correctos con solo 288 flops adicionales por pasada lo que nos demuestra que es posible mejorar considerablemente la precisión aprovechando la factorización LU previamente calculada, sin necesidad de repetir todo el proceso de eliminación. 

Una de las principales conclusiones del análisis es que el costo adicional de incorporar estabilidad resulta pequeño en comparación con los errores que pueden producirse al trabajar sin pivoteo. Asimismo, se observó que la estimación basada en  κ·ε no representa adecuadamente la precisión obtenida en este caso, ya que no considera las características particulares de la matriz, especialmente las diferencias de escala entre sus filas. Finalmente, al considerar un sistema de mayor dimensión, como uno de tamaño (n=1000), el costo de la eliminación gaussiana sería aproximadamente 6.7e8 flops, mientras que la búsqueda del pivote tendría un costo de orden (O(n²)), por lo tanto, su impacto sería relativamente pequeño frente al costo total del procedimiento, manteniendo conveniente su utilización. Para este tipo de problemas sería recomendable emplear una factorización LU por bloques mediante bibliotecas optimizadas, realizar un equilibrado previo de filas y utilizar estimadores eficientes del número de condición, como los métodos de Hager y Higham, evitando calcular todos los valores singulares de la matriz.

## 10. Reproducibilidad

Para facilitar la revisión y ejecución del trabajo, se registran las herramientas, versiones y procedimientos utilizados durante el desarrollo del laboratorio.
**Entorno de ejecución:**
- Python 3.12.3.
- NumPy 2.5.4.
- mpmath 1.3.0.
Las versiones utilizadas durante la ejecución se registran por separado de las dependencias fijadas en el archivo `requirements.txt`, donde se especifican NumPy 2.2.6 y mpmath 1.3.0.

**Comando de ejecución:**
Desde el directorio raíz del repositorio se debe ejecutar:
`python src/main.py`
Este comando ejecuta el programa principal y genera los archivos con los resultados del procedimiento:
`resultados/solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.
De esta manera, es posible reproducir el procedimiento, revisar los resultados y comprobar el comportamiento de las dos estrategias de pivoteo implementadas.

## Fuentes

- Enunciado de la Prueba Parcial 3, «El pivote», asignatura Computación Numérica (INF-321), Universidad Católica del Maule, noviembre de 2026.
- Material de clases de la asignatura Computación Numérica, correspondiente a los contenidos de eliminación gaussiana, factorización, estrategias de pivoteo y condicionamiento de sistemas lineales.

