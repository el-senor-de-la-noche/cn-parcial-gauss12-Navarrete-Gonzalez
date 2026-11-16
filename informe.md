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


## 5. Análisis de mantisa y pivoteo

## 6. Condicionamiento

## 7. Resultados

## 8. Verificación

## 9. Conclusiones

## 10. Reproducibilidad