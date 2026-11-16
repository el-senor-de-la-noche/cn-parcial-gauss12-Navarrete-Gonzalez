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

## 4. Conteo de operaciones

## 5. Análisis de mantisa y pivoteo

## 6. Condicionamiento

## 7. Resultados

## 8. Verificación

## 9. Conclusiones

## 10. Reproducibilidad