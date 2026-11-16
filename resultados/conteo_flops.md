# Conteo de operaciones (flops)

Convención: 1 flop = 1 operación aritmética elemental (multiplicación, división, suma o resta). Comparaciones, valores absolutos e intercambios de filas se reportan aparte como sobrecosto de estabilidad. El barrido del paso k va de la columna k en adelante; el multiplicador m(i,k) se calcula una sola vez por fila y se reutiliza. Los totales provienen del contador instrumentado en `src/gauss.py`.

## Eliminación por paso (corrida con pivoteo parcial escalado)

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
| **Eliminación** | **704** | **638** | **1342** |
| Sustitución regresiva | 78 | 66 | 144 |
| **Total base** | **782** | **704** | **1486** |

## Comparación con la cota

| Corrida | Flops contados | Cota de referencia | Desviación | ¿Dentro de ±10 %? |
| --- | ---: | ---: | ---: | --- |
| Sin pivoteo | 1486 | 1486 | +0.00 % | sí |
| Pivoteo parcial escalado | 1486 | 1486 | +0.00 % | sí |

## Sobrecosto de estabilidad (pivoteo parcial escalado)

Se reporta aparte; no está incluido en el total base.

| Concepto | Cantidad |
| --- | ---: |
| Comparaciones para calcular las escalas s_i | 132 |
| Divisiones |a_ik| / s_i | 77 |
| Comparaciones de la búsqueda de pivote | 66 |
| Intercambios de filas realizados | 8 |
| **Sobrecosto total (operaciones de pivoteo)** | **275** |

## Costo del refinamiento iterativo (por pasada)

Reutiliza L y U (no se vuelve a factorizar): una sustitución progresiva y una regresiva por pasada.

- Sustitución progresiva: 132 flops
- Sustitución regresiva: 144 flops
- Actualización x + d: 12 flops
- **Total por pasada: 288 flops** (el residuo en precisión extendida se calcula aparte, con aritmética exacta).

## Indicadores de estabilidad (T2)

| Indicador | Sin pivoteo | Pivoteo parcial escalado |
| --- | ---: | ---: |
| Multiplicador máximo | 2.7526e+13 | 2.4656e+09 |
| Factor de crecimiento | 2.6412e+03 | 1.7832e+00 |
| log2(multiplicador máximo) | 44.6 | 31.2 |
| Bits de mantisa desplazados (máx., de 53) | 18.1 | 6.3 |
| Intercambios de filas | 0 | 8 |
| Error relativo final (norma inf, vs. exacto) | 2.080e-10 | 2.357e-13 |
| Dígitos correctos (mínimo entre las 12 incógnitas) | 7 | 10 |
