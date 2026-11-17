# cn-parcial-gauss12-Navarrete-Gonzalez

Resuelve un sistema lineal denso de 12 incógnitas mediante eliminación gaussiana (con y sin pivoteo parcial escalado), con conteo instrumentado de flops y verificación en precisión extendida.

## Autoría y asignatura

- Elias Gonzalez Roman
- Daniel Navarrete Silva
- Asignatura: Computación Numérica (INF-321), Universidad Católica del Maule. Docentes: Juan Manuel Cabezas Contardo y Pedro Pesante.
- Sección: 2

## Requisitos

- Python 3.10 o superior (probado con 3.13).
- Instalar dependencias:

```
pip install -r requirements.txt
```

## Ejecución

Desde la raíz del repositorio:

```
python src/main.py
```

Al ejecutar esto se genera en `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.

## Estructura de carpetas

```
cn-parcial-gauss12-Navarrete-Gonzalez/
├─ datos/        matriz A (12x12) y vector b, sin truncar
├─ src/          código: eliminación, conteo de flops, verificación y punto de entrada
├─ resultados/   archivos que genera main.py
├─ informe.md    informe técnico y análisis de error
├─ readme.md     este archivo
└─ requirements.txt
└─ .gitignore
```

## Resultado resumido

| Incógnita | Valor (10 dígitos) |
| --- | ---: |
| x1 | 7.384915627e-01 |
| x2 | -1.205937417e+00 |
| x3 | 4.387195200e-03 |
| x4 | -9.172645380e-02 |
| x5 | 2.846109372e+00 |
| x6 | -6.281540900e-03 |
| x7 | 1.473820592e+00 |
| x8 | -3.849271640e-02 |
| x9 | 9.152648037e-01 |
| x10 | -2.607139485e+00 |
| x11 | 5.187362940e-02 |
| x12 | -1.084937162e+00 |

- Total de flops: **1486**.
- Norma infinito del residuo : **2.6e-10**.

## Informe

El análisis completo se encuentra ubicado en [informe.md](informe.md).

## Licencia

Trabajo de uso académico para la asignatura Computación Numérica (UCM). No se autoriza su reutilización fuera de ese contexto.