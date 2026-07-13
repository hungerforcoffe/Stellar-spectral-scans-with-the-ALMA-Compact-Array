# Resumen comparativo — ajuste RADEX vs método actual

Ω = placeholder del haz de ACA: masas y N absolutas PROVISIONALES;
los cocientes (ratio, τ) y los chi² son mucho menos sensibles a Ω.

## Masa de gas (M☉)

| Fuente | Actual (τ-corr, Tex=20) | Modo A | Modo B | Modo B 3 líneas |
|---|---|---|---|---|
| IRAS_06531 | 7.89e-04 | 8.34e-04 | 2.21e-03 | 5.69e-04 |
| IRAS_07145 | 5.64e-03 | 4.54e-03 | 2.61e-02 | 7.12e-03 |
| IRAS_07180 | 3.01e-04 | 1.05e-04 | 9.08e-04 | 1.96e-04 |
| IRAS_08500 | 4.10e-02 | 1.50e-02 | 6.57e-02 | 5.19e-02 |
| V510_Pup | 8.77e-04 | 5.65e-04 | 8.92e-04 | 8.41e-04 |
| RV_Aqr | 2.55e-02 | 3.02e-02 | 7.65e-02 | — |
| RZ_Sgr | 8.46e-03 | 1.55e-03 | 7.71e-03 | — |
| UU_For | 7.71e-03 | 8.04e-03 | 4.34e-02 | — |
| R_Hor | 1.20e-03 | 1.36e-03 | 2.57e-03 | — |

## Opacidad y ratio isotópico

| Fuente | τ actual (analítico) | τ modo A | τ modo B | N12/N13 (A) | X Ramstedt | χ² B (2 lín.) |
|---|---|---|---|---|---|---|
| IRAS_06531 | 0.00 | 0.10 | 0.00 | 77.2 | 34.0 | 17.3 |
| IRAS_07145 | 0.52 | 0.04 | 0.65 | 23.8 | 34.0 | 0.0 |
| IRAS_07180 | 2.72 | 0.02 | 0.62 | 4.3 | 13.0 | 10.8 |
| IRAS_08500 | 3.18 | 0.37 | 3.74 | 4.4 | 13.0 | 0.0 |
| V510_Pup | 1.00 | 0.03 | 0.10 | 10.6 | 12.1 | 0.1 |
| RV_Aqr | 0.00 | 0.33 | 0.02 | 27.3 | 13.0 | 12.0 |
| RZ_Sgr | 6.25 | 0.29 | 7.38 | 4.5 | 26.0 | 0.0 |
| UU_For | 0.00 | 0.10 | 0.00 | 16.8 | 13.0 | 1.4 |
| R_Hor | 0.00 | 0.25 | 0.01 | 37.8 | 13.0 | 27.8 |

## Lectura

- **χ²(B) ≈ 0** (RZ Sgr, IRAS 08500, V510, IRAS 07145): el esquema
  '12CO saturado + X de Ramstedt' es consistente, y el τ que ajusta
  RADEX confirma de forma independiente el τ analítico del isotopólogo
  (RZ Sgr: 7.4 vs 6.25; IRAS 08500: 3.7 vs 3.18).
- **χ²(B) ≫ 1 con ratio observado > X** (R Hor 39>13, RV Aqr 27>13,
  UU For 19>13, IRAS 06531 85>34): matemáticamente imposible de ajustar
  con X fijo — el dato sugiere que el X real de esas fuentes es mayor
  que la mediana de Ramstedt (dispersión de factor ~15 en estrellas C).
  En el método actual estas fuentes daban τ≈0 por la misma razón.
- **IRAS 07180** (ratio 4.9 < X=13 pero χ²=16): para saturar el 12CO
  con un flujo tan débil hace falta una fuente más compacta que el haz
  → se resolverá al medir el tamaño real.
- Modo C: n(H₂) queda esencialmente degenerado (líneas casi
  termalizadas); no aporta constraint útil — esperado.
