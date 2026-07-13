# Ajuste RADEX de los flujos integrados de CO

Carpeta de trabajo para lo que pidió Luciano: *"instalar RADEX, tomar los
flujos integrados de cada línea y fuente, y ponerlos en un loop donde se
cambia la abundancia (o densidad columnar y algún otro parámetro, máximo
dos) hasta reproducir al mismo tiempo las dos líneas por fuente"*.

**Estado: ya funciona.** RADEX está instalado (vía `pythonradex`, ver §5),
validado contra el run de RADEX online de Luciano (§6), y el loop de ajuste
está implementado en `radex_fit.py`, corrido para las 9 fuentes con los
**beams reales de los cubos** (tomados de los notebooks) como ángulo sólido
(§4). Resultados en `resultados/` y tabla consolidada en
`resultados/resumen_comparativo.md`.

---

## 1. Qué hay en esta carpeta

| Archivo | Qué es |
|---|---|
| `radex_fit.py` | El script del loop de ajuste (3 modos, ver §3) |
| `observaciones.csv` | Flujos integrados por fuente y línea + parámetros por fuente (editable) |
| `comparar_resultados.py` | Junta los modos A/B/C con el método actual → `resultados/resumen_comparativo.{csv,md}` |
| `test_validacion.py` | Reproduce el run de RADEX online de Luciano (IRAS 06531) |
| `data/co.dat`, `data/13co.dat` | Datos moleculares LAMDA (Leiden) |
| `resultados/` | Tablas de ajuste (CSV) y mapas de chi² (PNG); corridas con `--fuente` van a `*_parcial.csv` |

Los flujos de `observaciones.csv` vienen de `../../line_ratios.csv`; las
distancias, f_CO y X isotópicos de `../gas_mass_results.csv`; los anchos de
línea de los ajustes de los notebooks (dv = 2·V_e).

## 2. La idea clave: qué se ajusta y por qué

RADEX **no tiene "abundancia" como parámetro de entrada**. Sus inputs son:

- **N** — densidad columnar de la molécula (cm⁻²)  ← *esto es lo que se varía en el loop*
- **T_kin** — temperatura cinética del gas
- **n(H₂)** — densidad volumétrica del colisionador
- **dv** — ancho de línea
- T_bg — fondo (CMB, 2.73 K)

La abundancia f_CO = CO/H₂ entra **después**, al convertir la N(CO) ajustada
en N(H₂) y en masa: M_gas ∝ N(CO)/f_CO. O sea: el "loop en abundancia" del
correo se traduce en la práctica en un **loop en densidad columnar**; la
abundancia sigue siendo el supuesto externo de siempre (Ramstedt 2008).

Con **dos líneas por fuente** (¹²CO y ¹³CO de la misma transición) se pueden
ajustar **como máximo dos parámetros**. Las combinaciones útiles son:

## 3. Los tres modos del script

### Modo A (default, recomendado como primer paso): N(¹²CO) y N(¹³CO)
`python radex_fit.py`

T_kin y n(H₂) fijos (columnas del CSV, por defecto 20 K y 10⁵ cm⁻³). Cada
línea determina su columna por bisección (el flujo es monótono en N).
**Qué gana esto frente al método actual de Mass_gass:**
- la **opacidad τ sale de RADEX** de forma autoconsistente, sin asumir el X
  de Ramstedt (reemplaza el τ analítico de Myers/Goldsmith-Langer);
- la **T_ex se calcula** (no-LTE) en lugar de asumirse igual a 20 K;
- el cociente ajustado **N(¹²CO)/N(¹³CO) es un resultado**, comparable con
  el X de Ramstedt — si sale parecido al ratio de flujos observado, las
  líneas son delgadas; si sale mayor, hay saturación del ¹²CO.

### Modo B: N(¹²CO) y T_kin (ajuste conjunto de las dos líneas)
`python radex_fit.py --modo B`

Se fija N(¹³CO) = N(¹²CO)/X (X de Ramstedt) y n(H₂), y se busca el mínimo
de chi² en una grilla (N, T_kin) — mapa guardado como PNG con contornos
1/2/3σ. Aquí las dos líneas se ajustan *simultáneamente* de verdad: el
¹³CO (delgado) fija N, y el ¹²CO (más saturado) responde a T_kin.
**Advertencia:** si las dos líneas son ópticamente delgadas, el cociente del
modelo queda clavado en ≈X y no hay forma de reproducir un ratio observado
distinto de X — el chi² mínimo queda alto (le pasa a R Hor: ratio observado
39 vs X=13, chi²≈28). Eso no es un fallo del código: es el dato diciendo que
X no es 13 para esa fuente, o que la calibración relativa tiene un problema.

### Modo C: N(¹²CO) y n(H₂)
`python radex_fit.py --modo C`

Igual que B pero variando la densidad a T_kin fijo. Útil como test de
sensibilidad: a las densidades típicas de CSE (≥10⁵ cm⁻³) el CO(1-0) y
(2-1) están casi termalizados (n_crit ~ 10³⁻⁴ cm⁻³), así que n(H₂) casi no
apalanca — esperable que quede degenerado (contornos abiertos en el mapa).

### Extra: tres líneas para las fuentes del sur
`python radex_fit.py --modo B --tres-lineas`

Las 5 fuentes con ¹³CO(2-1) (los 4 IRAS + V510) **también tienen ¹²CO(1-0)**
→ 3 líneas y 2 parámetros = ajuste sobredeterminado, el más informativo de
todos. El cociente ¹²CO(2-1)/(1-0) es un termómetro/densitómetro real. Es
la opción que recomendaría discutir con Luciano como resultado principal
para esas 5 fuentes.

## 4. El punto delicado: ángulo sólido (unidades)

RADEX predice **intensidades** (equivalente a K km/s); las observaciones
están en **Jy km/s integrados en una región**. Para compararlos hace falta
el ángulo sólido Ω de la región emisora:

    F_pred [Jy km/s] ∝ I_RADEX × Ω

- En `observaciones.csv` la columna `theta_fwhm_arcsec` (FWHM gaussiano;
  Ω = 1.133·θ²) usa los **beams reales de los cubos** tomados de los
  headers FITS que están en los notebooks (BMAJ/BMIN por fuente y
  transición, θ = √(BMAJ·BMIN) — la misma convención de los notebooks;
  con la media geométrica, Ω = 1.133·θ² es exactamente el Ω del haz
  elíptico). Estado: `BEAM_NOTEBOOK`. Para el ¹³CO se usa el beam del
  ¹²CO de la misma transición (diferencia ~4%, despreciable).
- **Qué significa usar el beam:** se asume que la emisión llena ~el haz.
  Si la fuente es más compacta, la N ajustada es un promedio en el haz
  (límite inferior de la N real) y el τ verdadero es mayor. Ojo: en el
  régimen ópticamente delgado la **masa no depende de Ω** (M ∝ N·Ω y
  N ∝ F/Ω se cancelan) — solo la corrección por saturación depende del
  tamaño. Refinamiento futuro: medir el tamaño emisor real en CARTA
  (momento 0 → gaussiana 2D) para las fuentes resueltas.

## 5. Instalación (ya hecha)

RADEX clásico es Fortran y no compila fácil en Windows; en su lugar usamos
**pythonradex 2.0.2** (`pip install pythonradex`), una reimplementación
publicada del mismo cálculo (probabilidad de escape / LVG) que lee los
mismos archivos LAMDA y tiene un modo `"LVG sphere RADEX"` que usa
exactamente las ecuaciones del RADEX original. Requiere solo
numpy/scipy/numba (ya instalados con Python 3.14).

## 6. Validación contra el RADEX online de Luciano

`python test_validacion.py` reproduce el run guardado en
`Working folder/parameters_RADEX/IRAS_06531_parameters.txt`
(CO, T_kin=130 K, n=10⁶ cm⁻³, N=10¹⁷ cm⁻², dv=57.22 km/s, LVG):

| Transición | T_ex (K) | τ | T_R (K) |
|---|---|---|---|
| 1–0 RADEX online | 133.699 | 2.263e-3 | 0.2941 |
| 1–0 pythonradex | 133.721 | 2.263e-3 | 0.2941 |
| 2–1 RADEX online | 132.464 | 8.580e-3 | 1.0830 |
| 2–1 pythonradex | 132.383 | 8.585e-3 | 1.0833 |

Dos convenciones que hubo que igualar (ya incorporadas al script):
1. RADEX normaliza τ con perfil gaussiano de FWHM=dv → factor
   1.0645=√(π/4ln2) respecto al perfil rectangular LVG puro. El test lo
   muestra con y sin el factor. En `radex_fit.py` se usa la convención
   rectangular autoconsistente de pythonradex (diferencia ~6% en τ,
   irrelevante frente a las demás incertidumbres).
2. Al flujo del modelo se le resta el CMB absorbido, (1−e^−τ)·I_CMB·Ω,
   porque la observación interferométrica (continuo restado) mide eso.

## 7. Qué revisar/llenar en `observaciones.csv`

1. **`theta_fwhm_arcsec`** — ya usa los beams reales de los notebooks (§4).
   Refinamiento opcional: reemplazar por el tamaño emisor medido (momento 0
   → gaussiana 2D) en las fuentes resueltas y poner `theta_estado=MEDIDO`.
2. **dv de V510 Pup** — verificado en el notebook: V_e=101.54±7.67 (1-0) y
   178.06±8.17 km/s (2-1) son ajustes formalmente buenos, pero anchos tan
   grandes (outflow post-AGB) merecen discusión: ¿ajustar solo la componente
   central o todo el perfil con alas?
3. **dv de IRAS 08500 (2-1)** — verificado: el fit del notebook da
   V_e=12.343±0.149, independiente del 1-0 (11.078±0.946). El problema de
   valores duplicados de las notas de reunión era de una tabla externa,
   no de estos fits. Ya está bien en el CSV.
4. **`err_frac`** — ahora 10% (¹²CO) y 15% (¹³CO) genéricos; poner los
   errores reales (calibración ACA ~5-10% + ruido del ajuste) cuando estén.
5. **X_iso de V510** — no hay valor publicado (post-AGB); está puesto el
   ratio observado 12.13 como provisional (equivale a suponerlo delgado).

## 8. Resultados (Ω = haz real de cada cubo)

Los cuatro modos ya corrieron para las 9 fuentes; la tabla consolidada
(modos A/B/B-3líneas/C vs el método actual de `../gas_mass_results.csv`)
está en **`resultados/resumen_comparativo.md`** (regenerar con
`python comparar_resultados.py`). Lo esencial:

- **Las masas del modo A salen del mismo orden** que el método actual en
  todas las fuentes — verificación cruzada satisfactoria.
- **El modo B valida el τ analítico del isotopólogo de forma
  independiente**: en las fuentes donde el ratio observado < X, RADEX
  ajusta perfecto (χ²≈0) y encuentra τ(¹²CO) muy parecido al analítico —
  RZ Sgr: 7.4 vs 6.25; IRAS 08500: 3.7 vs 3.18. Fuerte argumento para el
  paper a favor del método del punto 1.4 de las notas de reunión.
- **Donde el ratio observado > X, el modo B no puede ajustar** (χ²=12–28:
  R Hor 39>13, RV Aqr 27>13, UU For 19>13, IRAS 06531 85>34): con X fijo
  el ratio del modelo nunca supera X. No es fallo del código: o el X real
  de esas fuentes es mayor que la mediana de Ramstedt (la dispersión en
  estrellas C es de factor ~15), o hay un problema de calibración
  relativa. Son las mismas fuentes que daban τ≈0 en el método analítico.
- Con el Ω del haz, en el modo A **todas las líneas salen delgadas**
  (τ<0.4) y el N12/N13 ajustado ≈ ratio de flujos observado. Si la
  emisión real es más compacta que el haz, τ sube y N es un límite
  inferior (la masa en régimen delgado NO depende de Ω, ver §4).
- **Modo C**: n(H₂) queda degenerado (líneas casi termalizadas a
  n≳10⁴ cm⁻³) — como se esperaba, no es un parámetro útil aquí.
- **Modo B con 3 líneas** (fuentes del sur): el cociente ¹²CO(2-1)/(1-0)
  entre bandas depende de que la emisión llene ambos haces por igual
  (Ω difiere ~4-5× entre bandas); las T_kin que salen hay que leerlas
  con esa cautela. Con tamaños emisores medidos por banda será el ajuste
  más potente de todos.

## 9. Limitaciones a declarar

- Modelo de **una zona homogénea** (una T, una n, una N por fuente): las
  CSE tienen gradientes; lo ajustado son valores "característicos". Es el
  mismo nivel de aproximación que Cerrigone y que el método actual.
- Geometría LVG de esfera en expansión — la adecuada para CSE.
- Reparto orto/para de H₂ térmico a T_kin (convención de RADEX).
- Sin polvo ni solapamiento de líneas.
- La masa final sigue dependiendo linealmente de 1/f_CO (supuesto externo).
