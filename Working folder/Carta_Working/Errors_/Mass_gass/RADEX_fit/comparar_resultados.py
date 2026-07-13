# -*- coding: utf-8 -*-
"""
Junta los resultados de los modos A/B/C de radex_fit.py con el método
actual de Mass_gass (../gas_mass_results.csv) y genera:
    resultados/resumen_comparativo.csv  (tabla completa)
    resultados/resumen_comparativo.md   (tabla legible para la reunión)

Correr después de cada actualización de observaciones.csv:
    python radex_fit.py            # modo A
    python radex_fit.py --modo B
    python radex_fit.py --modo B --tres-lineas
    python radex_fit.py --modo C
    python comparar_resultados.py
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "resultados")

# método actual (nombres de fuente con espacio -> con guion bajo)
actual = pd.read_csv(os.path.join(HERE, "..", "gas_mass_results.csv"))
actual["fuente"] = actual["Source"].str.replace(" ", "_")

def cargar(nombre):
    path = os.path.join(RES, nombre)
    return pd.read_csv(path) if os.path.exists(path) else None

mA = cargar("ajuste_modoA.csv")
mB = cargar("ajuste_modoB.csv")
mB3 = cargar("ajuste_modoB_3lineas.csv")
mC = cargar("ajuste_modoC.csv")

def col_linea(df, fuente, prefijo, molecula):
    """Valor de la columna '<prefijo>_<molecula>_<trans>...' que exista."""
    fila = df[df["fuente"] == fuente]
    if fila.empty:
        return np.nan
    fila = fila.iloc[0]
    for c in df.columns:
        if c.startswith(f"{prefijo}_{molecula}"):
            v = fila[c]
            if pd.notna(v):
                return v
    return np.nan

filas = []
for fuente in mA["fuente"]:
    a = mA[mA["fuente"] == fuente].iloc[0]
    act = actual[actual["fuente"] == fuente]
    fila = {
        "fuente": fuente,
        # método actual (analítico, tau del isotopólogo, Tex=20 K asumida)
        "M_gas_actual": act["M_gas_Msun"].iloc[0] if not act.empty else np.nan,
        "tau_actual": act["tau"].iloc[0] if not act.empty else np.nan,
        # modo A
        "A_N12_cm2": col_linea(mA, fuente, "N", "12CO"),
        "A_tau12": col_linea(mA, fuente, "tau", "12CO"),
        "A_ratio_N12_N13": a.get("ratio_N12_N13", np.nan),
        "X_Ramstedt": a.get("X_iso_Ramstedt", np.nan),
        "A_M_gas": a.get("M_gas_Msun", np.nan),
    }
    if mB is not None:
        b = mB[mB["fuente"] == fuente]
        if not b.empty:
            b = b.iloc[0]
            fila.update({"B_N12_cm2": b["N_12CO_cm2"], "B_Tkin_K": b["Tkin_K"],
                         "B_tau12": col_linea(mB, fuente, "tau", "12CO"),
                         "B_chi2": b["chi2_min"], "B_M_gas": b["M_gas_Msun"]})
    if mB3 is not None:
        b3 = mB3[(mB3["fuente"] == fuente)]
        # solo cuenta como "3 lineas" si tiene el 2-1 Y el 1-0 a la vez
        def _tiene(pref):
            return any(c.startswith(pref) and pd.notna(b3.iloc[0][c])
                       for c in mB3.columns)
        if not b3.empty and _tiene("F_mod_12CO_2-1") and _tiene("F_mod_12CO_1-0"):
            b3 = b3.iloc[0]
            fila.update({"B3_N12_cm2": b3["N_12CO_cm2"],
                         "B3_Tkin_K": b3["Tkin_K"],
                         "B3_chi2": b3["chi2_min"],
                         "B3_M_gas": b3["M_gas_Msun"]})
    if mC is not None:
        c_ = mC[mC["fuente"] == fuente]
        if not c_.empty:
            c_ = c_.iloc[0]
            fila.update({"C_n_H2_cm3": c_["n_H2_cm3"], "C_chi2": c_["chi2_min"]})
    filas.append(fila)

df = pd.DataFrame(filas)
df.to_csv(os.path.join(RES, "resumen_comparativo.csv"), index=False)

def f(x, fmt="{:.2e}"):
    return fmt.format(x) if pd.notna(x) else "—"

lineas_md = [
    "# Resumen comparativo — ajuste RADEX vs método actual",
    "",
    "Ω = placeholder del haz de ACA: masas y N absolutas PROVISIONALES;",
    "los cocientes (ratio, τ) y los chi² son mucho menos sensibles a Ω.",
    "",
    "## Masa de gas (M☉)",
    "",
    "| Fuente | Actual (τ-corr, Tex=20) | Modo A | Modo B | Modo B 3 líneas |",
    "|---|---|---|---|---|",
]
for _, r in df.iterrows():
    lineas_md.append(
        f"| {r['fuente']} | {f(r['M_gas_actual'])} | {f(r['A_M_gas'])} "
        f"| {f(r.get('B_M_gas'))} | {f(r.get('B3_M_gas'))} |")

lineas_md += [
    "",
    "## Opacidad y ratio isotópico",
    "",
    "| Fuente | τ actual (analítico) | τ modo A | τ modo B | N12/N13 (A) "
    "| X Ramstedt | χ² B (2 lín.) |",
    "|---|---|---|---|---|---|---|",
]
for _, r in df.iterrows():
    lineas_md.append(
        f"| {r['fuente']} | {f(r['tau_actual'], '{:.2f}')} "
        f"| {f(r['A_tau12'], '{:.2f}')} | {f(r.get('B_tau12'), '{:.2f}')} "
        f"| {f(r['A_ratio_N12_N13'], '{:.1f}')} "
        f"| {f(r['X_Ramstedt'], '{:.1f}')} | {f(r.get('B_chi2'), '{:.1f}')} |")

lineas_md += [
    "",
    "## Lectura",
    "",
    "- **χ²(B) ≈ 0** (RZ Sgr, IRAS 08500, V510, IRAS 07145): el esquema",
    "  '12CO saturado + X de Ramstedt' es consistente, y el τ que ajusta",
    "  RADEX confirma de forma independiente el τ analítico del isotopólogo",
    "  (RZ Sgr: 7.4 vs 6.25; IRAS 08500: 3.7 vs 3.18).",
    "- **χ²(B) ≫ 1 con ratio observado > X** (R Hor 39>13, RV Aqr 27>13,",
    "  UU For 19>13, IRAS 06531 85>34): matemáticamente imposible de ajustar",
    "  con X fijo — el dato sugiere que el X real de esas fuentes es mayor",
    "  que la mediana de Ramstedt (dispersión de factor ~15 en estrellas C).",
    "  En el método actual estas fuentes daban τ≈0 por la misma razón.",
    "- **IRAS 07180** (ratio 4.9 < X=13 pero χ²=16): para saturar el 12CO",
    "  con un flujo tan débil hace falta una fuente más compacta que el haz",
    "  → se resolverá al medir el tamaño real.",
    "- Modo C: n(H₂) queda esencialmente degenerado (líneas casi",
    "  termalizadas); no aporta constraint útil — esperado.",
]

with open(os.path.join(RES, "resumen_comparativo.md"), "w",
          encoding="utf-8") as fh:
    fh.write("\n".join(lineas_md) + "\n")

print(df.to_string(index=False,
                   float_format=lambda x: f"{x:.3g}"))
print("\nGuardado: resultados/resumen_comparativo.csv y .md")
