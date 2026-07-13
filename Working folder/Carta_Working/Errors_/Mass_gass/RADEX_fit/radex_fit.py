# -*- coding: utf-8 -*-
"""
Ajuste de flujos integrados de CO con RADEX (via pythonradex, validado
contra RADEX online en test_validacion.py).

Idea: para cada fuente tenemos 2 líneas (12CO y 13CO de la misma transición,
más 12CO 1-0 extra en las fuentes del sur). Se varían como máximo 2
parámetros del modelo hasta minimizar la diferencia entre el flujo integrado
que predice RADEX y el observado (chi2).

OJO: RADEX no usa "abundancia" como input. Sus parámetros son:
    N (densidad columnar), Tkin, n(H2), dv (ancho de línea).
La abundancia f_CO = CO/H2 entra DESPUÉS, para convertir N(CO) -> N(H2) -> masa.
Por eso el loop varía N, no la abundancia.

Modos (elegir qué 2 parámetros se ajustan):
  A (default): Tkin y n(H2) fijos -> ajusta N(12CO) y N(13CO), una por línea
               (biseccion; cada línea determina su columna). Salida clave:
               ratio isotópico N12/N13 (comparar con X de Ramstedt), tau y
               Tex calculados por RADEX (reemplaza el supuesto Tex=20 K).
  B: n(H2) fijo y N13 = N12/X_iso -> ajusta N(12CO) y Tkin (grilla chi2 conjunta).
  C: Tkin fijo y N13 = N12/X_iso -> ajusta N(12CO) y n(H2) (grilla chi2 conjunta).

Uso:
    python radex_fit.py                     # modo A, todas las fuentes
    python radex_fit.py --modo B
    python radex_fit.py --modo B --fuente RZ_Sgr --fuente R_Hor
    python radex_fit.py --modo B --tres-lineas   # incluye 12CO(1-0) en el chi2
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from scipy import constants
from scipy.optimize import brentq

from pythonradex import radiative_transfer, helpers

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RESULTADOS = os.path.join(HERE, "resultados")

ARCSEC = np.pi / (180.0 * 3600.0)          # rad
M_H2 = 2.016 * 1.66054e-27                 # kg
MSUN = 1.98892e30                          # kg
PC = 3.0857e16                             # m
FACTOR_HE = 1.36                           # M_gas = 1.36 M_H2 (igual que en Mass_gass)

DATAFILES = {"12CO": os.path.join(DATA, "co.dat"),
             "13CO": os.path.join(DATA, "13co.dat")}
TRANS_IDX = {"1-0": 0, "2-1": 1}           # orden de las transiciones en el LAMDA

# rangos de búsqueda
LOGN_MIN, LOGN_MAX = 12.0, 20.5            # log10 N(12CO o 13CO) [cm-2]
TKIN_MIN, TKIN_MAX = 5.0, 300.0            # K
LOGNH2_MIN, LOGNH2_MAX = 3.0, 9.0          # log10 n(H2) [cm-3]


def omega_gauss(theta_fwhm_arcsec):
    """Ángulo sólido [sr] de una fuente gaussiana circular de FWHM dado."""
    th = theta_fwhm_arcsec * ARCSEC
    return np.pi * th**2 / (4.0 * np.log(2.0))


def split_ortho_para(n_H2_cm3, Tkin):
    """Reparte n(H2) en orto/para con la razón térmica (convención RADEX)."""
    opr = min(3.0, 9.0 * np.exp(-170.6 / Tkin))
    fo = opr / (1.0 + opr)
    n_m3 = n_H2_cm3 * 1e6
    return {"para-H2": (1.0 - fo) * n_m3, "ortho-H2": fo * n_m3}


class Modelo:
    """Un objeto RADEX (pythonradex Source) por molécula y ancho de línea."""

    _cache = {}

    def __new__(cls, molecula, dv_kms):
        key = (molecula, round(float(dv_kms), 3))
        if key not in cls._cache:
            obj = super().__new__(cls)
            obj._init(molecula, dv_kms)
            cls._cache[key] = obj
        return cls._cache[key]

    def _init(self, molecula, dv_kms):
        self.molecula = molecula
        self.src = radiative_transfer.Source(
            datafilepath=DATAFILES[molecula],
            geometry="LVG sphere RADEX",
            line_profile_type="rectangular",
            width_v=dv_kms * 1e3,
            warn_negative_tau=False,
        )
        self.cmb = helpers.generate_CMB_background(z=0)
        self._bg_set = False

    def flujo(self, N_cm2, Tkin, n_H2_cm3, trans, omega_sr):
        """Flujo integrado predicho [Jy km/s] + (tau, Tex) de la transición.

        Al flujo de línea de pythonradex se le resta el fondo CMB absorbido,
        (1-e^-tau)*I_CMB, porque la observación interferométrica con continuo
        restado mide (S - I_CMB)(1 - e^-tau)*Omega.
        """
        kwargs = dict(
            N=N_cm2 * 1e4,
            Tkin=float(Tkin),
            collider_densities=split_ortho_para(n_H2_cm3, Tkin),
            T_dust=0, tau_dust=0,
        )
        if not self._bg_set:
            kwargs["ext_background"] = self.cmb
            self._bg_set = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.src.update_parameters(**kwargs)
            self.src.solve_radiative_transfer()
            i = TRANS_IDX[trans]
            nu0 = float(self.src.emitting_molecule.nu0[i])
            tau = float(self.src.tau_nu0_individual_transitions[i])
            Tex = float(self.src.Tex[i])
            F = float(self.src.frequency_integrated_emission(
                "flux", transitions=i, solid_angle=omega_sr))          # W/m2
        dnu = nu0 * self.src.emitting_molecule.width_v / constants.c
        F -= float(self.cmb(nu0)) * (1.0 - np.exp(-tau)) * dnu * omega_sr
        F_jykms = F * 1e26 * (constants.c / 1e3) / nu0
        return F_jykms, tau, Tex


def cargar_observaciones(path, incluir_1_0=False):
    obs = pd.read_csv(path)
    obs["molecula"] = obs["linea"].str.split("_").str[0]
    obs["trans"] = obs["linea"].str.split("_").str[1]
    fuentes = {}
    for fuente, g in obs.groupby("fuente", sort=False):
        # la transición "principal" es la que tiene 12CO y 13CO a la vez
        trans_principal = None
        for t in ("2-1", "1-0"):
            mols = set(g.loc[g["trans"] == t, "molecula"])
            if {"12CO", "13CO"} <= mols:
                trans_principal = t
                break
        if trans_principal is None:
            print(f"[AVISO] {fuente}: no tiene 12CO y 13CO de una misma "
                  "transición; se omite.")
            continue
        sel = g[g["trans"] == trans_principal]
        if incluir_1_0 and trans_principal == "2-1":
            extra = g[(g["trans"] == "1-0") & (g["molecula"] == "12CO")]
            sel = pd.concat([sel, extra])
        fuentes[fuente] = sel.reset_index(drop=True)
    return fuentes


def chi2_fuente(filas, N12_cm2, Tkin, n_H2_cm3, X_iso):
    """chi2 de todas las líneas de una fuente para un punto (N12, Tkin, n)."""
    chi2 = 0.0
    detalles = []
    for _, fila in filas.iterrows():
        N = N12_cm2 if fila["molecula"] == "12CO" else N12_cm2 / X_iso
        m = Modelo(fila["molecula"], fila["dv_kms"])
        try:
            F, tau, Tex = m.flujo(N, Tkin, n_H2_cm3, fila["trans"],
                                  omega_gauss(fila["theta_fwhm_arcsec"]))
        except Exception:
            return np.inf, []
        sigma = fila["err_frac"] * fila["flujo_Jy_kms"]
        chi2 += ((F - fila["flujo_Jy_kms"]) / sigma) ** 2
        detalles.append((fila["linea"], F, fila["flujo_Jy_kms"], tau, Tex))
    return chi2, detalles


def masa_gas(N12_cm2, f_CO, theta_fwhm_arcsec, D_pc):
    """M_gas [Msun] a partir de la N(12CO) promediada en Omega."""
    N_H2 = N12_cm2 / f_CO * 1e4                     # m-2
    omega = omega_gauss(theta_fwhm_arcsec)
    return N_H2 * omega * (D_pc * PC) ** 2 * M_H2 * FACTOR_HE / MSUN


# ----------------------------------------------------------------- modo A
def ajustar_modo_A(fuente, filas):
    """Tkin y n fijos; una N por línea (bisección sobre log10 N)."""
    Tkin = float(filas["Tkin_fijo_K"].iloc[0])
    n_H2 = float(filas["n_H2_cm3"].iloc[0])
    res = {"fuente": fuente, "modo": "A", "Tkin_K": Tkin, "n_H2_cm3": n_H2}
    print(f"\n=== {fuente} (modo A: Tkin={Tkin} K, n(H2)={n_H2:.1e} cm-3) ===")
    for _, fila in filas.iterrows():
        m = Modelo(fila["molecula"], fila["dv_kms"])
        omega = omega_gauss(fila["theta_fwhm_arcsec"])
        F_obs = fila["flujo_Jy_kms"]

        def dif(logN):
            return m.flujo(10**logN, Tkin, n_H2, fila["trans"], omega)[0] - F_obs

        try:
            lo, hi = dif(LOGN_MIN), dif(LOGN_MAX)
            if lo > 0 or hi < 0:
                raise ValueError(
                    f"flujo observado fuera del rango del modelo "
                    f"[{lo + F_obs:.3g}, {hi + F_obs:.3g}] Jy km/s — "
                    "revisar theta/dv/Tkin")
            logN = brentq(dif, LOGN_MIN, LOGN_MAX, xtol=1e-4)
        except ValueError as e:
            print(f"  {fila['linea']:9s}: SIN SOLUCION ({e})")
            res[f"N_{fila['linea']}_cm2"] = np.nan
            continue
        N = 10**logN
        F, tau, Tex = m.flujo(N, Tkin, n_H2, fila["trans"], omega)
        print(f"  {fila['linea']:9s}: N = {N:.3e} cm-2 | tau = {tau:.3f} | "
              f"Tex = {Tex:.1f} K | F_modelo = {F:.2f} vs F_obs = {F_obs:.2f}")
        res[f"N_{fila['linea']}_cm2"] = N
        res[f"tau_{fila['linea']}"] = tau
        res[f"Tex_{fila['linea']}_K"] = Tex

    lineas = list(filas["linea"])
    l12 = [l for l in lineas if l.startswith("12CO")][0]
    l13 = [l for l in lineas if l.startswith("13CO")][0]
    N12, N13 = res.get(f"N_{l12}_cm2"), res.get(f"N_{l13}_cm2")
    if N12 and N13 and np.isfinite(N12) and np.isfinite(N13):
        res["ratio_N12_N13"] = N12 / N13
        res["X_iso_Ramstedt"] = filas["X_iso"].iloc[0]
        print(f"  -> N(12CO)/N(13CO) ajustado = {N12 / N13:.1f}  "
              f"(X Ramstedt = {filas['X_iso'].iloc[0]})")
    if N12 and np.isfinite(N12):
        fila12 = filas[filas["linea"] == l12].iloc[0]
        res["M_gas_Msun"] = masa_gas(N12, fila12["f_CO"],
                                     fila12["theta_fwhm_arcsec"], fila12["D_pc"])
        print(f"  -> M_gas = {res['M_gas_Msun']:.3e} Msun "
              f"(f_CO={fila12['f_CO']}, PROVISIONAL si theta es placeholder)")
    return res


# ------------------------------------------------------------- modos B y C
def ajustar_grilla(fuente, filas, modo, n_grid=(40, 30), sufijo=""):
    """Grilla chi2 en (N12, Tkin) [modo B] o (N12, n_H2) [modo C]."""
    X_iso = float(filas["X_iso"].iloc[0])
    Tkin_fijo = float(filas["Tkin_fijo_K"].iloc[0])
    n_fijo = float(filas["n_H2_cm3"].iloc[0])

    logN = np.linspace(LOGN_MIN + 1, LOGN_MAX - 0.5, n_grid[0])
    if modo == "B":
        eje2 = np.geomspace(TKIN_MIN, TKIN_MAX, n_grid[1])
        nombre2, unidad2 = "Tkin", "K"
    else:
        eje2 = np.geomspace(10**LOGNH2_MIN, 10**LOGNH2_MAX, n_grid[1])
        nombre2, unidad2 = "n(H2)", "cm-3"

    print(f"\n=== {fuente} (modo {modo}: ajusta N(12CO) y {nombre2}; "
          f"X_iso={X_iso}) ===")
    chi2 = np.full((len(logN), len(eje2)), np.inf)
    for i, lN in enumerate(logN):
        for j, p2 in enumerate(eje2):
            Tkin = p2 if modo == "B" else Tkin_fijo
            n_H2 = n_fijo if modo == "B" else p2
            chi2[i, j], _ = chi2_fuente(filas, 10**lN, Tkin, n_H2, X_iso)

    i0, j0 = np.unravel_index(np.nanargmin(chi2), chi2.shape)
    # refinamiento local con una grilla fina alrededor del mínimo
    dlN = logN[1] - logN[0]
    dl2 = np.log10(eje2[1] / eje2[0])
    logN_f = np.linspace(logN[i0] - dlN, logN[i0] + dlN, 15)
    eje2_f = np.geomspace(eje2[j0] * 10**-dl2, eje2[j0] * 10**dl2, 15)
    chi2_f = np.full((15, 15), np.inf)
    for i, lN in enumerate(logN_f):
        for j, p2 in enumerate(eje2_f):
            Tkin = p2 if modo == "B" else Tkin_fijo
            n_H2 = n_fijo if modo == "B" else p2
            chi2_f[i, j], _ = chi2_fuente(filas, 10**lN, Tkin, n_H2, X_iso)
    i1, j1 = np.unravel_index(np.nanargmin(chi2_f), chi2_f.shape)

    N12 = 10**logN_f[i1]
    p2 = eje2_f[j1]
    Tkin = p2 if modo == "B" else Tkin_fijo
    n_H2 = n_fijo if modo == "B" else p2
    c2min, detalles = chi2_fuente(filas, N12, Tkin, n_H2, X_iso)

    print(f"  mejor ajuste: N(12CO) = {N12:.3e} cm-2, "
          f"{nombre2} = {p2:.4g} {unidad2}, chi2 = {c2min:.2f} "
          f"({len(filas)} lineas, {len(filas) - 2} g.d.l.)")
    for linea, F, F_obs, tau, Tex in detalles:
        print(f"    {linea:9s}: F_modelo = {F:8.2f} vs F_obs = {F_obs:8.2f} "
              f"Jy km/s | tau = {tau:.3f} | Tex = {Tex:.1f} K")

    fila12 = filas[filas["molecula"] == "12CO"].iloc[0]
    Mg = masa_gas(N12, fila12["f_CO"], fila12["theta_fwhm_arcsec"],
                  fila12["D_pc"])
    print(f"  -> M_gas = {Mg:.3e} Msun (PROVISIONAL si theta es placeholder)")

    guardar_mapa_chi2(fuente, modo + sufijo, logN, eje2, chi2, np.log10(N12),
                      p2, nombre2, unidad2)

    res = {"fuente": fuente, "modo": modo, "N_12CO_cm2": N12,
           "chi2_min": c2min, "M_gas_Msun": Mg, "X_iso": X_iso}
    res["Tkin_K" if modo == "B" else "n_H2_cm3"] = p2
    res["n_H2_cm3" if modo == "B" else "Tkin_K"] = n_fijo if modo == "B" else Tkin_fijo
    for linea, F, F_obs, tau, Tex in detalles:
        res[f"F_mod_{linea}"] = F
        res[f"tau_{linea}"] = tau
        res[f"Tex_{linea}_K"] = Tex
    return res


def guardar_mapa_chi2(fuente, modo, logN, eje2, chi2, logN_best, p2_best,
                      nombre2, unidad2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c2 = np.where(np.isfinite(chi2), chi2, np.nanmax(chi2[np.isfinite(chi2)]))
    c2min = np.nanmin(c2)
    fig, ax = plt.subplots(figsize=(7, 5))
    pc = ax.pcolormesh(logN, np.log10(eje2), np.log10(c2).T, shading="auto",
                       cmap="viridis")
    fig.colorbar(pc, ax=ax, label="log10 chi2")
    # contornos 1, 2 y 3 sigma para 2 parámetros
    ax.contour(logN, np.log10(eje2), (c2 - c2min).T,
               levels=[2.30, 6.17, 11.8], colors="w", linewidths=1)
    ax.plot(logN_best, np.log10(p2_best), "r*", ms=14,
            label="mejor ajuste")
    ax.set_xlabel("log10 N(12CO) [cm-2]")
    ax.set_ylabel(f"log10 {nombre2} [{unidad2}]")
    ax.set_title(f"{fuente} — modo {modo} (contornos: 1, 2, 3 sigma)")
    ax.legend()
    out = os.path.join(RESULTADOS, f"chi2_{fuente}_modo{modo}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  mapa chi2 guardado en resultados/{os.path.basename(out)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--modo", choices=["A", "B", "C"], default="A")
    ap.add_argument("--fuente", action="append",
                    help="ajustar solo esta(s) fuente(s)")
    ap.add_argument("--tres-lineas", action="store_true",
                    help="incluir el 12CO(1-0) extra en el chi2 (fuentes 2-1)")
    ap.add_argument("--obs", default=os.path.join(HERE, "observaciones.csv"))
    args = ap.parse_args()

    os.makedirs(RESULTADOS, exist_ok=True)
    fuentes = cargar_observaciones(args.obs, incluir_1_0=args.tres_lineas)
    if args.fuente:
        fuentes = {k: v for k, v in fuentes.items() if k in args.fuente}
        if not fuentes:
            raise SystemExit("ninguna fuente coincide con --fuente")

    hay_placeholder = any((v["theta_estado"] == "BEAM_PLACEHOLDER").any()
                          for v in fuentes.values())
    if hay_placeholder:
        print("*" * 72)
        print("AVISO: theta_fwhm_arcsec es un PLACEHOLDER (haz de ACA) para al")
        print("menos una fuente. Las N y masas absolutas NO son confiables hasta")
        print("medir el tamano real de la region emisora (ver README, paso 1).")
        print("*" * 72)

    resultados = []
    for fuente, filas in fuentes.items():
        if args.modo == "A":
            resultados.append(ajustar_modo_A(fuente, filas))
        else:
            resultados.append(ajustar_grilla(
                fuente, filas, args.modo,
                sufijo="_3lineas" if args.tres_lineas else ""))

    df = pd.DataFrame(resultados)
    sufijo = "_3lineas" if args.tres_lineas else ""
    if args.fuente:
        sufijo += "_parcial"
    out = os.path.join(RESULTADOS, f"ajuste_modo{args.modo}{sufijo}.csv")
    df.to_csv(out, index=False)
    print(f"\nResultados guardados en resultados/{os.path.basename(out)}")


if __name__ == "__main__":
    main()
