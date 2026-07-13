# -*- coding: utf-8 -*-
"""
Test de validación: reproducir con pythonradex el run de RADEX online
que está en Working folder/parameters_RADEX/IRAS_06531_parameters.txt

Parámetros de ese run:
    Molécula: CO | Tbg = 2.73 K | Tkin = 130 K | n(H2) = 1e6 cm-3
    N(CO) = 1e17 cm-2 | dv = 57.22 km/s | método: LVG

Resultados esperados (RADEX online):
    Transición   Tex (K)    tau         T_R (K)
    1-0          133.699    2.263e-03   0.2941
    2-1          132.464    8.580e-03   1.083

Notas de convención (verificadas):
  * RADEX normaliza tau con un perfil gaussiano de FWHM = dv, lo que
    introduce un factor 1.0645 = sqrt(pi/(4 ln2)) respecto al perfil
    rectangular que usa pythonradex en LVG. Para imitar a RADEX hay que
    usar width_v = dv * 1.0645.
  * T_R de RADEX es la temperatura Rayleigh-Jeans de (S - I_bg)(1 - e^-tau),
    es decir, con el fondo CMB restado.
"""
import os
import numpy as np
from pythonradex import radiative_transfer, helpers

HERE = os.path.dirname(os.path.abspath(__file__))
CO_DAT = os.path.join(HERE, "data", "co.dat")

GAUSS_FWHM_FACTOR = np.sqrt(np.pi / (4 * np.log(2)))  # 1.0645, convención RADEX

# --- parámetros del run de referencia (unidades SI: pythonradex usa m, no cm)
Tkin = 130.0                      # K
n_H2 = 1e6 * 1e6                  # cm-3 -> m-3
N_CO = 1e17 * 1e4                 # cm-2 -> m-2
dv = 57.22 * 1e3                  # km/s -> m/s

# RADEX reparte "H2" en orto/para con la razón térmica: opr = 9*exp(-170.6/Tkin)
opr = min(3.0, 9.0 * np.exp(-170.6 / Tkin))
frac_ortho = opr / (1.0 + opr)
collider_densities = {"para-H2": n_H2 * (1 - frac_ortho),
                      "ortho-H2": n_H2 * frac_ortho}

def correr(width_v):
    src = radiative_transfer.Source(
        datafilepath=CO_DAT,
        geometry="LVG sphere RADEX",      # mismas ecuaciones que el RADEX original
        line_profile_type="rectangular",  # LVG exige perfil rectangular
        width_v=width_v,
    )
    src.update_parameters(
        N=N_CO,
        Tkin=Tkin,
        collider_densities=collider_densities,
        ext_background=helpers.generate_CMB_background(z=0),
        T_dust=0, tau_dust=0,
    )
    src.solve_radiative_transfer()
    return src

esperado = {0: (133.699, 2.263e-3, 0.2941), 1: (132.464, 8.580e-3, 1.083)}
nombres = {0: "1-0", 1: "2-1"}

for etiqueta, factor in [("width_v = dv (rectangular puro)", 1.0),
                         ("width_v = dv*1.0645 (imitando RADEX)", GAUSS_FWHM_FACTOR)]:
    src = correr(dv * factor)
    print(f"\n--- {etiqueta} ---")
    print(f"{'Trans':6s} {'Tex(py)':>9s} {'Tex(RADEX)':>11s} {'tau(py)':>11s} "
          f"{'tau(RADEX)':>11s} {'TR(py)':>9s} {'TR(RADEX)':>10s}")
    ok = True
    for i in (0, 1):
        nu0 = src.emitting_molecule.nu0[i]
        Tex = src.Tex[i]
        tau = src.tau_nu0_individual_transitions[i]
        # T_R al estilo RADEX: RJ de (S - I_bg)(1 - e^-tau)
        S = helpers.B_nu(T=Tex, nu=nu0)
        I_bg = helpers.generate_CMB_background(z=0)(nu0)
        TR = helpers.RJ_brightness_temperature(
            specific_intensity=(S - I_bg) * (1 - np.exp(-tau)), nu=nu0)
        e_Tex, e_tau, e_TR = esperado[i]
        print(f"{nombres[i]:6s} {Tex:9.3f} {e_Tex:11.3f} {tau:11.3e} "
              f"{e_tau:11.3e} {float(TR):9.4f} {e_TR:10.4f}")
        if abs(Tex - e_Tex) / e_Tex > 0.02 or abs(tau - e_tau) / e_tau > 0.02 \
           or abs(TR - e_TR) / e_TR > 0.02:
            ok = False
    print("=> coincide con RADEX online dentro del 2%:", "SI" if ok else "NO")
