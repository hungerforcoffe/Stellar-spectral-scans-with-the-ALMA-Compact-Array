import numpy as np
from scipy import interpolate, integrate
import matplotlib.pyplot as plt
import os
import sys


ARCHIVO_SWS    = "71602537_sws.tbl.txt"          # Cambia por el nombre real de tu archivo SWS
ARCHIVO_LWS    = "71602537_pws.tbl.txt"          # Cambia por el nombre real de tu archivo LWS
ARCHIVO_FILTRO = "Transmition_curve_4.5.tbl.txt"   # Curva de transmisión IRAC 4.5 µm


DISTANCIA_KPC = 1.15835              # ← CAMBIA ESTE VALOR por la distancia real

T_REF    = 2000.0                # Temperatura [K]
R_REF_CM = 5e13                  # Radio [cm]

FILTRO_MIN_UM = 3.9
FILTRO_MAX_UM = 5.2


h  = 6.62607015e-27   # erg·s
c  = 2.99792458e10    # cm/s
kB = 1.380649e-16     # erg/K
KPC_CM = 3.085677581e21  # cm por kpc
JY_PER_CGS = 1e23     # 1 Jy = 1e-23 erg/s/cm²/Hz  →  factor de conversión



def cargar_filtro_irac45(archivo):

    if not os.path.exists(archivo):
        print(f"[ERROR] No se encontró el archivo del filtro: '{archivo}'")
        print("  → Asegúrate de que el archivo esté en el mismo directorio")
        print("    que este script y que el nombre coincida con ARCHIVO_FILTRO.")
        sys.exit(1)

    data   = np.loadtxt(archivo, comments="#")
    wav_um = data[:, 0]   # ya en µm
    resp   = data[:, 1]   # electrons/photon

    # Normalizar al máximo → transmisión efectiva [0, 1]
    trans = resp / resp.max()

    print(f"[Filtro] Cargado: '{archivo}'  "
          f"| Rango: {wav_um.min():.3f} – {wav_um.max():.3f} µm  "
          f"| {len(wav_um)} puntos")
    return wav_um, trans


def cargar_espectro(archivo):

    data = np.loadtxt(archivo, comments="#")
    wav  = data[:, 0]   # µm
    flux = data[:, 1]   # Jy
    err  = data[:, 2]   # Jy
    return wav, flux, err


def unir_sws_lws(archivo_sws, archivo_lws):

    wav_s, flux_s, err_s = cargar_espectro(archivo_sws)
    wav_l, flux_l, err_l = cargar_espectro(archivo_lws)

    wav  = np.concatenate([wav_s,  wav_l])
    flux = np.concatenate([flux_s, flux_l])
    err  = np.concatenate([err_s,  err_l])

    idx  = np.argsort(wav)
    wav, flux, err = wav[idx], flux[idx], err[idx]

    mascara = np.isfinite(flux) & np.isfinite(wav) & np.isfinite(err)
    return wav[mascara], flux[mascara], err[mascara]



def fotometria_sintetica(wav_spec, flux_spec, wav_filt, trans_filt):

    wav_min = max(wav_spec.min(), wav_filt.min())
    wav_max = min(wav_spec.max(), wav_filt.max())

    if wav_min >= wav_max:
        raise ValueError("El espectro y el filtro no se solapan. "
                         "Revisa el rango de longitudes de onda.")

    mascara_spec = (wav_spec >= wav_min) & (wav_spec <= wav_max)
    mascara_filt = (wav_filt >= wav_min) & (wav_filt <= wav_max)

    wav_comun = wav_spec[mascara_spec]   # Usamos la grilla del espectro

    interp_filt = interpolate.interp1d(
        wav_filt[mascara_filt], trans_filt[mascara_filt],
        kind="linear", bounds_error=False, fill_value=0.0
    )
    trans_en_spec = interp_filt(wav_comun)

    flux_en_spec = flux_spec[mascara_spec]

    # Integración numérica (trapecio)
    numerador   = np.trapezoid(flux_en_spec * trans_en_spec, wav_comun)
    denominador = np.trapezoid(trans_en_spec, wav_comun)

    if denominador == 0:
        raise ValueError("La integral de transmisión es cero. "
                         "Verifica la curva del filtro.")

    flux_sintetico = numerador / denominador
    return flux_sintetico, wav_comun, flux_en_spec, trans_en_spec

def planck_Jy(wav_um, T_K, R_cm, D_kpc):

    wav_cm = wav_um * 1e-4           # µm → cm
    nu     = c / wav_cm              # Hz

    # Función de Planck B_ν [erg/s/cm²/Hz/sr]
    B_nu = (2 * h * nu**3 / c**2) / (np.exp(h * nu / (kB * T_K)) - 1)

    D_cm = D_kpc * KPC_CM
    flujo_cgs = np.pi * B_nu * (R_cm / D_cm)**2   # erg/s/cm²/Hz

    flujo_Jy = flujo_cgs * JY_PER_CGS             # Jy
    return flujo_Jy



def calcular_W(flux_target_Jy, T_ref, R_ref_cm, D_kpc, wav_ref_um=4.6):

    flux_ref_Jy = planck_Jy(wav_ref_um, T_ref, R_ref_cm, D_kpc)
    W = flux_target_Jy / flux_ref_Jy
    return W, flux_ref_Jy


def graficar(wav_spec, flux_spec,
             wav_filt, trans_filt,
             wav_comun, flux_comun, trans_comun,
             flux_sintetico):


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8),
                                   gridspec_kw={"height_ratios": [3, 1]})
    fig.suptitle("Fotometría Sintética – IRAC 4.5 µm", fontsize=14, fontweight="bold")

    # Panel superior: espectro + punto fotométrico
    ax1.plot(wav_spec, flux_spec, color="steelblue", lw=0.8,
             label="Espectro ISO SWS+LWS", zorder=2)
    ax1.axvspan(FILTRO_MIN_UM, FILTRO_MAX_UM, alpha=0.15, color="orange",
                label="Región del filtro IRAC 4.5 µm")
    ax1.plot(4.5, flux_sintetico, "r*", ms=14, zorder=5,
             label=f"Flujo sintético = {flux_sintetico:.3f} Jy")
    ax1.set_xlim(1, 30)
    ax1.set_ylabel("Flujo (Jy)", fontsize=11)
    ax1.set_yscale("linear")
    ax1.legend(fontsize=9)
    ax1.grid(True, ls="--", alpha=0.4)
    ax1.set_title("Espectro completo", fontsize=10)

    ax2b = ax2.twinx()
    ax2.plot(wav_comun, flux_comun, color="steelblue", lw=1.2,
             label="Flujo en región del filtro")
    ax2b.fill_between(wav_filt, trans_filt, alpha=0.3, color="orange")
    ax2b.plot(wav_filt, trans_filt, color="darkorange", lw=1.5,
              label="Transmisión IRAC 4.5 µm")
    ax2.set_xlabel("Longitud de onda (µm)", fontsize=11)
    ax2.set_ylabel("Flujo (Jy)", fontsize=10, color="steelblue")
    ax2b.set_ylabel("Transmisión", fontsize=10, color="darkorange")
    ax2b.set_ylim(0, 1.2)
    ax2.set_xlim(FILTRO_MIN_UM, FILTRO_MAX_UM)
    ax2.grid(True, ls="--", alpha=0.4)

    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper right")

    plt.tight_layout()
    plt.savefig("fotometria_sintetica.png", dpi=150, bbox_inches="tight")
    print("[Figura] Guardada como 'fotometria_sintetica.png'")
    plt.show()



def main():
    print("=" * 65)
    print("  FOTOMETRÍA SINTÉTICA ISO → Filtro IRAC 4.5 µm")
    print("=" * 65)

    wav_filt, trans_filt = cargar_filtro_irac45(ARCHIVO_FILTRO)

    print(f"\n[Espectros] Cargando SWS: {ARCHIVO_SWS}")
    print(f"[Espectros] Cargando LWS: {ARCHIVO_LWS}")
    wav_spec, flux_spec, err_spec = unir_sws_lws(ARCHIVO_SWS, ARCHIVO_LWS)
    print(f"[Espectros] Rango combinado: {wav_spec.min():.2f} – "
          f"{wav_spec.max():.2f} µm  | {len(wav_spec)} puntos")

    print("\n[Paso 3] Calculando flujo sintético en filtro IRAC 4.5 µm...")
    flux_target, wav_comun, flux_comun, trans_comun = fotometria_sintetica(
        wav_spec, flux_spec, wav_filt, trans_filt
    )
    print(f"  → Flux_target = {flux_target:.4f} Jy")
    print(f"     (Cerrigone reporta ≈ 11.38 Jy para su fuente)")

    print(f"\n[Paso 4] Calculando flujo del cuerpo negro de referencia...")
    print(f"  T = {T_REF} K | R = {R_REF_CM:.1e} cm | D = {DISTANCIA_KPC} kpc")
    W, flux_ref = calcular_W(flux_target, T_REF, R_REF_CM, DISTANCIA_KPC)
    print(f"  → Flux_ref   = {flux_ref:.4e} Jy  (a 4.6 µm)")

    print(f"\n[Paso 5] Calculando W...")
    print(f"  → W = Flux_target / Flux_ref = {flux_target:.4f} / {flux_ref:.4e}")
    print(f"  → W = {W:.4f}")

    print("\n" + "=" * 65)
    print("  RESUMEN FINAL")
    print("=" * 65)
    print(f"  Flux sintético IRAC 4.5 µm : {flux_target:.4f}  Jy")
    print(f"  Flux cuerpo negro ref.      : {flux_ref:.4e}  Jy")
    print(f"  Distancia usada             : {DISTANCIA_KPC}  kpc")
    print(f"  W                           : {W:.4f}")
    print("=" * 65)

    print("\n[Figura] Generando gráfica...")
    graficar(wav_spec, flux_spec,
             wav_filt, trans_filt,
             wav_comun, flux_comun, trans_comun,
             flux_target)

    return flux_target, flux_ref, W


if __name__ == "__main__":
    flux_target, flux_ref, W = main()