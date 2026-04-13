import numpy as np
import pandas as pd

stars = {
    'IRAS_06531': {
        'T_mb_peak': 0.329,  
        'V_exp':     27.312,  
        'D':         1158.35, 
        'BMAJ_deg':  4.972155506402E-03,
        'BMIN_deg':  2.351931261700E-03,
        'W':         0.7572,
        's_J':       1.0,
        'f_CO':      1e-3,    
        'type':      'C star',
        'transition':'CO(1-0)'
    },
    'IRAS_07145': {
        'T_mb_peak': 0.182,  
        'V_exp':     28.618,    
        'D':         4144.8,   
        'BMAJ_deg':  3.618474997656E-03,
        'BMIN_deg':  2.157910122644E-03,
        'W':         7.6518,
        's_J':       1.0,
        'f_CO':      1e-3,     
        'type':      'C star',
        'transition':'CO(1-0)'
    },
    'IRAS_07180': {
        'T_mb_peak': 0.058,   
        'V_exp':     15.547,    
        'D':         1168.22,   
        'BMAJ_deg':  3.618474997656E-03,
        'BMIN_deg':  2.157910122644E-03,
        'W':         0.5282,
        's_J':       1.0,
        'f_CO':      1e-3,     
        'type':      'OH/IR star',
        'transition':'CO(1-0)'
    },
    'RV_Aqr': {
        'T_mb_peak': 2.5315,   
        'V_exp':     14.92,    
        'D':         635.24,   
        'BMAJ_deg':  4.608768552314E-03,
        'BMIN_deg':  2.453349310141E-03,
        'W':         0.3586,
        's_J':       1.0,
        'f_CO':      1e-3,     
        'type':      'C star',
        'transition':'CO(1-0)'
    },
    'RZ_Sgr': {
        'T_mb_peak': 2.534,     
        'V_exp':     6.506,
        'D':         730,
        'BMAJ_deg':  3.718631930136E-03,
        'BMIN_deg':  2.511840595147E-03,
        'W':         0.5680,
        's_J':       1.0,
        'f_CO':      6e-4,    
        'type':      'S star',
        'transition':'CO(1-0)'
    },
    'UU_For': {
        'T_mb_peak': 1.175,
        'V_exp':     13.391,
        'D':         614.76,
        'BMAJ_deg':  4.293047298298E-03,
        'BMIN_deg':  2.455889138257E-03,
        'W':         0.1662,
        's_J':       1.0,
        'f_CO':      2e-4,     
        'type':      'Mira M',
        'transition':'CO(1-0)'
    },
    'R_Hor': {
        'T_mb_peak': 3.561,
        'V_exp':     4.619,
        'D':         231.29,
        'BMAJ_deg':  3.667337556605E-03,
        'BMIN_deg':  2.591346856427E-03,
        'W':         0.0432,
        's_J':       1.0,
        'f_CO':      2e-4,     
        'type':      'Mira M',
        'transition':'CO(1-0)'
    },
    'IRAS_08500': {
        'T_mb_peak': 1.064,
        'V_exp':     12.476,
        'D':         1526.72,
        'BMAJ_deg':  4.608768552314E-03,
        'BMIN_deg':  2.453349310141E-03,
        'W':         1.0670,
        's_J':       1.0,
        'f_CO':      2e-4,     
        'type':      'Mira M',
        'transition':'CO(1-0)'
    },
    'V510_Pup': {
        'T_mb_peak': 0.074,   
        'V_exp':     152.250,    
        'D':         436.26,   
        'BMAJ_deg':  4.571202342152E-03,
        'BMIN_deg':  2.611943497287E-03,
        'W':         0.0431,
        's_J':       1.0,
        'f_CO':      4e-4,     
        'type':      'Post-AGB F5Ie',
        'transition':'CO(1-0)'
    },
}

# ============================================================
# CÁLCULO
# ============================================================
def calc_mdot(params):
    T_mb   = params['T_mb_peak']
    V_exp  = params['V_exp']
    D      = params['D']
    BMAJ   = params['BMAJ_deg'] * 3600  # deg -> arcsec
    BMIN   = params['BMIN_deg'] * 3600
    theta  = np.sqrt(BMAJ * BMIN)
    W      = params['W']
    s_J    = params['s_J']
    f_CO   = params['f_CO']

    log_term = np.log10(W / 0.04)
    T_term   = T_mb / (log_term * s_J)
    Dtheta   = D * theta

    Mdot_thin = (4.55e-19
                 * (T_term ** (5/6))
                 * (f_CO ** (-1))
                 * (V_exp ** (11/6))
                 * (Dtheta ** (5/3)))

    Mdot_thick = (1.4 * T_mb * (V_exp**2) * (D**2) * (theta**2)
                  / (2e19 * (f_CO**0.85) * s_J))

    return Mdot_thin, Mdot_thick, theta, Dtheta, log_term

# ============================================================
# CONSTRUIR TABLA
# ============================================================
rows = []
for name, params in stars.items():
    # saltar si faltan datos
    if any(v is None for v in params.values()):
        print(f"  {name}: faltan parámetros, saltando...")
        continue

    Mdot_thin, Mdot_thick, theta, Dtheta, logW = calc_mdot(params)

    rows.append({
        'Source':          name,
        'Type':            params['type'],
        'Transition':      params['transition'],
        'D_pc':            params['D'],
        'T_mb_peak_K':     params['T_mb_peak'],
        'V_exp_kms':       params['V_exp'],
        'theta_arcsec':    round(theta, 4),
        'D_theta':         round(Dtheta, 2),
        'W':               round(params['W'], 4),
        'log_W_0.04':      round(logW, 4),
        'f_CO':            params['f_CO'],
        's_J':             params['s_J'],
        'Mdot_thin_Msunyr':  f'{Mdot_thin:.4e}',
        'Mdot_thick_Msunyr': f'{Mdot_thick:.4e}',
    })

df = pd.DataFrame(rows)
print(df.to_string(index=False))

# ============================================================
# EXPORTAR
# ============================================================
df.to_csv('mass_loss_CO10.csv', index=False)
print("\nGuardado: mass_loss_CO10.csv")