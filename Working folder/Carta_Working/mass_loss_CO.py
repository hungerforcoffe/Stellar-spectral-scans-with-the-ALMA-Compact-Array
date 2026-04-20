import numpy as np
import pandas as pd

stars = {
    'IRAS_06531_10': {
        'name':      'IRAS_06531',
        'T_mb_peak': 0.3379,  
        'V_exp':     26.67,  
        'D':         1158.35, 
        'BMAJ_deg':  4.972155506402E-03,
        'BMIN_deg':  2.351931261700E-03,
        'W':         0.7572,
        's_J':       1.0,
        'f_CO':      1e-3,    
        'type':      'C star',
        'transition':'CO(1-0)'
    },
    'IRAS_06531_21': {
        'name':      'IRAS_06531',
        'T_mb_peak': 1.4107,  
        'V_exp':     25.55,  
        'D':         1158.35, 
        'BMAJ_deg':  1.809633214088E-03,
        'BMIN_deg':  1.284642819447E-03,
        'W':         0.7572,
        's_J':       0.6,
        'f_CO':      1e-3,    
        'type':      'C star',
        'transition':'CO(2-1)'
    },
    'IRAS_07145_10': {
        'name':      'IRAS_07145',
        'T_mb_peak': 0.1788,  
        'V_exp':     28.13,    
        'D':         4144.8,   
        'BMAJ_deg':  3.618474997656E-03,
        'BMIN_deg':  2.157910122644E-03,
        'W':         7.6518,
        's_J':       1.0,
        'f_CO':      1e-3,     
        'type':      'C star',
        'transition':'CO(1-0)'
    },
    'IRAS_07145_21': {
        'name':      'IRAS_07145',
        'T_mb_peak': 0.5994,  
        'V_exp':     28.20,    
        'D':         4144.8,   
        'BMAJ_deg':  1.945121659795E-03,
        'BMIN_deg':  1.168968128587E-03,
        'W':         7.6518,
        's_J':       0.6,
        'f_CO':      1e-3,     
        'type':      'C star',
        'transition':'CO(2-1)'
    },
    'IRAS_07180_10': {
        'name':      'IRAS_07180',
        'T_mb_peak': 0.0424,   
        'V_exp':     20.87,    
        'D':         1168.22,   
        'BMAJ_deg':  3.618474997656E-03,
        'BMIN_deg':  2.157910122644E-03,
        'W':         0.5282,
        's_J':       1.0,
        'f_CO':      1e-3,     
        'type':      'OH/IR star',
        'transition':'CO(1-0)'
    },
    'IRAS_07180_21': {
        'name':      'IRAS_07180',
        'T_mb_peak': 0.3485,   
        'V_exp':     12.73,    
        'D':         1168.22,   
        'BMAJ_deg':  1.945121659795E-03,
        'BMIN_deg':  1.168968128587E-03,
        'W':         0.5282,
        's_J':       0.6,
        'f_CO':      1e-3,     
        'type':      'OH/IR star',
        'transition':'CO(2-1)'
    },
    'IRAS_08500_10': {
        'name':      'IRAS_08500',
        'T_mb_peak': 1.2623,
        'V_exp':     11.08,
        'D':         1526.72,
        'BMAJ_deg':  4.608768552314E-03,
        'BMIN_deg':  2.453349310141E-03,
        'W':         1.0670,
        's_J':       1.0,
        'f_CO':      2e-4,     
        'type':      'Mira M',
        'transition':'CO(1-0)'
    },
    'IRAS_08500_21': {
        'name':      'IRAS_08500',
        'T_mb_peak': 1.2623,
        'V_exp':     11.08,
        'D':         1526.72,
        'BMAJ_deg':  2.031578940321E-03,
        'BMIN_deg':  1.273627257867E-03,
        'W':         1.0670,
        's_J':       0.6,
        'f_CO':      2e-4,     
        'type':      'Mira M',
        'transition':'CO(2-1)'
    },
    'RV_Aqr': {
        'name':      'RV_Aqr',
        'T_mb_peak': 4.7257,   
        'V_exp':     12.80,    
        'D':         635.24,   
        'BMAJ_deg':  4.608768552314E-03,
        'BMIN_deg':  2.453349310141E-03,
        'W':         0.3586,
        's_J':       1.0,
        'f_CO':      2e-4,     
        'type':      'C star',
        'transition':'CO(1-0)'
    },
    'RZ_Sgr': {
        'name':      'RZ_Sgr',
        'T_mb_peak': 4.8113,     
        'V_exp':     6.42,
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
        'name':      'UU_For',
        'T_mb_peak': 1.7344,
        'V_exp':     12.74,
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
        'name':      'R_Hor',
        'T_mb_peak': 6.4034,
        'V_exp':     6.99,
        'D':         231.29,
        'BMAJ_deg':  3.667337556605E-03,
        'BMIN_deg':  2.591346856427E-03,
        'W':         0.0432,
        's_J':       1.0,
        'f_CO':      2e-4,     
        'type':      'Mira M',
        'transition':'CO(1-0)'
    },
    'V510_Pup_10': {
        'name':      'V510_Pup',
        'T_mb_peak': 0.0515,   
        'V_exp':     101.54,    
        'D':         436.26,   
        'BMAJ_deg':  4.571202342152E-03,
        'BMIN_deg':  2.611943497287E-03,
        'W':         0.0431,
        's_J':       1.0,
        'f_CO':      4e-4,     
        'type':      'Post-AGB F5Ie',
        'transition':'CO(1-0)'
    },
    'V510_Pup_21': {
        'name':      'V510_Pup',
        'T_mb_peak': 0.6814,   
        'V_exp':     178.06,    
        'D':         436.26,   
        'BMAJ_deg':  1.941144005032E-03,
        'BMIN_deg':  1.210112294358E-03, 
        'W':         0.0431,
        's_J':       0.6,
        'f_CO':      4e-4,     
        'type':      'Post-AGB F5Ie',
        'transition':'CO(2-1)'
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
        'Source':          params['name'],
        'Type':            params['type'],
        'Transition':      params['transition'],
        'D_pc':            params['D'],
        'T_mb_peak_K':     params['T_mb_peak'],
        'V_exp_kms':       params['V_exp'],
        'theta_arcsec':    round(theta, 4),
        'W':               round(params['W'], 4),
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
df.to_csv('mass_loss_CO10_CO21.csv', index=False)
print("\nGuardado: mass_loss_CO10_CO21.csv")