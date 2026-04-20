import numpy as np
import pandas as pd

stars = {
    'IRAS_06531': {
        'name':                     'IRAS_06531',
        'CO(1-0) Jy·km/s':         27.9920 ,
        '13CO(1-0) Jy·km/s':       None   ,
        'CO(2-1) Jy·km/s':         87.8595 ,
        '13CO(2-1) Jy·km/s':       1.0361 ,
        'Ratio R(2-1)/(1-0)':      3.1387,
        'Ratio R(2-1)/13CO(2-1)':  84.6838,
        'Ratio R(1-0)/13CO(1-0)':  None,
        'type':                     'C star',
    },
    'IRAS_07145': {
        'name':                     'IRAS_07145',
        'CO(1-0) Jy·km/s':         10.7001 ,
        '13CO(1-0) Jy·km/s':       None  ,
        'CO(2-1) Jy·km/s':         38.4563 ,
        '13CO(2-1) Jy·km/s':       1.4299 ,
        'Ratio R(2-1)/(1-0)':      3.5940,
        'Ratio R(2-1)/13CO(2-1)':  26.7066,
        'Ratio R(1-0)/13CO(1-0)':  None,
        'type':                     'C star',
    },
    'IRAS_07180': {
        'name':                     'IRAS_07180',
        'CO(1-0) Jy·km/s':         1.0541 ,
        '13CO(1-0) Jy·km/s':       None   ,
        'CO(2-1) Jy·km/s':         11.3168 ,
        '13CO(2-1) Jy·km/s':       2.2871 ,
        'Ratio R(2-1)/(1-0)':      10.7362,
        'Ratio R(2-1)/13CO(2-1)':  4.9481,
        'Ratio R(1-0)/13CO(1-0)':  None,
        'type':                     'OH/IR star',
    },
    'IRAS_08500': {
        'name':                     'IRAS_08500',
        'CO(1-0) Jy·km/s':         44.9409 ,
        '13CO(1-0) Jy·km/s':       None   ,
        'CO(2-1) Jy·km/s':         158.2523 ,
        '13CO(2-1) Jy·km/s':       35.8345 ,
        'Ratio R(2-1)/(1-0)':      3.5213,
        'Ratio R(2-1)/13CO(2-1)':  4.4162,
        'Ratio R(1-0)/13CO(1-0)':  None,
        'type':                     'Mira M',
    },
    'RV_Aqr': {
        'name':                     'RV_Aqr',
        'CO(1-0) Jy·km/s':         197.9703 ,
        '13CO(1-0) Jy·km/s':       7.3052   ,
        'CO(2-1) Jy·km/s':         None,
        '13CO(2-1) Jy·km/s':       None,
        'Ratio R(2-1)/(1-0)':      None,
        'Ratio R(2-1)/13CO(2-1)':  None,
        'Ratio R(1-0)/13CO(1-0)':  27.0999,
        'type':                     'C star',
    },
    'RZ_Sgr': {
        'name':                     'RZ_Sgr',
        'CO(1-0) Jy·km/s':         70.5720  ,
        '13CO(1-0) Jy·km/s':       15.1140   ,
        'CO(2-1) Jy·km/s':         None,
        '13CO(2-1) Jy·km/s':       None,
        'Ratio R(2-1)/(1-0)':      None,
        'Ratio R(2-1)/13CO(2-1)':  None,
        'Ratio R(1-0)/13CO(1-0)':  4.6693,
        'type':                     'S star',
    },
    'UU_For': {
        'name':                     'UU_For',
        'CO(1-0) Jy·km/s':         63.8513   ,
        '13CO(1-0) Jy·km/s':       3.3920    ,
        'CO(2-1) Jy·km/s':         None,
        '13CO(2-1) Jy·km/s':       None,
        'Ratio R(2-1)/(1-0)':      None,
        'Ratio R(2-1)/13CO(2-1)':  None,
        'Ratio R(1-0)/13CO(1-0)':  18.8238,
        'type':                     'Mira M',
    },
    'R_Hor': {
        'name':                     'R_Hor',
        'CO(1-0) Jy·km/s':         70.4809   ,
        '13CO(1-0) Jy·km/s':       1.8016    ,
        'CO(2-1) Jy·km/s':         None,
        '13CO(2-1) Jy·km/s':       None,
        'Ratio R(2-1)/(1-0)':      None,
        'Ratio R(2-1)/13CO(2-1)':  None,
        'Ratio R(1-0)/13CO(1-0)':  39.1210,
        'type':                     'Mira M',
    },
    'V510_Pup': {
        'name':                     'V510_Pup',
        'CO(1-0) Jy·km/s':         15.6387   ,
        '13CO(1-0) Jy·km/s':       None    ,
        'CO(2-1) Jy·km/s':         173.9363 ,
        '13CO(2-1) Jy·km/s':       14.3408 ,
        'Ratio R(2-1)/(1-0)':      11.1222,
        'Ratio R(2-1)/13CO(2-1)':  12.1288,
        'Ratio R(1-0)/13CO(1-0)':  None,
        'type':                     'Post-AGB F5Ie',
    },
}

# ============================================================
# CONSTRUIR TABLA (solo fuentes con datos completos)
# ============================================================
rows = []
for name, params in stars.items():
    rows.append({
        'Source':                   params['name'],
        'Type':                     params['type'],
        'CO(1-0) Jy·km/s':         params['CO(1-0) Jy·km/s'],
        '13CO(1-0) Jy·km/s':       params['13CO(1-0) Jy·km/s'],
        'CO(2-1) Jy·km/s':         params['CO(2-1) Jy·km/s'],
        '13CO(2-1) Jy·km/s':       params['13CO(2-1) Jy·km/s'],
        'Ratio R(2-1)/(1-0)':      params['Ratio R(2-1)/(1-0)'],
        'Ratio R(2-1)/13CO(2-1)':  params['Ratio R(2-1)/13CO(2-1)'],
        'Ratio R(1-0)/13CO(1-0)':  params['Ratio R(1-0)/13CO(1-0)'],
    })

df = pd.DataFrame(rows)
print(df.to_string(index=False))

df.to_csv('line_ratios.csv', index=False)
print("\nGuardado: line_ratios.csv")