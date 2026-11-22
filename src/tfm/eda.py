"""
Análisis exploratorio (EDA) por donante y control de calidad por célula.

Regla de oro: las variables del donante se cuentan POR DONANTE (una fila por donante), nunca por célula; si no,
un donante con muchas células pesaría más.

    ordinal(serie)                       Braak/Thal/CERAD/intervalos -> número ordenado
    exact_age(development_stage)         "88-year-old stage" -> 88 (y marca los "80 and over")
    association_table(donantes, ...)     contraste de cada variable frente al ADNC y correlación con el CPS
    mad_outliers(obs, ...)               células atípicas por tipo celular (k MAD en escala log)
    donor_medians(obs, ...)              mediana por (donante, tipo) de las métricas de QC
    composition(obs, ...)                proporción de cada subtipo dentro de cada tipo celular, por donante
"""
import re

import numpy as np
import pandas as pd
from scipy import stats

ROMANOS = {"0": 0, "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}
CERAD = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}
NO_DATO = {"Reference", "nan", "None", "", "Unknown", "unknown"}


def ordinal(s):
    """Convierte categorías ordenadas a números. Devuelve NaN si no sabe ordenar el valor.

    'Braak IV' -> 4, 'Thal 3' -> 3, 'Moderate' -> 2, '78 to 89 years old' -> 78, '3.2 to 5.9 hours' -> 3.2.
    """
    def uno(v):
        v = str(v).strip()
        if v in NO_DATO:
            return np.nan
        if v in CERAD:
            return CERAD[v]
        if v == "Not Identified":                       # LATE-NC sin patología
            return 0.0
        m = re.match(r"^LATE Stage (\d)$", v)
        if m:
            return float(m.group(1))
        m = re.match(r"^(?:Braak|Thal)\s+([0IVX\d]+)$", v)
        if m:
            x = m.group(1)
            return float(ROMANOS.get(x, x)) if x in ROMANOS or x.isdigit() else np.nan
        m = re.match(r"^(?:Less than\s+)?(\d+(?:\.\d+)?)", v)
        return float(m.group(1)) if m else np.nan
    return s.astype(str).map(uno).astype(float)


def exact_age(dev_stage):
    """'88-year-old stage' -> 88. '80 year-old and over stage' -> 80 con marca de censura."""
    txt = dev_stage.astype(str)
    edad = txt.str.extract(r"(\d+)")[0].astype(float)
    censurada = txt.str.contains("and over", case=False)
    return edad, censurada


def _cramers_v(tabla):
    chi2 = stats.chi2_contingency(tabla, correction=False)[0]
    n = tabla.values.sum()
    k = min(tabla.shape) - 1
    return float(np.sqrt(chi2 / (n * k))) if n and k else np.nan


def association_table(donantes, columnas, label, orden, cps=None):
    """Una fila por variable: tipo de contraste frente al ADNC (sin Reference), estadístico, p y tamaño del efecto,
    y correlación de Spearman con el CPS si la variable se puede ordenar.

    - Variables con orden (numéricas o ordinales): Kruskal-Wallis entre clases + Spearman con la clase ordinal.
    - Variables sin orden (sexo, APOE4, etnia...): chi-cuadrado y V de Cramér. Con pocos donantes por casilla el
      p-valor del chi-cuadrado es aproximado; se indica la casilla esperada mínima.
    """
    d = donantes[donantes[label].isin(orden)]
    rango = d[label].map({c: i for i, c in enumerate(orden)})
    filas = []
    for c in columnas:
        if c not in d.columns:
            continue
        num = pd.to_numeric(d[c], errors="coerce") if pd.api.types.is_numeric_dtype(d[c]) else ordinal(d[c])
        fila = {"variable": c, "n_donantes": int(d[c].astype(str).map(lambda v: v not in NO_DATO).sum())}
        if num.notna().sum() >= 0.8 * len(d) and num.nunique() > 2:
            grupos = [num[d[label] == k].dropna() for k in orden]
            h, p = stats.kruskal(*grupos)
            rho, prho = stats.spearmanr(rango[num.notna()], num.dropna())
            fila.update(tipo="ordenada", contraste="Kruskal-Wallis", estadistico=h, p=p,
                        efecto=f"rho con clase = {rho:.2f}", rho_clase=rho, p_rho_clase=prho)
            if cps is not None:
                r2, p2 = stats.spearmanr(cps.reindex(d.index)[num.notna()], num.dropna())
                fila.update(rho_cps=r2, p_rho_cps=p2)
        else:
            valido = d[c].astype(str).map(lambda v: v not in NO_DATO)
            tabla = pd.crosstab(d.loc[valido, label], d.loc[valido, c].astype(str))
            if tabla.shape[1] < 2:
                continue
            chi2, p, _, esperado = stats.chi2_contingency(tabla, correction=False)
            fila.update(tipo="categórica", contraste="chi-cuadrado", estadistico=chi2, p=p,
                        efecto=f"V de Cramér = {_cramers_v(tabla):.2f}", v_cramer=_cramers_v(tabla),
                        esperado_min=float(esperado.min()))
        filas.append(fila)
    return pd.DataFrame(filas).set_index("variable")


def mad_outliers(obs, metricas, group, k=3.0):
    """Marca como atípicas las células a más de k MAD de la mediana de su tipo celular, en escala log1p.

    Es el criterio adaptativo habitual en single-cell: cada tipo celular tiene su propio rango normal (las neuronas
    tienen muchos más recuentos que la glía), así que un umbral fijo para todos no sirve.
    Devuelve un DataFrame booleano células x métricas (True = atípica) y una tabla con los límites usados.
    """
    marcas, limites = {}, []
    for m in metricas:
        x = np.log1p(obs[m].astype(float))
        med = x.groupby(obs[group], observed=True).transform("median")
        mad = (x - med).abs().groupby(obs[group], observed=True).transform("median") * 1.4826
        marcas[m] = (x - med).abs() > k * mad
        for t, g in x.groupby(obs[group], observed=True):
            md = g.median()
            md_ = (g - md).abs().median() * 1.4826
            limites.append({"tipo": t, "metrica": m, "minimo": float(np.expm1(md - k * md_)),
                            "maximo": float(np.expm1(md + k * md_))})
    return pd.DataFrame(marcas, index=obs.index), pd.DataFrame(limites)


def donor_medians(obs, metricas, donor, group):
    """Mediana de cada métrica por (donante, tipo celular)."""
    return obs.groupby([donor, group], observed=True)[list(metricas)].median().reset_index()


def composition(obs, donor, group, sub):
    """Proporción de cada `sub` (p. ej. Supertype) dentro de cada `group`, por donante. Formato largo."""
    n = obs.groupby([donor, group, sub], observed=True).size().rename("n").reset_index()
    n["proporcion"] = n["n"] / n.groupby([donor, group], observed=True)["n"].transform("sum")
    return n
