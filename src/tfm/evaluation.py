"""
Evaluación con validación cruzada por donante.

La unidad de análisis es el donante: cada donante aparece entero en train o en test, nunca en ambos.

    donor_splits(y, n_splits, n_repeats, seed)    particiones estratificadas por clase, repetidas
    cross_validate(X, y, make_model, splits)      métricas por repetición y predicciones fuera de pliegue
    technical_features(pb, cfg)                   matriz donante x (log n_celulas, log total_counts) por tipo
    summarize(resultados)                         media y desviación de cada métrica

Métricas (las clases están desequilibradas: Low 12, Intermediate 21, High 42):
    balanced_accuracy   media del acierto por clase; el azar da 1/3
    f1_macro            media del F1 de cada clase
    kappa_cuadratico    acuerdo ordinal: confundir Low con High penaliza más que con Intermediate; el azar da 0
"""
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, cohen_kappa_score, f1_score
from sklearn.model_selection import RepeatedStratifiedKFold


def donor_splits(y, n_splits=5, n_repeats=10, seed=0):
    """Lista de (repetición, train_idx, test_idx). `y` tiene una etiqueta por donante."""
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    return [(i // n_splits, tr, te) for i, (tr, te) in enumerate(rskf.split(np.zeros(len(y)), y))]


def _metrics(y_true, y_pred, orden):
    rank = {c: i for i, c in enumerate(orden)}
    return {
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "f1_macro": f1_score(y_true, y_pred, average="macro", labels=orden, zero_division=0),
        "kappa_cuadratico": cohen_kappa_score([rank[v] for v in y_true], [rank[v] for v in y_pred],
                                              weights="quadratic"),
    }


def cross_validate(X, y, make_model, splits, orden):
    """Entrena un modelo nuevo en cada pliegue. Las métricas se calculan por repetición, juntando las
    predicciones de sus pliegues (cada donante se predice una vez por repetición)."""
    X = np.asarray(X)
    y = np.asarray(y)
    pred = {}
    for rep, tr, te in splits:
        modelo = make_model()
        modelo.fit(X[tr], y[tr])
        pred.setdefault(rep, np.empty(len(y), dtype=object))[te] = modelo.predict(X[te])
    filas = [{"repeticion": rep, **_metrics(y, p, orden)} for rep, p in sorted(pred.items())]
    return pd.DataFrame(filas), pred


def summarize(resultados):
    """{nombre: DataFrame de cross_validate} -> tabla con media ± desviación de cada métrica."""
    filas = []
    for nombre, df in resultados.items():
        fila = {"modelo": nombre}
        for m in ["balanced_accuracy", "f1_macro", "kappa_cuadratico"]:
            fila[f"{m}_media"] = df[m].mean()
            fila[f"{m}_sd"] = df[m].std()
        filas.append(fila)
    return pd.DataFrame(filas).set_index("modelo")


def partial_spearman(x, y, Z):
    """Correlación de Spearman parcial entre x e y controlando por las columnas de Z.

    Se pasa todo a rangos, se quita a x y a y la parte que explica Z (regresión lineal sobre los rangos) y se
    correlacionan los residuos. Devuelve (rho, p) con n - 2 - k grados de libertad.
    """
    from scipy import stats

    x, y = stats.rankdata(np.asarray(x, float)), stats.rankdata(np.asarray(y, float))
    Z = np.asarray(Z, float).reshape(len(x), -1)
    Zr = np.column_stack([np.ones(len(x))] + [stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    rx = x - Zr @ np.linalg.lstsq(Zr, x, rcond=None)[0]
    ry = y - Zr @ np.linalg.lstsq(Zr, y, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    gl = len(x) - 2 - Z.shape[1]
    t = r * np.sqrt(gl / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), gl))


def donor_covariates(obs, cfg, mt_col="pct_counts_mt"):
    """Covariables biológicas y de calidad por donante (sin Reference), para controlarlas:

    apoe4 (1 = portador), sexo_mujer (1 = mujer), edad (exacta; los "80 o más" valen 80) y edad_80_o_mas
    (1 = dato censurado), y la mediana del % mitocondrial de cada tipo celular (mt_<tipo>).
    """
    from .data import get_donor_table
    from .eda import exact_age

    C = cfg["columns"]
    don, lab, ct = C["donor"], C["label"], C["cell_type"]
    d = get_donor_table(obs, don)
    d = d[d[lab].isin(cfg["labels"]["classification"])]
    edad, cens = exact_age(d["development_stage"])
    cov = pd.DataFrame({
        "apoe4": (d["APOE4 status"].astype(str) == "Y").astype(float),
        "sexo_mujer": (d["sex"].astype(str) == "female").astype(float),
        "edad": edad.clip(upper=80).astype(float),
        "edad_80_o_mas": cens.astype(float),
    }, index=d.index.astype(str))
    mt = obs.groupby([obs[don].astype(str), ct], observed=True)[mt_col].median().unstack(ct)
    mt.columns = [f"mt_{c}" for c in mt.columns]
    mt = mt.reindex(cov.index)
    cov = cov.join(mt.fillna(mt.median()))
    return cov


def technical_features(pb, cfg, columnas=("n_celulas", "total_counts")):
    """Una fila por donante (sin Reference) con log10 de cada variable técnica en cada tipo celular.
    Los pares que faltan se rellenan con 0 (log10(1)). Devuelve (X, y)."""
    C = cfg["columns"]
    don, lab, ct = C["donor"], C["label"], C["cell_type"]
    obs = pb.obs[pb.obs[lab].isin(cfg["labels"]["classification"])]
    partes = []
    for c in columnas:
        t = obs.pivot_table(index=don, columns=ct, values=c, aggfunc="sum", observed=True)
        t = np.log10(t.fillna(0).clip(lower=1))
        t.columns = [f"{c}_{x}" for x in t.columns]
        partes.append(t)
    X = pd.concat(partes, axis=1)
    y = obs.groupby(don, observed=True)[lab].first().astype(str).reindex(X.index)
    return X, y
