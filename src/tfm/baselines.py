"""
Baselines: predecir el ADNC del donante a partir de sus pseudobulk.

Diseño
    - Un modelo por tipo celular, con los donantes que tienen ese tipo (tras igualar profundidad no todos lo tienen).
    - Fusión tardía: para cada donante se promedian las probabilidades de los modelos de los tipos que tiene.
    - Todo lo que mira los datos (genes expresados, genes más variables, escalado) se ajusta solo con los
      donantes de entrenamiento de cada pliegue, dentro de un Pipeline de scikit-learn.

    LogCPMTopVar                 transformador: recuentos -> log-CPM de los genes expresados más variables
    build_type_data(pb, ...)     por tipo celular: recuentos y variables técnicas, una fila por donante
    run_baselines(...)           validación cruzada por donante de varios modelos en cada tipo
    fuse(...)                    fusión tardía de las probabilidades por donante
"""
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.base import BaseEstimator, TransformerMixin

from .evaluation import _metrics


class LogCPMTopVar(BaseEstimator, TransformerMixin):
    """Recibe [recuentos | n_tech columnas técnicas] y devuelve [log-CPM de n_genes genes | técnicas].

    En fit (solo con los donantes de entrenamiento):
        1. CPM y genes expresados: CPM >= min_cpm en al menos min_frac de las muestras.
        2. De esos, los n_genes con mayor varianza de log-CPM (no usa la etiqueta).
    """

    def __init__(self, n_genes=2000, min_cpm=1.0, min_frac=0.5, n_tech=0):
        self.n_genes = n_genes
        self.min_cpm = min_cpm
        self.min_frac = min_frac
        self.n_tech = n_tech

    def _split(self, X):
        X = np.asarray(X, dtype=np.float64)
        if self.n_tech:
            return X[:, :-self.n_tech], X[:, -self.n_tech:]
        return X, None

    @staticmethod
    def _logcpm(C):
        tot = C.sum(1, keepdims=True)
        tot[tot == 0] = 1
        return np.log1p(C / tot * 1e6), C / tot * 1e6

    def fit(self, X, y=None):
        C, _ = self._split(X)
        L, cpm = self._logcpm(C)
        expresado = (cpm >= self.min_cpm).mean(0) >= self.min_frac
        var = np.where(expresado, L.var(0), -1)
        n = min(self.n_genes, int(expresado.sum()))
        self.genes_ = np.sort(np.argsort(var)[::-1][:n])
        self.n_expresados_ = int(expresado.sum())
        return self

    def transform(self, X):
        C, T = self._split(X)
        L, _ = self._logcpm(C)
        out = L[:, self.genes_]
        return np.hstack([out, T]) if T is not None else out


def build_type_data(pb_ig, pb_orig, cfg):
    """Por tipo celular: (recuentos donantes x genes, técnicas donantes x 2). También la etiqueta de cada donante.

    Las técnicas (log10 de n_celulas y de total_counts ORIGINALES, antes de igualar) salen de pb_orig: son
    la información que el modelo trampa usa y la que queremos comprobar si la expresión supera.
    """
    C = cfg["columns"]
    don, lab, ct = C["donor"], C["label"], C["cell_type"]
    datos = {}
    for t in sorted(pb_ig.obs[ct].unique()):
        sub = pb_ig[pb_ig.obs[ct] == t]
        X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
        cuentas = pd.DataFrame(X, index=sub.obs[don].astype(str).values, columns=sub.var_names)
        orig = pb_orig.obs[pb_orig.obs[ct] == t].set_index(don)
        orig.index = orig.index.astype(str)
        tec = np.log10(orig.loc[cuentas.index, ["n_celulas", "total_counts"]].clip(lower=1).astype(float))
        tec.columns = ["log_n_celulas", "log_total_counts_original"]
        datos[t] = (cuentas, tec)
    y = pb_ig.obs.groupby(pb_ig.obs[don].astype(str), observed=True)[lab].first().astype(str)
    return datos, y


def run_baselines(datos, y, splits, modelos, orden, verbose=True):
    """Validación cruzada por donante de cada (tipo, modelo).

    modelos: {nombre: (fabrica, conjunto)} donde fabrica(n_tech) devuelve un Pipeline nuevo y conjunto es
             "expresion", "tecnicas" o "expresion+tecnicas".
    splits: lista (repeticion, train_idx, test_idx) sobre el índice de `y` (todos los donantes).

    Devuelve (métricas por repetición, tipo y modelo; probabilidades fuera de pliegue
              {(modelo, tipo): array repeticiones x donantes x clases, NaN si el donante no tiene el tipo}).
    """
    donantes = y.index.to_numpy()
    n_rep = max(r for r, _, _ in splits) + 1
    filas, probas = [], {}
    t0 = time.time()
    for t, (cuentas, tec) in datos.items():
        presente = np.isin(donantes, cuentas.index)
        for nombre, (fabrica, conjunto) in modelos.items():
            M, n_tech = _feature_matrix(cuentas, tec, conjunto)
            M = M.reindex(donantes)                       # NaN en donantes sin este tipo
            P = np.full((n_rep, len(donantes), len(orden)), np.nan)
            for rep, tr, te in splits:
                tr, te = tr[presente[tr]], te[presente[te]]
                modelo = fabrica(n_tech)
                modelo.fit(M.values[tr], y.values[tr])
                p = modelo.predict_proba(M.values[te])
                P[rep][np.ix_(te, [orden.index(c) for c in modelo.classes_])] = p
            probas[(nombre, t)] = P
            for rep in range(n_rep):
                ok = presente
                pred = np.array(orden)[np.nanargmax(P[rep][ok], axis=1)]
                filas.append({"repeticion": rep, "tipo": t, "modelo": nombre, "n_donantes": int(ok.sum()),
                              **_metrics(y.values[ok], pred, orden)})
            if verbose:
                print(f"  {t:10s} {nombre:28s} {time.time() - t0:6.0f} s")
    return pd.DataFrame(filas), probas


def fuse(probas, y, orden, modelo, tipos=None):
    """Fusión tardía: media de las probabilidades de los tipos disponibles para cada donante."""
    claves = [(m, t) for (m, t) in probas if m == modelo and (tipos is None or t in tipos)]
    P = np.nanmean(np.stack([probas[k] for k in claves]), axis=0)     # repeticiones x donantes x clases
    filas = []
    for rep in range(P.shape[0]):
        pred = np.array(orden)[np.argmax(P[rep], axis=1)]
        filas.append({"repeticion": rep, "tipo": "fusión", "modelo": modelo, "n_donantes": len(y),
                      **_metrics(y.values, pred, orden)})
    return pd.DataFrame(filas), P


def _feature_matrix(cuentas, tec, conjunto):
    if conjunto == "tecnicas":
        return tec, 0
    if conjunto == "expresion":
        return cuentas, 0
    return pd.concat([cuentas, tec], axis=1), tec.shape[1]


def restrict(datos, y, clases):
    """Se queda solo con los donantes de `clases` (p. ej. Low y High)."""
    y2 = y[y.isin(clases)]
    datos2 = {t: (c.loc[c.index.isin(y2.index)], tec.loc[tec.index.isin(y2.index)]) for t, (c, tec) in datos.items()}
    return datos2, y2


def auc_by_rep(P, y, positiva, orden):
    """AUC ROC por repetición para un problema binario, con los donantes que tienen predicción."""
    from sklearn.metrics import roc_auc_score

    k = orden.index(positiva)
    out = []
    for rep in range(P.shape[0]):
        ok = ~np.isnan(P[rep][:, k])
        out.append(roc_auc_score((y.values[ok] == positiva).astype(int), P[rep][ok, k]))
    return np.array(out)


def run_regression(datos, objetivo, splits, modelos, verbose=True):
    """Como run_baselines pero para un objetivo continuo (p. ej. el CPS).

    Métricas por repetición: correlación de Spearman y R² entre predicción fuera de pliegue y valor real.
    Devuelve (métricas, predicciones {(modelo, tipo): array repeticiones x donantes, NaN si falta el tipo}).
    """
    from scipy.stats import spearmanr
    from sklearn.metrics import r2_score

    donantes = objetivo.index.to_numpy()
    n_rep = max(r for r, _, _ in splits) + 1
    filas, preds = [], {}
    t0 = time.time()
    for t, (cuentas, tec) in datos.items():
        presente = np.isin(donantes, cuentas.index)
        for nombre, (fabrica, conjunto) in modelos.items():
            M, n_tech = _feature_matrix(cuentas, tec, conjunto)
            M = M.reindex(donantes)
            Y = np.full((n_rep, len(donantes)), np.nan)
            for rep, tr, te in splits:
                tr, te = tr[presente[tr]], te[presente[te]]
                modelo = fabrica(n_tech)
                modelo.fit(M.values[tr], objetivo.values[tr])
                Y[rep, te] = modelo.predict(M.values[te])
            preds[(nombre, t)] = Y
            for rep in range(n_rep):
                filas.append({"repeticion": rep, "tipo": t, "modelo": nombre, "n_donantes": int(presente.sum()),
                              "spearman": spearmanr(objetivo.values[presente], Y[rep, presente])[0],
                              "r2": r2_score(objetivo.values[presente], Y[rep, presente])})
            if verbose:
                print(f"  {t:10s} {nombre:28s} {time.time() - t0:6.0f} s")
    return pd.DataFrame(filas), preds


def fuse_regression(preds, objetivo, modelo):
    """Media de las predicciones de los tipos disponibles para cada donante."""
    from scipy.stats import spearmanr
    from sklearn.metrics import r2_score

    Y = np.nanmean(np.stack([v for (m, _), v in preds.items() if m == modelo]), axis=0)
    filas = [{"repeticion": rep, "tipo": "fusión", "modelo": modelo, "n_donantes": len(objetivo),
              "spearman": spearmanr(objetivo.values, Y[rep])[0], "r2": r2_score(objetivo.values, Y[rep])}
             for rep in range(Y.shape[0])]
    return pd.DataFrame(filas), Y


def summarize_by(df, keys=("tipo", "modelo"), metricas=("balanced_accuracy", "f1_macro", "kappa_cuadratico")):
    """Media y desviación de las métricas agrupando por `keys`."""
    metricas = list(metricas)
    g = df.groupby(list(keys))
    out = g[metricas].mean().add_suffix("_media").join(g[metricas].std().add_suffix("_sd"))
    out["n_donantes"] = g["n_donantes"].first()
    return out[[c for m in metricas for c in (f"{m}_media", f"{m}_sd")] + ["n_donantes"]]
