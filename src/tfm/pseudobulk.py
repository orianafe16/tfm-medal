"""
Pseudobulk: de células a (donante, tipo celular).

    build_pseudobulk(path, ...)       suma los recuentos de raw por (donante, tipo) leyendo el .h5ad a trozos
    log_cpm(pb)                       normaliza a recuentos por millón y aplica log1p
    expressed_genes(pb, ...)          genes expresados en cada tipo celular
    technical_by_class(pb, cfg)       ¿cambian n_celulas o los recuentos totales con la clase?
    gene_symbol_index(pb)             símbolo -> posición, para buscar genes por nombre

Por qué sumar recuentos: la suma de recuentos de varias células sigue siendo un recuento, y se puede
normalizar como un experimento de RNA-seq bulk. Promediar valores log-normalizados mezclaría escalas.
"""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from .data import _read_rows, clean_gene_name, get_donor_table, read_obs_var


def build_pseudobulk(path, donor_col, group_col, source="raw/X", chunk_cells=4000,
                     donor_expected=None, verbose=True):
    """Devuelve un AnnData con una fila por (donante, tipo celular) y los recuentos sumados.

    obs: donante, tipo, n_celulas, total_counts, n_genes_detectados y los metadatos constantes del donante.
    var: genes de `source` (IDs de ENSEMBL), con feature_name, símbolo limpio y si están en X (HVG).
    """
    import anndata as ad

    obs, var, raw_var, _, _ = read_obs_var(path)
    genes = raw_var if source.startswith("raw") else var

    claves = obs[donor_col].astype(str) + "|" + obs[group_col].astype(str)
    codigos, grupos = pd.factorize(claves, sort=True)
    n_grupos = len(grupos)

    acumulado = None
    with h5py.File(path, "r") as f:
        n_cel = obs.shape[0]
        for a in range(0, n_cel, chunk_cells):
            b = min(a + chunk_cells, n_cel)
            M = _read_rows(f, source, a, b)                       # células a..b x genes
            G = sp.csr_matrix((np.ones(b - a), (codigos[a:b], np.arange(b - a))), shape=(n_grupos, b - a))
            parte = (G @ M.astype(np.float64)).tocsr()           # suma por grupo
            acumulado = parte if acumulado is None else acumulado + parte
            if verbose:
                print(f"\r  células {b:,}/{n_cel:,}", end="")
    if verbose:
        print()

    partes = pd.Series(grupos).str.split("|", n=1, expand=True)
    pb_obs = pd.DataFrame({donor_col: partes[0].values, group_col: partes[1].values},
                          index=pd.Index(grupos, name="pseudobulk"))
    pb_obs["n_celulas"] = np.bincount(codigos, minlength=n_grupos)
    pb_obs["total_counts"] = np.asarray(acumulado.sum(1)).ravel()
    pb_obs["n_genes_detectados"] = np.diff(acumulado.indptr)

    donantes = get_donor_table(obs, donor_col, expected=donor_expected)
    meta = donantes.drop(columns="n_celulas")
    meta.index = meta.index.astype(str)
    pb_obs = pb_obs.join(meta, on=donor_col)
    for c in pb_obs.columns:                     # h5ad no admite columnas object mezcladas
        if pb_obs[c].dtype == object:
            pb_obs[c] = pb_obs[c].astype(str)

    pb_var = genes.copy()
    pb_var.index = pb_var.index.astype(str)
    nombre = pb_var["feature_name"].astype(str) if "feature_name" in pb_var else pb_var.index.to_series()
    pb_var["simbolo"] = nombre.map(clean_gene_name).values
    if source.startswith("raw") and "feature_name" in var and "feature_name" in raw_var:
        pb_var["en_X_hvg"] = raw_var["feature_name"].astype(str).isin(set(var["feature_name"].astype(str))).values
    for c in pb_var.columns:
        if pb_var[c].dtype == object:
            pb_var[c] = pb_var[c].astype(str)

    pb = ad.AnnData(X=acumulado.astype(np.float32), obs=pb_obs, var=pb_var)
    pb.uns["pseudobulk"] = {"fuente": source, "agrupado_por": [donor_col, group_col], "operacion": "suma de recuentos"}
    return pb


def log_cpm(pb, layer=None):
    """log1p(CPM): cada pseudobulk se escala a 1 millón de recuentos y se aplica log(1 + x). Devuelve denso."""
    X = pb.X if layer is None else pb.layers[layer]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    tot = X.sum(1, keepdims=True)
    tot[tot == 0] = 1
    return np.log1p(X / tot * 1e6)


def expressed_genes(pb, group_col, min_cpm=1.0, min_frac=0.5, min_cells=0):
    """Por tipo celular: genes con CPM >= min_cpm en al menos min_frac de sus pseudobulk.

    Devuelve un DataFrame booleano genes x tipos. Es la idea de filterByExpr (edgeR): quitar genes que
    apenas se detectan, que solo añaden ruido. No usa la etiqueta, así que no hay fuga de información.
    """
    X = pb.X.toarray() if sp.issparse(pb.X) else np.asarray(pb.X)
    tot = X.sum(1, keepdims=True)
    tot[tot == 0] = 1
    cpm = X / tot * 1e6
    res = {}
    for t in sorted(pb.obs[group_col].unique()):
        m = ((pb.obs[group_col] == t) & (pb.obs["n_celulas"] >= min_cells)).to_numpy()
        res[t] = (cpm[m] >= min_cpm).mean(0) >= min_frac if m.any() else np.zeros(X.shape[1], bool)
    return pd.DataFrame(res, index=pb.var_names)


def technical_by_class(pb, cfg, exclude=True):
    """Para cada tipo celular: medianas de n_celulas y total_counts por clase y correlación de Spearman
    entre la clase ordinal (Low < Intermediate < High) y esas variables. Contado por pseudobulk."""
    C = cfg["columns"]
    lab, ct = C["label"], C["cell_type"]
    orden = cfg["labels"]["classification"]
    obs = pb.obs[pb.obs[lab].isin(orden)] if exclude else pb.obs
    rango = obs[lab].map({c: i for i, c in enumerate(orden)})
    filas = []
    for t, g in obs.groupby(ct, observed=True):
        r = rango.loc[g.index]
        fila = {ct: t}
        for v in ["n_celulas", "total_counts"]:
            for c in orden:
                fila[f"{v}_mediana_{c}"] = float(g.loc[g[lab] == c, v].median())
            rho, p = stats.spearmanr(r, g[v])
            fila[f"{v}_rho"], fila[f"{v}_p"] = float(rho), float(p)
        filas.append(fila)
    return pd.DataFrame(filas).set_index(ct)


def pca(M, n_genes=2000, n_comps=10, seed=0):
    """PCA sobre una matriz muestras x genes (log-CPM): se queda con los n_genes de mayor varianza,
    estandariza cada gen y proyecta. Devuelve (componentes, varianza_explicada, índices_de_genes)."""
    from sklearn.decomposition import PCA

    var = M.var(0)
    idx = np.argsort(var)[::-1][:min(n_genes, int((var > 0).sum()))]
    Z = M[:, idx]
    Z = (Z - Z.mean(0)) / Z.std(0)
    n_comps = min(n_comps, Z.shape[0] - 1, Z.shape[1])
    modelo = PCA(n_components=n_comps, random_state=seed).fit(Z)
    return modelo.transform(Z), modelo.explained_variance_ratio_, idx


def r2_by(pcs, grupos):
    """R² de cada componente explicado por una variable: categórica (ANOVA) o numérica (correlación²)."""
    g = pd.Series(np.asarray(grupos))
    out = []
    for i in range(pcs.shape[1]):
        y = pd.Series(pcs[:, i])
        if pd.api.types.is_numeric_dtype(g):
            out.append(float(np.corrcoef(y, g.astype(float))[0, 1] ** 2))
        else:
            medias = y.groupby(g.values).transform("mean")
            out.append(float(1 - ((y - medias) ** 2).sum() / ((y - y.mean()) ** 2).sum()))
    return np.array(out)


def gene_symbol_index(pb):
    """Diccionario símbolo -> posición del primer gen con ese símbolo."""
    s = pd.Series(np.arange(pb.n_vars), index=pb.var["simbolo"].astype(str).values)
    return s[~s.index.duplicated()].to_dict()


# --------------------------------------------------------------------------- igualar profundidad
def _eligible(pb, cfg, min_cells):
    C = cfg["columns"]
    return (pb.obs[C["label"]].isin(cfg["labels"]["classification"]) & (pb.obs["n_celulas"] >= min_cells)).to_numpy()


def depth_targets(pb, cfg, quantile, min_cells):
    """Profundidad objetivo de cada tipo celular: el cuantil `quantile` de total_counts entre los
    pseudobulk que se usarán para clasificar (sin Reference y con al menos min_cells células)."""
    ct = cfg["columns"]["cell_type"]
    ok = _eligible(pb, cfg, min_cells)
    return pb.obs[ok].groupby(ct, observed=True)["total_counts"].quantile(quantile)


def thinning_losses(pb, cfg, quantiles, min_cells):
    """Cuántos pseudobulk de cada clase se perderían con cada cuantil (por profundidad o por min_cells)."""
    C = cfg["columns"]
    lab, ct = C["label"], C["cell_type"]
    orden = cfg["labels"]["classification"]
    obs = pb.obs[pb.obs[lab].isin(orden)]
    filas = []
    for q in quantiles:
        objetivo = depth_targets(pb, cfg, q, min_cells)
        queda = (obs["n_celulas"] >= min_cells) & (obs["total_counts"] >= obs[ct].map(objetivo).astype(float))
        fila = {"cuantil": q, "profundidad_mediana_objetivo": float(objetivo.median())}
        for c in orden:
            m = obs[lab] == c
            fila[f"pierde_{c}"] = int((m & ~queda).sum())
            fila[f"queda_{c}"] = int((m & queda).sum())
        filas.append(fila)
    return pd.DataFrame(filas).set_index("cuantil")


def downsample_to_depth(pb, cfg, quantile=0.05, min_cells=10, seed=0):
    """Iguala la profundidad dentro de cada tipo celular por submuestreo binomial.

    Para cada pseudobulk con total T >= objetivo D de su tipo, cada recuento x pasa a Binomial(x, D/T):
    es como si se hubieran secuenciado solo D lecturas. El total queda en ~D y las proporciones entre
    genes se conservan en esperanza. Se descartan Reference, los pseudobulk con < min_cells células
    y los que no llegan a D. Devuelve (pb_igualado, informe por tipo y clase).
    """
    import anndata as ad

    C = cfg["columns"]
    lab, ct = C["label"], C["cell_type"]
    objetivo = depth_targets(pb, cfg, quantile, min_cells)
    tot = pb.obs["total_counts"].to_numpy()
    D = pb.obs[ct].map(objetivo).astype(float).to_numpy()
    queda = _eligible(pb, cfg, min_cells) & (tot >= D)

    sub = pb[queda].copy()
    X = sub.X.tocsr().astype(np.float64) if sp.issparse(sub.X) else sp.csr_matrix(sub.X, dtype=np.float64)
    p = np.repeat(D[queda] / tot[queda], np.diff(X.indptr))
    rng = np.random.default_rng(seed)
    X.data = rng.binomial(X.data.astype(np.int64), np.clip(p, 0, 1)).astype(np.float64)
    X.eliminate_zeros()

    sub = ad.AnnData(X=X.astype(np.float32), obs=sub.obs.copy(), var=sub.var.copy(), uns=dict(sub.uns))
    sub.obs["total_counts_original"] = sub.obs["total_counts"]
    sub.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
    sub.obs["n_genes_detectados"] = np.diff(X.indptr)
    sub.uns["igualado"] = {"cuantil": quantile, "min_cells": min_cells, "semilla": seed,
                           "objetivo_por_tipo": {k: float(v) for k, v in objetivo.items()}}

    orden = cfg["labels"]["classification"]
    obs = pb.obs
    informe = pd.DataFrame({
        "objetivo": objetivo,
        **{f"queda_{c}": obs[queda & (obs[lab] == c).to_numpy()].groupby(ct, observed=True).size() for c in orden},
        **{f"pierde_{c}": obs[~queda & (obs[lab] == c).to_numpy()].groupby(ct, observed=True).size() for c in orden},
    }).fillna(0)
    for c in informe.columns[1:]:
        informe[c] = informe[c].astype(int)
    return sub, informe


def pca_r2_by_type(pb, cfg, min_cpm=1.0, min_frac=0.5, n_genes=2000, n_comps=5, seed=0):
    """PCA dentro de cada tipo celular (sin Reference) y R² de cada PC con ADNC, log10(n_celulas) y
    log10(total_counts). Devuelve (tabla larga de R², {tipo: (pcs, var_explicada, obs)})."""
    C = cfg["columns"]
    lab, ct = C["label"], C["cell_type"]
    orden = cfg["labels"]["classification"]
    L = log_cpm(pb)
    expr = expressed_genes(pb, ct, min_cpm=min_cpm, min_frac=min_frac)
    filas, salida = [], {}
    for t in sorted(pb.obs[ct].unique()):
        m = ((pb.obs[ct] == t) & pb.obs[lab].isin(orden)).to_numpy()
        if m.sum() < 5:
            continue
        sub = pb.obs[m]
        p, v, _ = pca(L[np.ix_(m, expr[t].to_numpy())], n_genes=n_genes, n_comps=n_comps, seed=seed)
        salida[t] = (p, v, sub)
        variables = {"ADNC": sub[lab].astype(str),
                     "log n_celulas": np.log10(sub["n_celulas"].clip(lower=1)),
                     "log total_counts": np.log10(sub["total_counts"].clip(lower=1))}
        for nombre, valores in variables.items():
            for i, r2 in enumerate(r2_by(p, valores)):
                filas.append({"tipo": t, "variable": nombre, "PC": f"PC{i + 1}", "R2": r2, "var_explicada": float(v[i])})
    return pd.DataFrame(filas), salida
