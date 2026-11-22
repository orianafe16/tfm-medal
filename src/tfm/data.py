"""
Carga y descripción del subset del SEA-AD.

Funciones principales:
    load_adata(path, backed)        carga el .h5ad, en modo backed si no cabe con holgura en RAM
    describe_matrices(path)         X, raw y layers: forma, formato, RAM y si son recuentos o normalizados
    check_normalization(path)       comprueba que X = log1p(raw / total * T) y estima T
    classify_obs_columns(obs, ...)  agrupa las columnas de obs según configs/config.yaml
    get_donor_table(adata, ...)     una fila por donante con los metadatos constantes dentro del donante
    check_gene_names(adata)         duplicados, IDs de ENSEMBL, genes mitocondriales, X ⊂ raw
    check_expectations(obs, ...)    compara con lo que se espera del subset (config: expected)

Casi todo lee el fichero con h5py a trozos, así que funciona sin cargar la matriz en memoria.
"""
import re
import warnings

import h5py
import numpy as np
import pandas as pd

try:  # anndata >= 0.11
    from anndata.io import read_elem
except ImportError:  # anndata 0.9-0.10
    from anndata.experimental import read_elem

GB = 1024**3
ENSEMBL = re.compile(r"^ENS[A-Z]*G\d{6,}")


# --------------------------------------------------------------------------- utilidades h5py
def _matrix_paths(f):
    rutas = ["X"]
    if "raw" in f and "X" in f["raw"]:
        rutas.append("raw/X")
    if "layers" in f:
        rutas += [f"layers/{k}" for k in f["layers"].keys()]
    return rutas


def _matrix_info(f, ruta):
    """Forma, formato, tipo de dato, no-ceros y GB en RAM de una matriz, sin leerla."""
    obj = f[ruta]
    if isinstance(obj, h5py.Dataset):
        forma = tuple(int(v) for v in obj.shape)
        n = int(np.prod(forma))
        return dict(matriz=ruta, formato="densa", n_filas=forma[0], n_cols=forma[1],
                    dtype=str(obj.dtype), nnz=n, gb_ram=n * obj.dtype.itemsize / GB)
    forma = tuple(int(v) for v in obj.attrs["shape"])
    nbytes = sum(obj[k].size * obj[k].dtype.itemsize for k in ("data", "indices", "indptr"))
    return dict(matriz=ruta, formato=str(obj.attrs.get("encoding-type", "dispersa")),
                n_filas=forma[0], n_cols=forma[1], dtype=str(obj["data"].dtype),
                nnz=int(obj["data"].size), gb_ram=nbytes / GB)


def _sample_rows(n_rows, n_blocks=20, block=50, seed=0):
    """Índices de filas en bloques contiguos repartidos por el fichero (leer bloques es rápido)."""
    rng = np.random.default_rng(seed)
    starts = np.unique(rng.integers(0, max(1, n_rows - block), n_blocks))
    return [(int(s), int(min(s + block, n_rows))) for s in starts]


def _read_rows(f, ruta, start, stop):
    """Devuelve las filas [start, stop) como scipy.sparse.csr_matrix."""
    import scipy.sparse as sp

    obj = f[ruta]
    if isinstance(obj, h5py.Dataset):
        return sp.csr_matrix(obj[start:stop])
    enc = obj.attrs.get("encoding-type", "")
    n_cols = int(obj.attrs["shape"][1])
    if enc != "csr_matrix":
        raise NotImplementedError(f"{ruta} está en formato {enc}; solo se leen por filas matrices CSR o densas")
    ptr = obj["indptr"][start:stop + 1]
    data = obj["data"][ptr[0]:ptr[-1]]
    ind = obj["indices"][ptr[0]:ptr[-1]]
    return sp.csr_matrix((data, ind, ptr - ptr[0]), shape=(stop - start, n_cols))


def read_obs_var(path):
    """Lee solo obs, var y raw.var (tablas pequeñas), sin tocar las matrices."""
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        raw_var = read_elem(f["raw"]["var"]) if "raw" in f and "var" in f["raw"] else None
        obsm = {k: tuple(int(v) for v in (f["obsm"][k].shape if isinstance(f["obsm"][k], h5py.Dataset)
                                         else f["obsm"][k].attrs["shape"]))
                for k in f["obsm"].keys()} if "obsm" in f else {}
        uns = list(f["uns"].keys()) if "uns" in f else []
    return obs, var, raw_var, obsm, uns


# --------------------------------------------------------------------------- carga
def estimate_ram_gb(path):
    """GB que ocuparían en RAM todas las matrices del fichero."""
    with h5py.File(path, "r") as f:
        return sum(_matrix_info(f, r)["gb_ram"] for r in _matrix_paths(f))


def load_adata(path, backed="auto", margin=2.0):
    """Carga el .h5ad.

    backed:
        "auto" -> backed solo si las matrices x `margin` no caben en la RAM disponible
        "r"    -> siempre backed (X y raw se quedan en disco; obs y var en memoria)
        None   -> carga todo en memoria
    """
    import anndata as ad

    if backed == "auto":
        necesario = estimate_ram_gb(path) * margin
        try:
            import psutil
            libre = psutil.virtual_memory().available / GB
        except ImportError:
            libre = float("inf")
        backed = "r" if necesario > libre else None
        print(f"Matrices x{margin}: {necesario:.1f} GB | RAM libre: {libre:.1f} GB -> "
              f"{'modo backed' if backed else 'carga completa'}")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        return ad.read_h5ad(path, backed=backed)


# --------------------------------------------------------------------------- matrices
def _diagnose(vals, row_sums, expm1_sums):
    vals = vals[vals != 0]
    if vals.size == 0:
        return "vacía"
    if (vals < 0).any():
        return "escalada (tiene negativos): no son recuentos ni log-normalizado"
    frac_int = np.mean(np.isclose(vals, np.round(vals)))
    if frac_int > 0.999:
        return "recuentos crudos (enteros)"
    cv = lambda v: float(np.std(v) / np.mean(v)) if np.mean(v) > 0 else np.nan
    if cv(row_sums) < 0.02:
        return f"normalizado SIN log (suma por célula constante ~{np.median(row_sums):.0f})"
    if vals.max() < 20:
        if cv(expm1_sums) < 0.02:
            return f"log-normalizado: log1p(cuentas normalizadas a ~{np.median(expm1_sums):.0f} por célula)"
        return "log-normalizado (decimales, máx < 20); la suma no es constante porque no están todos los genes"
    return "decimales sin patrón estándar: revisar"


def describe_matrices(path):
    """Una fila por matriz (X, raw/X, layers/*) con forma, formato, RAM y diagnóstico de una muestra."""
    filas = []
    with h5py.File(path, "r") as f:
        for ruta in _matrix_paths(f):
            info = _matrix_info(f, ruta)
            info["densidad"] = info["nnz"] / (info["n_filas"] * info["n_cols"])
            bloques = [_read_rows(f, ruta, a, b) for a, b in _sample_rows(info["n_filas"])]
            import scipy.sparse as sp
            m = sp.vstack(bloques).tocsr()
            vals = m.data
            me = m.copy()
            with np.errstate(over="ignore"):  # en recuentos, expm1 desborda: es esperable y no se usa
                me.data = np.expm1(me.data)
            filas.append({**info,
                          "min_no_cero": float(vals[vals != 0].min()) if vals.size else np.nan,
                          "max": float(vals.max()) if vals.size else np.nan,
                          "frac_enteros": float(np.mean(np.isclose(vals, np.round(vals)))) if vals.size else np.nan,
                          "suma_celula_mediana": float(np.median(np.asarray(m.sum(1)).ravel())),
                          "diagnostico": _diagnose(vals, np.asarray(m.sum(1)).ravel(),
                                                   np.asarray(me.sum(1)).ravel()),
                          "celulas_muestreadas": m.shape[0]})
    return pd.DataFrame(filas).set_index("matriz")


def clean_gene_name(s):
    """'SIMBOLO_ENSG...' -> 'SIMBOLO'; 'ENSG...12.3' -> 'ENSG...12' (quita la versión)."""
    s = str(s)
    if "_" in s and not s.startswith("ENS"):
        s = s.split("_")[0]
    return re.sub(r"^(ENS[A-Z]*G\d+)\.\d+$", r"\1", s)


def map_genes_x_to_raw(var, raw_var):
    """Posición en raw de cada gen de X (-1 si no se encuentra).

    Prueba varias formas de emparejar (X y raw pueden usar símbolos o IDs de ENSEMBL) y se queda con la
    que más genes encuentra. Devuelve (posiciones, nombre_de_la_estrategia).
    """
    def lookup(claves_x, claves_raw):
        pos = pd.Series(np.arange(len(claves_raw)), index=pd.Index(claves_raw))
        pos = pos[~pos.index.duplicated()]
        return pos.reindex(pd.Index(claves_x)).fillna(-1).astype(int).to_numpy()

    ix, rx = var.index.astype(str), raw_var.index.astype(str)
    # En caso de empate gana la primera: de más exacta a menos
    candidatas = {"índice X = índice raw": (ix, rx)}
    if "feature_name" in raw_var:
        fr = raw_var["feature_name"].astype(str)
        if "feature_name" in var:
            candidatas["X.feature_name = raw.feature_name"] = (var["feature_name"].astype(str), fr)
        candidatas["índice X (limpio) = raw.feature_name (limpio)"] = (ix.map(clean_gene_name), fr.map(clean_gene_name))
    candidatas["índice X (limpio) = índice raw (limpio)"] = (ix.map(clean_gene_name), rx.map(clean_gene_name))
    mejor = max(candidatas, key=lambda k: (lookup(*candidatas[k]) >= 0).sum())
    return lookup(*candidatas[mejor]), mejor


def gene_tables(path):
    """Tablas de genes de X y de raw con sus columnas, para inspeccionarlas."""
    _, var, raw_var, _, _ = read_obs_var(path)
    var = var.copy()
    var.index.name = "gen_X"
    if raw_var is not None:
        pos, estrategia = map_genes_x_to_raw(var, raw_var)
        var["id_en_raw"] = np.where(pos >= 0, raw_var.index.to_numpy()[np.maximum(pos, 0)], None)
        raw_var = raw_var.copy()
        raw_var.index.name = "gen_raw"
        raw_var["en_X"] = raw_var.index.isin(var["id_en_raw"].dropna())
        var.attrs["estrategia"] = estrategia
    return var, raw_var


def check_normalization(path, x="X", raw="raw/X"):
    """Comprueba si X se obtuvo de raw como log1p(raw / total_celula * T) y estima T.

    Para cada entrada no nula: T = expm1(X) * total_raw_celula / raw. Si X viene de raw con esa
    normalización, T sale igual en todas las entradas y células (p. ej. 10.000).
    """
    obs, var, raw_var, _, _ = read_obs_var(path)
    if raw_var is None:
        return {"conclusion": "no hay raw: no se puede comprobar"}
    cols_raw, estrategia = map_genes_x_to_raw(var, raw_var)
    validos = cols_raw >= 0
    cols_raw = cols_raw[validos]
    if validos.sum() == 0:
        return {"genes_X_encontrados_en_raw": 0, "conclusion": "ningún gen de X se encuentra en raw: revisar nombres"}

    T_celula, coinciden = [], []
    with h5py.File(path, "r") as f:
        for a, b in _sample_rows(obs.shape[0], n_blocks=10, block=50, seed=1):
            mx = _read_rows(f, x, a, b).toarray()[:, validos]
            mr_full = _read_rows(f, raw, a, b)
            total = np.asarray(mr_full.sum(1)).ravel()
            mr = mr_full[:, cols_raw].toarray()
            for i in range(mx.shape[0]):
                nz = (mx[i] > 0) & (mr[i] > 0)
                coinciden.append(np.mean((mx[i] > 0) == (mr[i] > 0)))
                if nz.sum() >= 5 and total[i] > 0:
                    T = np.expm1(mx[i, nz]) * total[i] / mr[i, nz]
                    T_celula.append((np.median(T), np.std(T) / np.mean(T)))
    T_celula = np.array(T_celula)
    res = {
        "emparejamiento_genes": estrategia,
        "genes_X_encontrados_en_raw": int(validos.sum()),
        "genes_X": int(len(var)),
        "frac_ceros_coinciden_X_raw": float(np.mean(coinciden)),
        "celulas_comprobadas": int(len(T_celula)),
        "T_mediana": float(np.median(T_celula[:, 0])) if len(T_celula) else np.nan,
        "T_cv_entre_celulas": float(np.std(T_celula[:, 0]) / np.mean(T_celula[:, 0])) if len(T_celula) else np.nan,
        "T_cv_dentro_celula_mediana": float(np.median(T_celula[:, 1])) if len(T_celula) else np.nan,
    }
    ok = (res["frac_ceros_coinciden_X_raw"] > 0.99 and res["T_cv_entre_celulas"] < 0.02
          and res["T_cv_dentro_celula_mediana"] < 0.02)
    res["conclusion"] = (f"X = log1p(raw / total_celula x {res['T_mediana']:.0f}), con el total sobre todos los genes de raw"
                         if ok else "X no encaja con log1p(raw normalizado por célula): revisar")
    return res


# --------------------------------------------------------------------------- obs
def _donor_level(obs, donor_col):
    """True para las columnas con un único valor dentro de cada donante."""
    n = obs.groupby(donor_col, observed=True).nunique(dropna=False)
    return (n.max() <= 1).reindex(obs.columns, fill_value=True)


def classify_obs_columns(obs, groups, donor_col):
    """Tabla con una fila por columna de obs: grupo (según config), tipo, nº de valores, nivel y ejemplos."""
    grupo_de = {c: g for g, cols in groups.items() for c in cols}
    nivel = _donor_level(obs, donor_col)
    filas = []
    for c in obs.columns:
        s = obs[c]
        ej = [str(v) for v in pd.unique(s.dropna())[:3]]
        filas.append({"columna": c, "grupo": grupo_de.get(c, "sin_clasificar"), "dtype": str(s.dtype),
                      "n_valores": int(s.nunique()), "n_nulos": int(s.isna().sum()),
                      "nivel": "donante" if nivel[c] else "célula", "ejemplos": "; ".join(ej)})
    t = pd.DataFrame(filas)
    orden = list(groups) + ["sin_clasificar"]
    t["grupo"] = pd.Categorical(t["grupo"], categories=orden, ordered=True)
    sobran = [c for c in grupo_de if c not in obs.columns]
    if sobran:
        print(f"Columnas del config que no están en obs: {sobran}")
    if (t["grupo"] == "sin_clasificar").any():
        warnings.warn(f"Columnas sin clasificar: {t.loc[t.grupo == 'sin_clasificar', 'columna'].tolist()}")
    return t.sort_values(["grupo", "columna"]).reset_index(drop=True)


def get_donor_table(adata, donor_col="donor_id", expected=None):
    """Una fila por donante con las columnas de obs que son constantes dentro de cada donante.

    adata: AnnData (en memoria o backed) o directamente el DataFrame obs.
    expected: columnas que deberían ser constantes por donante. Si alguna no lo es, avisa
    y la deja fuera. La lista de columnas que varían queda en tabla.attrs["no_constantes"].
    """
    obs = adata if isinstance(adata, pd.DataFrame) else adata.obs
    n = obs.groupby(donor_col, observed=True).nunique(dropna=False)
    constantes = [c for c in obs.columns if c != donor_col and n[c].max() <= 1]
    tabla = obs.groupby(donor_col, observed=True)[constantes].first()
    tabla.insert(0, "n_celulas", obs.groupby(donor_col, observed=True).size())

    problemas = []
    for c in expected or []:
        if c not in obs.columns:
            problemas.append({"columna": c, "problema": "no está en obs", "n_donantes": np.nan, "ejemplo": ""})
        elif n[c].max() > 1:
            malos = n.index[n[c] > 1]
            ej = obs.loc[obs[donor_col] == malos[0], c].astype(str).value_counts().to_dict()
            problemas.append({"columna": c, "problema": "varía dentro del donante",
                              "n_donantes": int(len(malos)), "ejemplo": f"{malos[0]}: {ej}"})
    tabla.attrs["no_constantes"] = problemas
    if problemas:
        warnings.warn("Columnas que deberían ser de donante y no lo son: "
                      + ", ".join(f"{p['columna']} ({p['problema']})" for p in problemas))
    return tabla


# --------------------------------------------------------------------------- genes
def _gene_summary(names, nombre):
    names = pd.Index(names.astype(str))
    dup = names[names.duplicated(keep=False)]
    return {
        "matriz": nombre,
        "n_genes": len(names),
        "n_nombres_repetidos": int(dup.nunique()),
        "n_filas_con_nombre_repetido": int(len(dup)),
        "ejemplos_repetidos": "; ".join(dup.unique()[:8]),
        "n_ensembl": int(names.str.match(ENSEMBL).sum()),
        "ejemplos_ensembl": "; ".join(names[names.str.match(ENSEMBL)][:5]),
        "n_mitocondriales_MT-": int(names.str.upper().str.startswith("MT-").sum()),
        "n_con_guion_bajo": int(names.str.contains("_").sum()),
    }


def check_gene_names(adata=None, var=None, raw_var=None):
    """Resumen de los nombres de gen de X y raw: duplicados, ENSEMBL sin símbolo, MT-, X ⊂ raw."""
    if adata is not None:
        var = adata.var
        raw_var = adata.raw.var if adata.raw is not None else None
    filas = [_gene_summary(var.index, "X (índice)")]
    if "feature_name" in var:
        filas.append(_gene_summary(var["feature_name"], "X (feature_name)"))
    extra = {"columnas_var": list(var.columns)}
    if raw_var is not None:
        filas.append(_gene_summary(raw_var.index, "raw (índice)"))
        if "feature_name" in raw_var:
            filas.append(_gene_summary(raw_var["feature_name"], "raw (feature_name)"))
        extra["columnas_raw_var"] = list(raw_var.columns)
        pos, estrategia = map_genes_x_to_raw(var, raw_var)
        extra["emparejamiento_X_raw"] = estrategia
        extra["genes_X_no_en_raw"] = int((pos < 0).sum())
        extra["ejemplos_X_no_en_raw"] = list(var.index[pos < 0][:10].astype(str))
    for col in ("feature_biotype", "feature_is_filtered", "highly_variable"):
        if col in var:
            extra[f"X_{col}"] = {str(k): int(v) for k, v in var[col].value_counts().items()}
    return pd.DataFrame(filas).set_index("matriz"), extra


# --------------------------------------------------------------------------- comprobaciones
def check_expectations(obs, cfg, n_genes_x=None, n_genes_raw=None):
    """Compara el fichero con lo esperado (config: expected). Devuelve una tabla esperado/observado."""
    C, E = cfg["columns"], cfg["expected"]
    lab, don, ct = C["label"], C["donor"], C["cell_type"]
    filas = [("Células totales", E["n_cells"], len(obs)),
             ("Donantes", E["n_donors"], obs[don].nunique()),
             ("Tipos celulares", E["n_cell_types"], obs[ct].nunique())]
    for clase, n in obs[lab].value_counts().reindex(cfg["labels"]["all"]).items():
        filas.append((f"Células {clase}", E["cells_per_class"], int(0 if pd.isna(n) else n)))
    if n_genes_x is not None:
        filas.append(("Genes en X", E["n_genes_X"], n_genes_x))
    if n_genes_raw is not None:
        filas.append(("Genes en raw", E["n_genes_raw"], n_genes_raw))
    etiquetas_por_donante = obs.groupby(don, observed=True)[lab].nunique()
    filas.append(("Donantes con más de un ADNC", 0, int((etiquetas_por_donante > 1).sum())))
    t = pd.DataFrame(filas, columns=["comprobacion", "esperado", "observado"])
    t["ok"] = (t["esperado"] == t["observado"]).map({True: "sí", False: "NO"})
    return t


def confounder_tables(donor_table, cfg, cols, exclude=True):
    """Cruza cada confusor con la clase, contando DONANTES (no células). Devuelve {columna: tabla}."""
    lab = cfg["columns"]["label"]
    dt = donor_table
    if exclude:
        dt = dt[~dt[lab].isin(cfg["labels"]["exclude"])]
    orden = [c for c in cfg["labels"]["all"] if c in set(dt[lab])]
    out = {}
    for c in cols:
        if c in dt.columns:
            out[c] = pd.crosstab(dt[lab].astype(str), dt[c].astype(str)).reindex(orden)
    return out


def donors_per_class(donor_table, cfg):
    lab = cfg["columns"]["label"]
    t = donor_table.groupby(lab, observed=True).agg(donantes=("n_celulas", "size"),
                                                   celulas=("n_celulas", "sum"),
                                                   celulas_por_donante_mediana=("n_celulas", "median"),
                                                   celulas_por_donante_min=("n_celulas", "min"),
                                                   celulas_por_donante_max=("n_celulas", "max"))
    return t.reindex([c for c in cfg["labels"]["all"] if c in t.index])
