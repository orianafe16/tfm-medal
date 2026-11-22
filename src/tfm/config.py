"""Lectura de configs/config.yaml y rutas del proyecto."""
from pathlib import Path

import yaml


def find_project_root(start=None):
    """Sube carpetas desde `start` hasta encontrar configs/config.yaml.

    Así el código funciona igual desde la raíz, desde notebooks/ o desde cualquier subcarpeta.
    """
    p = Path(start or Path.cwd()).resolve()
    for d in [p, *p.parents]:
        if (d / "configs" / "config.yaml").exists():
            return d
    raise FileNotFoundError("No encuentro configs/config.yaml en esta carpeta ni en las superiores")


def load_config(root=None):
    """Devuelve el config como diccionario, con las rutas convertidas a absolutas."""
    root = Path(root) if root else find_project_root()
    with open(root / "configs" / "config.yaml", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    cfg["root"] = root
    cfg["paths"] = {k: root / v for k, v in cfg["paths"].items()}
    return cfg
