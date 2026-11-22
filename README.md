# TFM · ADNC a partir de snRNA-seq del SEA-AD

Predecir el estadio de Alzheimer (ADNC: Low / Intermediate / High) de cada donante a partir de snRNA-seq del
SEA-AD (MTG): primero con modelos sencillos sobre pseudobulk y después con GNNs.
Grupo MEDAL (CTB-UPM), tutora Lucía Prieto.

## Estructura

```
Data/raw/          subset_ADNC_GNN.h5ad (no se modifica nunca)
Data/interim/      ficheros intermedios
Data/processed/    datos listos para modelar (p. ej. pseudobulk)
notebooks/         solo llaman a src/tfm y enseñan resultados
src/tfm/           lógica reutilizable (config.py, data.py)
configs/           config.yaml: rutas, semilla, columnas clave, clasificación de obs
results/figures/   figuras
results/tables/    tablas que generan los notebooks
docs/              datos.md (descripción del dataset), decisiones.md (registro de decisiones)
```

`Data/`, `results/` y `notes/` (notas de trabajo personales) no se suben a git.

## Cómo ejecutar un notebook

En VS Code el kernel puede caerse al importar h5py: arranca sin activar el entorno y Windows carga otra
DLL de HDF5. Desde una terminal con el entorno activado funciona:

```powershell
conda activate tfm-medal
cd C:\Users\ori\MEDAL\Code
jupyter nbconvert --to notebook --execute --inplace notebooks\00_carga_y_ficha.ipynb
```

Después se abre el notebook en VS Code para leer las salidas. Para trabajar de forma interactiva:
`jupyter lab` desde esa misma terminal.
