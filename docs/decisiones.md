# Decisiones

Una línea por decisión: fecha · qué · por qué.

- 2026-11-22 · Trabajar con `subset_ADNC_GNN.h5ad` (24.000 células, 80 donantes) para entender los datos y los primeros modelos · subset de desarrollo proporcionado por el grupo MEDAL, manejable en un portátil; un subset mayor se construirá más adelante si hace falta.
- 2026-11-22 · Excluir los 5 donantes Reference de la clasificación principal (quedan 75: Low 12, Intermediate 21, High 42) · son otra cohorte, de 29–60 años, y no tienen datos de patología, APOE4 ni PMI comparables (`docs/datos.md`, §5).
- 2026-11-22 · Hacer el pseudobulk sumando recuentos de `raw`, no desde `X` · `X` es log1p(CP10k) y solo tiene 6.000 genes; `raw` son los recuentos originales con los 36.412 genes (comprobado: X = log1p(raw / total × 10.000)).
- 2026-11-22 · Emparejar los genes de `X` y `raw` por `feature_name` · `X` usa símbolos y `raw` IDs de ENSEMBL; por `feature_name` se encuentran los 6.000 genes y se distinguen los dos `LINC03021`.
- 2026-11-22 · Particiones siempre por donante, estratificadas por clase · la etiqueta es del donante; mezclar células de un donante en train y test infla los resultados.
