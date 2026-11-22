# Decisiones

Una línea por decisión: fecha · qué · por qué.

- 2026-11-22 · Trabajar con `subset_ADNC_GNN.h5ad` (24.000 células, 80 donantes) para entender los datos y los primeros modelos · subset de desarrollo proporcionado por el grupo MEDAL, manejable en un portátil; un subset mayor se construirá más adelante si hace falta.
- 2026-11-22 · Excluir los 5 donantes Reference de la clasificación principal (quedan 75: Low 12, Intermediate 21, High 42) · son otra cohorte, de 29–60 años, y no tienen datos de patología, APOE4 ni PMI comparables (`docs/datos.md`, §5).
- 2026-11-22 · Hacer el pseudobulk sumando recuentos de `raw`, no desde `X` · `X` es log1p(CP10k) y solo tiene 6.000 genes; `raw` son los recuentos originales con los 36.412 genes (comprobado: X = log1p(raw / total × 10.000)).
- 2026-11-22 · Emparejar los genes de `X` y `raw` por `feature_name` · `X` usa símbolos y `raw` IDs de ENSEMBL; por `feature_name` se encuentran los 6.000 genes y se distinguen los dos `LINC03021`.
- 2026-11-22 · Particiones siempre por donante, estratificadas por clase · la etiqueta es del donante; mezclar células de un donante en train y test infla los resultados.
- 2026-11-22 · Construir el pseudobulk con los 36.412 genes de `raw` y filtrar después los genes expresados por tipo celular (CPM ≥ 1 en ≥ 50 % de los pseudobulk) · los 6.000 HVG de `X` excluyen genes de interés (APP, MAPT, TREM2, GSK3B) y se seleccionaron sobre otro conjunto de células; el filtro por tipo no usa la etiqueta.
- 2026-11-22 · Guardar todos los pseudobulk con su `n_celulas` y aplicar el umbral mínimo de células al modelar · un umbral de 10 eliminaría 40 pseudobulk de High y solo 1 de Low: la decisión depende de la clase y debe evaluarse.
- 2026-11-22 · Mantener a los donantes Reference en el pseudobulk, marcados, y excluirlos al modelar · permite usarlos como control aparte.
- 2026-11-22 · Igualar la profundidad de los pseudobulk por submuestreo binomial dentro de cada tipo celular, con objetivo = cuantil 5 % de `total_counts` y mínimo de 10 células · sin igualar, la profundidad domina la PC1 (R² hasta 0,84) y está ligada a la clase; el cuantil 5 % elimina ese efecto perdiendo 67 de 446 pseudobulk (62 de High).
- 2026-11-22 · Usar el modelo de control técnico (log n_celulas y log total_counts por tipo; balanced accuracy 0,66, kappa 0,71) como listón mínimo de cualquier modelo de expresión · las variables técnicas predicen el ADNC por el diseño del subset.
- 2026-11-22 · Calibrar el "azar" de cada métrica con permutaciones de la etiqueta entre donantes (repitiendo todo el proceso de CV) · con ~75 donantes la CV sin señal no da el valor teórico (AUC ≈ 0,43, Spearman ≈ −0,14 en datos sintéticos).
- 2026-11-22 · No filtrar más células: usar el QC del SEA-AD tal cual · los umbrales por tipo (≥ 1.000–2.000 genes, < 5 % mitocondrial) ya están aplicados; el umbral fijo del repositorio quitaría hasta el 62 % de las neuronas, y el criterio MAD solo marca diferencias mínimas de % mitocondrial (límites 0,4–1,6 %).
- 2026-11-22 · Añadir APOE4 y el % mitocondrial medio por donante a las variables que se controlan antes de interpretar genes · APOE4 es el único confusor asociado al ADNC tras corrección (p BH = 0,047) y el % mitocondrial tiende a ser mayor en High en los 6 tipos.

