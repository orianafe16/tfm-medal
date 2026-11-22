# Descripción de los datos: `subset_ADNC_GNN.h5ad`

Generada el 22 nov 2026 con `notebooks/00_carga_y_ficha.ipynb`. Las tablas de origen están en `results/tables/`.

## 1. Lo esencial

- **24.000 células × 6.000 genes** en `X` (log-normalizado) y **36.412 genes** en `raw` (recuentos crudos).
- **80 donantes**: Reference 5 · Low 12 · Intermediate 21 · High 42. Sin Reference quedan **75**.
- **6 tipos celulares** (`cell_type_abbr`): AST, L2IT-GLU, MGC, OLG, SST-GABA, VIP-GABA.
- **Balanceado por células, no por donantes:** exactamente 1.000 células por tipo celular y clase (6.000 por clase).
- El fichero coincide con la descripción del subset en las 10 comprobaciones (`comprobaciones.csv`).

## 2. Estructura y matrices

| Parte | Contenido |
|---|---|
| `X` | CSR float32, 24.000 × 6.000, 17,6 % no ceros, 0,19 GB en RAM. Decimales, máx 6,2 |
| `raw.X` | CSR float32, 24.000 × 36.412, 11,1 % no ceros, 0,73 GB en RAM. Enteros, máx ~3.800, mediana 7.483 cuentas por célula |
| `layers` | ninguna |
| `obsm` | `X_scVI` (20 dimensiones) y `X_umap` (2), calculados por el SEA-AD sobre el dataset completo |
| `obs` | 46 columnas |
| `var` (X) | `feature_name`, `feature_type`, `feature_length`, `n_cells`, `highly_variable`, `means`, `dispersions`... |

Todo junto ocupa ~0,9 GB en RAM, así que con 15,7 GB se puede cargar entero sin modo backed.

**Normalización comprobada:** para 500 células, `X = log1p(raw / total_célula × T)` con **T ≈ 10.000** (9.992, variación entre células del 0,04 %). Los ceros de `X` y de `raw` coinciden al 100 %.

Conclusiones:

- `X` es la normalización estándar CP10k + log1p. El total se calcula sobre los 36.412 genes, no sobre los 6.000, y después se seleccionan los genes.
- `raw` son los recuentos originales: **el pseudobulk se hará desde `raw`**, sumando recuentos.

## 3. Columnas de `obs`

Clasificación en `configs/config.yaml` y `obs_columnas.csv`. No ha quedado ninguna sin clasificar.

| Grupo | Nº | Columnas | Nivel |
|---|---|---|---|
| Identificador | 2 | `donor_id`, `observation_joinid` | donante / célula |
| Etiqueta | 1 | `ADNC` | donante |
| **Patología prohibida** | 8 | Braak, Thal, CERAD, CPS, `Cognitive status`, `disease` (+ ontología), `condition_celltype` | donante (`condition_celltype`: célula) |
| Confusor | 13 | edad (`Age at death`, `development_stage`), sexo, APOE4, PMI, años de educación, etnia, Lewy, LATE-NC, microinfartos | donante |
| Técnica | 11 | `assay`, `Specimen ID`, `Neurotypical reference`, `is_primary_data`, tejido, organismo... | mezcla |
| Calidad por célula | 5 | `n_counts`/`total_counts`, `n_genes`, % mitocondrial (dos versiones) | célula |
| Tipo celular | 6 | `Class` (3), `Subclass` (11), `Supertype` (78), `cell_type`, `cell_type_abbr` (6) | célula |

Notas:

- **`condition_celltype`** es `ADNC_tipo` (p. ej. `High_MGC`): lleva la etiqueta dentro.
- **`Cognitive status` y `disease`** (demencia sí o no) son consecuencias clínicas de la patología. Aun así no coinciden con el ADNC: hay 4 donantes Low con demencia y 11 High sin ella.
- **`Specimen ID`** tiene 206 valores para 80 donantes, así que hay varias muestras físicas por donante. Es el candidato natural a "lote".
- **`assay`** varía dentro del donante: algunos donantes tienen núcleos de 10x 3' v3 y de 10x multiome.
- **`is_primary_data`** es False en parte de las células: CELLxGENE marca así las células que también aparecen en otro dataset. No afecta al análisis.
- `get_donor_table` confirma que **todas** las columnas que deberían ser de donante lo son, y `donantes.csv` queda con 80 filas y 29 columnas (28 metadatos + `n_celulas`).

## 4. Donantes y células

| Clase | Donantes | Células | Células por donante, mediana (mín–máx) |
|---|---|---|---|
| Reference | 5 | 6.000 | 1.242 (143–2.682) |
| Low | 12 | 6.000 | 516 (266–756) |
| Intermediate | 21 | 6.000 | 271 (153–501) |
| High | 42 | 6.000 | 135 (43–275) |

Como cada clase tiene las mismas células, **cuantos más donantes tiene una clase, menos células aporta cada uno**.

### Pares donante × tipo celular (base del pseudobulk)

| Clase | Pares con 0 células | Pares con < 10 | Mediana |
|---|---|---|---|
| Reference | 4 | 2 | 196,5 |
| Low | 0 | 1 | 74 |
| Intermediate | 0 | 2 | 42 |
| High | 4 | 40 | 21 |

- Pares vacíos de High: H21.33.009 (AST), H21.33.020 (SST-GABA), H21.33.034 (OLG) y H21.33.046 (OLG).
- **Riesgo:** los pseudobulk de High suman muchas menos células, así que serán más ruidosos y con menos cuentas que los de Low. Un modelo podría distinguir High por esa diferencia técnica y no por la biología. Habrá que fijar un mínimo de células por par, normalizar el pseudobulk por profundidad y comprobar que el número de células no predice la clase por sí solo.

## 5. Reference

Son 5 donantes (`H18.30.001`, `H18.30.002`, `H19.30.001`, `H19.30.002`, `H200.1023`), todos con `Neurotypical reference = True` y de **29 a 60 años**; el resto tiene 65 o más. Todas sus columnas de patología, APOE4, PMI y educación valen literalmente "Reference", y su CPS es 0.

No solo son otra cohorte: **no tienen datos comparables**. Quedan fuera de la clasificación principal.

## 6. Confusores (75 donantes, sin Reference, contando donantes)

| | Low (12) | Intermediate (21) | High (42) |
|---|---|---|---|
| APOE4 positivo | **0** | 5 | **20** |
| Mujeres | 7 | 10 | **29** |
| 65–77 años | 2 | **0** | 6 |
| 90+ años | 5 | 12 | 24 |
| PMI (3 tramos) | 4 / 4 / 4 | 8 / 6 / 7 | 14 / 19 / 9 |

- **APOE4 es el confusor más fuerte:** ningún Low lo tiene y casi la mitad de los High sí. Es biológicamente esperable, porque es el principal factor de riesgo genético. Pero `APOE` está entre los 6.000 genes de `X`, y un modelo podría apoyarse en él.
- **Sexo:** High tiene más mujeres (69 %).
- **Edad:** `development_stage` da la edad exacta, pero **41 de 75 donantes** aparecen como "80 year-old and over" (dato censurado). Para la mayoría solo sabemos que tienen 80 o más.
- Las copatologías (Lewy, LATE-NC) son más frecuentes en High, como se espera.

## 7. Genes

| | `X` (6.000) | `raw` (36.412) |
|---|---|---|
| Nombre (índice) | símbolo (`APOE`, `GFAP`...) | ID de ENSEMBL sin versión |
| `feature_name` | `SIMBOLO_ENSG...` o `ENSG....versión` | igual |
| Genes sin símbolo | **953 (16 %)**, como `ENSG00000176349.11` | 11.069 |
| Nombres repetidos | 1: `LINC03021` (son dos genes distintos, ENSG00000254319 y ENSG00000277526) | ninguno en el índice |
| Mitocondriales (`MT-`) | **13**, todos entre los HVG | 13 |

- La forma fiable de emparejar `X` con `raw` es **`feature_name`**: encuentra los 6.000 genes.
- Los 6.000 genes son todos HVG (`highly_variable = True`), calculados con el método de dispersión. `n_cells` llega a 27.344, más que las 24.000 células del subset, así que **los HVG se calcularon sobre un conjunto mayor**, no sobre este subset. Encaja con el pipeline del grupo: HVG sobre todas las células del fichero procesado.
- **Genes de Alzheimer:** `APOE` está en `X`; **APP, MAPT, TREM2 y GSK3B no están** entre los 6.000 HVG, pero sí en `raw`.
- Los 13 genes mitocondriales entre los HVG probablemente reflejan calidad técnica más que biología. En núcleos deberían tener muy poca señal.

## 8. Implicaciones para los siguientes pasos

1. **Pseudobulk desde `raw`**, sumando recuentos por donante y tipo celular, con un mínimo de células por par.
2. **Seleccionar genes dentro de la validación cruzada.** Los 6.000 HVG no usan la etiqueta, así que no hay fuga. Pero conviene decidir si se parte de ellos o de `raw` filtrado (para incluir APP, MAPT, TREM2...).
3. **Partición por donante**, estratificada por clase: con 12 donantes Low, cada pliegue de una validación cruzada de 5 tendrá 2–3.
4. **Controlar los confusores:** comprobar que un modelo con solo APOE4, sexo y edad no predice el ADNC tan bien como el de expresión. Valorar quitar `APOE` y los `MT-` de las características.
5. **Comprobar el efecto del número de células:** que `n_celulas` por par no prediga la clase por sí solo.

## 9. Limitaciones conocidas del subset

- Es un subset de desarrollo **balanceado artificialmente** (1.000 células por tipo celular y clase). Las proporciones de tipos celulares y de clases no reflejan las del SEA-AD completo, y los resultados obtenidos aquí deberán validarse en la distribución real.
- Los 6.000 HVG se seleccionaron sobre un conjunto de células mayor que el subset (`n_cells` llega a 27.344) e incluyen 953 genes sin símbolo y 13 mitocondriales.
- No consta cómo se muestrearon las células de cada donante dentro de cada combinación de clase y tipo celular.

## 10. Pseudobulk (`notebooks/01_pseudobulk.ipynb`)

Suma de los recuentos de `raw` por (donante, tipo celular): **472 pseudobulk × 36.412 genes**, guardados en
`Data/processed/pseudobulk.h5ad`. De los 480 pares posibles faltan 8 (donantes sin células de ese tipo).

| Clase | Pseudobulk | Células por pseudobulk (mediana) | Con < 10 células | Recuentos (mediana) |
|---|---|---|---|---|
| Reference | 26 | 216,5 | 2 | 1,70 M |
| Low | 72 | 74 | 1 | 0,93 M |
| Intermediate | 126 | 42 | 2 | 0,57 M |
| High | 248 | 21 | **40** | 0,22 M |

**Genes expresados** (CPM ≥ 1 en al menos la mitad de los pseudobulk del tipo): entre 13.230 (OLG) y 16.632
(L2IT-GLU); 19.633 en al menos un tipo y 11.060 en los seis. De los expresados en cada tipo, solo ~3.600–4.300
están entre los 6.000 HVG de `X`.

**Controles que salen bien:**
- La PCA con todos los pseudobulk separa los seis tipos celulares (PC1 32 %, PC2 19 %): glía frente a neuronas.
- Cada marcador destaca en su tipo: GFAP y AQP4 en AST, P2RY12 y CSF1R en MGC, MOG y PLP1 en OLG, SST en
  SST-GABA, VIP en VIP-GABA, SLC17A7 y CUX2 en L2IT-GLU. La anotación es coherente.
- TREM2 solo se expresa en microglía, y APOE sobre todo en astrocitos y microglía, como se describe en la literatura.

### Efecto técnico ligado a la clase

Como el subset tiene las mismas células por clase repartidas entre más o menos donantes, los pseudobulk de High
resumen menos células y menos recuentos. En **los seis tipos**, la clase (Low < Intermediate < High) se correlaciona
negativamente con el número de células (ρ de Spearman entre −0,50 y −0,73) y con los recuentos totales
(ρ entre −0,50 y −0,71); todas con p < 0,001.

En la PCA dentro de cada tipo celular (sin Reference), **la PC1 está dominada por la profundidad**:

| Tipo | Varianza de PC1 | R² con log(n_celulas) | R² con ADNC |
|---|---|---|---|
| AST | 17 % | 0,18 (PC2: 0,44) | 0,19 |
| L2IT-GLU | 40 % | 0,70 | 0,15 |
| MGC | 14 % | 0,68 | 0,23 |
| OLG | 20 % | 0,81 | 0,13 |
| SST-GABA | 29 % | 0,76 | 0,19 |
| VIP-GABA | 26 % | 0,67 | 0,20 |

La parte de PC1 que explica el ADNC puede deberse, al menos en parte, a su correlación con el número de células.
Por el mismo motivo, las diferencias de expresión entre clases que se ven en algunos genes (p. ej. TREM2, BIN1 y
CLU bajan de Low a High en microglía) **no se pueden interpretar todavía**: con menos recuentos, el log-CPM medio
de los genes poco o moderadamente expresados tiende a bajar.

**Consecuencia:** sin corregir la profundidad, un modelo podría predecir el ADNC a partir del número de células
de cada pseudobulk en lugar de la biología. Es una propiedad del diseño del subset (balanceado por células), no
de los donantes.
