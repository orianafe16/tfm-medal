# Resultados

## 1. Baselines sobre pseudobulk igualado (`notebooks/03_baselines.ipynb`, 22 nov 2026)

**Datos:** pseudobulk igualados en profundidad (`docs/datos.md`, §11): 74 donantes (Low 12, Intermediate 21,
High 41), 60–66 por tipo celular.

**Validación:** por donante, 5 pliegues estratificados × 10 repeticiones, mismas particiones para todos los
modelos. Genes expresados, 2.000 genes más variables y escalado se ajustan dentro de cada pliegue.

**Modelos por tipo celular:** trampa (log n_celulas y log recuentos originales del tipo, logística); regresión
logística L2 sobre log-CPM (C elegida por CV interna entre 0,001 y 1); random forest (500 árboles).

### Por tipo celular (balanced accuracy, media de 10 repeticiones; azar = 0,333)

| Tipo | Donantes | Trampa (técnicas) | Logística (expresión) | Random forest (expresión) |
|---|---|---|---|---|
| AST | 64 | 0,501 | 0,388 | 0,371 |
| L2IT-GLU | 63 | 0,657 | 0,385 | 0,357 |
| MGC | 63 | 0,531 | 0,350 | 0,333 |
| OLG | 60 | 0,462 | 0,327 | 0,351 |
| SST-GABA | 63 | 0,519 | 0,384 | 0,354 |
| VIP-GABA | 66 | 0,565 | 0,354 | 0,385 |

Kappa cuadrático de los modelos de expresión: entre −0,06 y 0,20 (azar = 0); trampa: entre 0,45 y 0,71.
Desviaciones entre repeticiones: 0,02–0,04.

### Global (74 donantes)

| Modelo | Balanced accuracy | Kappa cuadrático |
|---|---|---|
| Azar (siempre High) | 0,333 | 0,000 |
| Trampa global (técnicas de los 6 tipos) | 0,638 ± 0,023 | 0,685 |
| Fusión de trampas por tipo | 0,687 ± 0,016 | 0,760 |
| **Fusión logística (expresión)** | **0,385 ± 0,022** | **0,115** |
| **Fusión random forest (expresión)** | **0,378 ± 0,023** | **0,119** |

**Confusión de la fusión logística** (10 repeticiones sumadas): de 120 predicciones de donantes Low, 0 aciertos
(29 como Intermediate, 91 como High); Intermediate, 52 de 210; High, 372 de 410. El modelo predice casi siempre la
clase mayoritaria.

**Comprobación de fuga indirecta:** correlación de Spearman, dentro de cada clase, entre P(High) de la fusión
logística y log(n_celulas) del donante: Low −0,25 (p = 0,44), Intermediate −0,08 (p = 0,74), **High −0,44
(p = 0,004)**. Dentro de High, el modelo de expresión asigna más probabilidad de High a los donantes con menos
células: la poca señal que tiene parece técnica.

### Interpretación

1. Con este subset y estos modelos, **la expresión de los pseudobulk no predice el ADNC de forma útil**
   (balanced accuracy 0,33–0,39, kappa ≤ 0,20). La única información predictiva es la del diseño del muestreo.
2. Explicaciones posibles, no excluyentes:
   - **El balanceo elimina la señal más conocida.** Gabitto et al. (2024) describen sobre todo cambios de
     **composición** con la progresión (pérdida de interneuronas Sst, y después de neuronas excitadoras de capas
     superficiales y Pvalb). Con 1.000 células por tipo y clase, la composición queda fijada por diseño.
   - **Poca información por pseudobulk.** Igualar al cuantil 5 % reduce mucho la profundidad de Low e Intermediate
     (p. ej. Low pasa de ~0,9 M de recuentos de mediana a 70–320 mil según el tipo), y los pseudobulk de High
     resumen ~20 células.
   - **Pocos ejemplos frente a muchas variables:** ~50 donantes de entrenamiento por tipo (≈10 Low) y 2.000 genes.
     Los cambios de expresión dentro de un tipo celular asociados al ADNC son pequeños.
3. Estos resultados **no indican que no haya señal de ADNC en la expresión**: indican que este subset no permite
   medirla con modelos sencillos. Son la referencia que deberán superar los modelos siguientes (incluida la GNN).

## 2. Comprobaciones: Low frente a High y regresión del CPS (`notebooks/04_comprobaciones.ipynb`, 22 nov 2026)

Mismos datos, preprocesado dentro del pliegue y validación por donante (5 pliegues × 10 repeticiones) que en §1.
El CPS se usa solo como **objetivo**, nunca como característica.

**Referencia de azar por permutación.** Con pocos donantes, la validación cruzada sin señal no da el valor teórico
de azar (en datos sintéticos sin señal: AUC ≈ 0,43, Spearman ≈ −0,14). Por eso se barajó la etiqueta entre
donantes 20 veces, repitiendo todo el proceso. Con 20 permutaciones, el menor p-valor posible es 0,048.

### Low (12) frente a High (41): AUC ROC

| Tipo | Trampa (técnicas) | Logística (expresión) |
|---|---|---|
| AST | 0,867 | 0,704 ± 0,057 |
| L2IT-GLU | 0,977 | 0,598 ± 0,048 |
| MGC | 0,870 | 0,593 ± 0,064 |
| OLG | 0,897 | 0,468 ± 0,053 |
| SST-GABA | 0,895 | 0,636 ± 0,046 |
| VIP-GABA | 0,928 | 0,727 ± 0,041 |
| **Fusión** | **0,986** | **0,702 ± 0,033** |

Permutación (fusión de expresión): nulo medio 0,500, percentil 95 0,633; observado 0,702 (**p ≤ 0,048**). La
expresión separa Low de High por encima del azar, pero muy por debajo del control técnico, que separa los extremos
casi perfectamente.

### Regresión del CPS (74 donantes): Spearman entre CPS predicho y real

| Tipo | Trampa (técnicas) | Ridge (expresión) |
|---|---|---|
| **AST** | 0,099 | **0,559 ± 0,026** |
| L2IT-GLU | 0,355 | 0,443 ± 0,044 |
| **MGC** | 0,297 | **0,508 ± 0,029** |
| OLG | 0,451 | 0,480 ± 0,048 |
| SST-GABA | 0,400 | 0,447 ± 0,033 |
| VIP-GABA | 0,399 | 0,474 ± 0,051 |
| **Fusión** | 0,478 (R² 0,20) | **0,579 ± 0,029 (R² 0,25)** |

Permutación (fusión ridge): nulo medio −0,084, percentil 95 0,130; observado 0,579 (**p ≤ 0,048**).
El CPS se correlaciona con el número total de células del donante (Spearman −0,62, p = 5·10⁻⁹), lo que explica
que la trampa también lo prediga en parte.

### Interpretación

1. **Hay señal de progresión en la expresión.** Con un objetivo continuo, el modelo de expresión predice el CPS muy
   por encima del azar por permutación y **supera al control técnico** (0,58 frente a 0,48 en la fusión).
2. **Astrocitos y microglía destacan.** En AST la expresión predice el CPS (0,56) mientras que sus variables
   técnicas casi no lo hacen (0,10); en MGC, 0,51 frente a 0,30. En los demás tipos, expresión y trampa quedan
   cerca. Es coherente con la activación inflamatoria de astrocitos y microglía que Gabitto et al. (2024)
   describen en la fase temprana de la progresión, pero hay que confirmarlo.
3. **Por qué la clasificación de §1 no lo veía:** con 3 clases desequilibradas (12 Low), el modelo acaba
   prediciendo la clase mayoritaria; el CPS aprovecha la variación **dentro** de cada clase (p. ej. Low va de
   0,15 a 0,64) y los 74 donantes como valores continuos.
4. **Pendiente:** comprobar que la señal de expresión es independiente del número de células (expresión +
   técnicas frente a técnicas solas, o correlación parcial controlando log n_celulas). En AST la trampa no
   predice el CPS, así que ahí el riesgo es menor.
