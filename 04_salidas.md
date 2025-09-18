## 4) Salidas clave y consolidación final
```python
final_tsv = OUT_DIR / "metadatos_con_gbif.tsv"
if final_tsv.exists():
    df_final = pd.read_csv(final_tsv, sep="\t")
    df_final.to_csv(META_FILE, index=False)
    print("Escrito estándar:", META_FILE)
else:
    print("No se encontró", final_tsv)
```
**Resumen del bloque — Consolidación**  
**Qué hace:** escribe la última tabla enriquecida al nombre estándar `metadatos_enriquecidos.csv`.  
**Entradas:** `outputs/metadatos_con_gbif.tsv`.  
**Salidas:** `outputs/metadatos_enriquecidos.csv`.  
**Parámetros editables:** `META_BASE`, `META_FILE`.

# Resumen de metadatos generados (esperados)

Este resumen describe **qué columnas** deberían haberse añadido durante el flujo de enriquecimiento, asumiendo que todos los insumos locales y APIs estuvieron disponibles y accesibles.

## Núcleo geo y filtrado
- `latitude`, `longitude` (validadas y acotadas a rangos plausibles).

## Altitud (SRTM / DEM)
- `altitude_m` → **1 columna**  
*Requiere `data/dem/srtm.tif` (o ajustar `DEM_TIF`).*

## Clima (WorldClim 2.1, BIO1–BIO19)
- `bio1` … `bio19` → **19 columnas**  
*Requiere TIFFs `wc2.1_30s_bio_1..19.tif` en `data/clim/` (o ajustar `CLIM_DIR`).*

## Ecorregiones marinas (MEOW)
- `ECOREGION`, `PROVINCE`, `REALM` → **3 columnas**  
*Requiere shapefile `meow_ecos.shp` (+ `.shx/.dbf/.prj`) en `data/marine_ecoregions/MEOW/` (o ajustar `MEOW_SHP`).*

## Huella humana (HII 2020)
- `human_footprint_2020` → **1 columna**  
*Requiere `data/human/hii_2020.tif` (o ajustar `HII_TIF`).*

## Suelos (SoilGrids v2.0, 0–30 cm)
Propiedades: `bdod`, `cec`, `clay`, `sand`, `silt`, `phh2o`, `soc`, `nitrogen`, `cfvo`, `ocs`  
Profundidades: `0-5cm`, `5-15cm`, `15-30cm`  
Estadístico consolidado: `mean`  
- Columnas esperadas: **10 propiedades × 3 profundidades = 30 columnas**  
  Formato: `prop_depth_mean`, p.ej. `clay_0-5cm_mean`.  
> Nota: si la API expone `Q0.5` (mediana), el flujo puede intentar leerla, pero la salida consolidada se guarda en `*_mean` cuando hay valores.

## Biodiversidad (GBIF, radio 2 km)
- `gbif_richness_species_2km` (riqueza específica)  
- `gbif_shannon_H_2km` (diversidad de Shannon)  
- `gbif_pielou_J_2km` (equitatividad de Pielou)  
- `gbif_simpson_1minD_2km` (1 − índice de Simpson)  
- `gbif_chao1_2km` (riqueza estimada, Chao1)  
- `gbif_occurrences_effective` (número de ocurrencias contabilizadas)  
- `gbif_top_species_effective` (Top 5 especies con conteos)  
- `gbif_scope`, `gbif_radius_effective_m`, `gbif_spatial_mode`, `gbif_error` (diagnóstico)  
→ **≈ 10 columnas**

---

## Total **adicional** estimado (si todo estuvo disponible)
- Altitud: 1  
- Clima: 19  
- MEOW: 3  
- HII: 1  
- Suelos: 30  
- GBIF (métricas + diagnóstico): ~10  
**Suma mínima esperada:** **64 columnas** adicionales.

> El conteo real puede variar por: presencia/ausencia de insumos locales, errores de red o cambios en APIs. En tales casos se esperan `NaN` y el flujo debería continuar (salvo pasos marcados como obligatorios por el propio notebook).

## Utilidad
- **Contexto ambiental** (clima, altitud, ecorregión, presión antrópica) por muestra.  
- **Propiedades del suelo** (pH, textura, carbono/nitrógeno, densidad aparente) relevantes para funciones microbianas.  
- **Biodiversidad local** (riqueza/diversidad) para relacionar composición/función metagenómica con entorno biológico.  
- **Comparabilidad y reproducibilidad** al estandarizar metadatos para análisis multi-estudio, conservación y restauración.
