# 3) Limpieza y enriquecimiento de metadatos (centrado en coordenadas)

## 3.1 Altitud (SRTM)

```python
import elevation

INPUT_TSV = OUT_DIR / "metadatos_unificados.tsv"
OUTPUT_TSV= OUT_DIR / "metadatos_con_altitud.tsv"
BUFFER    = 1.0

df = pd.read_csv(INPUT_TSV, sep="\t")
lat_min, lat_max = df["latitude"].astype(float).min()-BUFFER, df["latitude"].astype(float).max()+BUFFER
lon_min, lon_max = df["longitude"].astype(float).min()-BUFFER, df["longitude"].astype(float).max()+BUFFER

DEM_TIF.parent.mkdir(parents=True, exist_ok=True)
elevation.clip(bounds=(lon_min, lat_min, lon_max, lat_max), output=str(DEM_TIF))
elevation.clean()

vals = []
with rasterio.open(DEM_TIF) as src:
    for _, r in df.iterrows():
        lon, lat = float(r["longitude"]), float(r["latitude"])
        x, y = rio_transform({'init': 'EPSG:4326'}, src.crs, [lon], [lat])
        v = next(src.sample([(x[0], y[0])]))[0]
        vals.append(float(v) if v is not None else np.nan)

df["altitude_m"] = vals
save_tsv(df, OUTPUT_TSV)
print("Altitud añadida ->", OUTPUT_TSV)
```

**Resumen del bloque — SRTM**  
**Qué hace:** recorta DEM y extrae altitud por coordenada.  
**Entradas:** `outputs/metadatos_con_coord.tsv`, DEM `srtm.tif`.  
**Salidas:** `outputs/metadatos_con_altitud.tsv`.  
**Parámetros editables:** `BUFFER`, `DEM_TIF`.
