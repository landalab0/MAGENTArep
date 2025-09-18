# 1) Descarga de metadatos (ENA + NCBI)

### 1.1 ENA Portal API
```python
# %%
ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"

def ena_query_for_mangrove() -> str:
    # scientific_name preferido + texto libre 'mangrove' como respaldo
    terms = '(scientific_name="mangrove metagenome" OR mangrove)'
    return f'(library_strategy="WGS" AND instrument_platform="ILLUMINA" AND {terms})'

ENA_FIELDS = [
    "run_accession","sample_accession","study_accession",
    "library_strategy","library_layout","instrument_platform",
    "collection_date","country",
    "location","lat","lon",
    "fastq_ftp","fastq_http","fastq_md5"
]

params = {
    "result": "read_run",
    "query":  ena_query_for_mangrove(),
    "fields": ",".join(ENA_FIELDS),
    "format": "tsv",
    "limit":  0
}
url = f"{ENA_BASE}?{urlencode(params)}"
print("URL ENA:", url)

r = requests.get(url, timeout=180, verify=CA_BUNDLE)
r.raise_for_status()

ena_tsv = META_DIR / "ena_mangrove.tsv"
with open(ena_tsv, "wb") as f:
    f.write(r.content)


ena_df = pd.read_csv(ena_tsv, sep="\t", dtype=str)
print("ENA metadatos filas:", len(ena_df), "->", ena_tsv)
```
## ENA Portal API — Manglar  

Consulta el portal de ENA (API) con una query específica para **metagenomas de manglar**. Luego los carga en un DataFrame de pandas para su posterior análisis.  

**Entradas:**  
- Endpoint público de ENA (`https://www.ebi.ac.uk/ena/portal/api/search`)  
- Parámetros de búsqueda (`params`) incluyendo filtro por:  
  - `library_strategy="WGS"`  
  - `instrument_platform="ILLUMINA"`  
  - `scientific_name="mangrove metagenome"` o término libre `mangrove`.  

**Salidas:**  
- Archivo `metadata/ena_mangrove.tsv`  
- DataFrame `ena_df` en memoria con los metadatos cargados.  

**Parámetros editables:**  
- `ena_query_for_mangrove()` (criterio de búsqueda en ENA)  
- `ENA_FIELDS` (campos a descargar, incluye coordenadas `lat`, `lon` y `location`)  
- `limit` (número de resultados, 0 = todos)  
- `timeout` (segundos de espera de la petición HTTP)

### 1.2 NCBI/SRA — EDirect

```python
# %%
query_ncbi = (
    '("WGS"[Strategy] AND "ILLUMINA"[Platform]) AND '
    '( ( mangrove[All Fields] OR mangroves[All Fields] ) '
    '  AND ( metagenome[All Fields] OR metagenomic[All Fields] OR metagenomics[All Fields] ) )'
)

sra_raw_csv = META_DIR / "sra_mangrove_raw_runinfo.csv"
sra_csv     = META_DIR / "sra_mangrove_wgs_illumina.csv"

cmd = 'esearch -db sra -query {q} | efetch -format runinfo > {out}'.format(
    q=shlex.quote(query_ncbi),
    out=shlex.quote(str(sra_raw_csv))
)
print("Comando EDirect:")
print(cmd)

ret = run_cmd(cmd)
need_fallback = (ret != 0) or (not sra_raw_csv.exists()) or (sra_raw_csv.stat().st_size == 0)
if need_fallback:
    print("[SRA][WARN] EDirect falló o sin filas. Intentando fallback RunInfo CSV…")
    url_fb = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?term=" + quote_plus(query_ncbi)
    print("[SRA][fallback] GET", url_fb)
    try:
        rr = requests.get(url_fb, timeout=180, verify=CA_BUNDLE)
        rr.raise_for_status()
        txt = rr.text
        if txt and "Run,ReleaseDate,LoadDate" in txt.splitlines()[0]:
            sra_raw_csv.write_text(txt, encoding="utf-8")
            print("[SRA][fallback] OK ->", sra_raw_csv)
        else:
            print("[SRA][fallback][WARN] Respuesta no luce como RunInfo CSV.")
    except Exception as e:
        print("[SRA][fallback][ERROR]", e)

postfilter_wgs_illumina_csv(sra_raw_csv, sra_csv)
if sra_csv.exists() and sra_csv.stat().st_size > 0:
    print("NCBI RunInfo (filtrado WGS+ILLUMINA) ->", sra_csv)
else:
    print("NCBI RunInfo no disponible.")
```
## NCBI/SRA — EDirect con fallback RunInfo (Manglar)  

Ejecuta una consulta en **NCBI SRA** usando **EDirect** (`esearch` + `efetch`) para obtener metadatos `RunInfo` de corridas relacionadas con **metagenomas de manglar**.  
- Refuerza el filtrado a **WGS + ILLUMINA** mediante un post-procesamiento.  
- Si EDirect falla o devuelve vacío, recurre a un **fallback** descargando el `RunInfo CSV` directamente desde el portal web de NCBI.  

**Entradas:**  
- Herramientas `esearch` y `efetch` disponibles en el `PATH`.  
- Endpoint web de NCBI SRA RunInfo como respaldo.  
- Query definida en `query_ncbi` que incluye:  
  - `"WGS"[Strategy] AND "ILLUMINA"[Platform]`  
  - `(mangrove OR mangroves)`  
  - `(metagenome OR metagenomic OR metagenomics)`  

**Salidas:**  
- `metadata/sra_mangrove_raw_runinfo.csv` → resultado bruto (EDirect o fallback).  
- `metadata/sra_mangrove_wgs_illumina.csv` → resultado filtrado solo con corridas WGS+ILLUMINA.  

**Parámetros editables:**  
- `query_ncbi` (criterio de búsqueda en NCBI).  
- Rutas de salida (`sra_raw_csv`, `sra_csv`).  
- `timeout` (segundos de espera en el fallback HTTP).  
- Número de reintentos en `run_cmd`.

  ### 1.3 Unificación y filtro técnico mínimo

```python
# %%
def normalize_ena(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={
        "run_accession":      "Run",
        "sample_accession":   "Sample",
        "study_accession":    "Study",
        "library_layout":     "LibraryLayout",
        "library_strategy":   "LibraryStrategy",
        "instrument_platform":"Platform",
        "lat":                "latitude",
        "lon":                "longitude",
        "location":           "geographic_location",
    })

def normalize_sra(df: pd.DataFrame) -> pd.DataFrame:
    ren = {}
    if "lat" in df.columns: ren["lat"] = "latitude"
    if "lon" in df.columns: ren["lon"] = "longitude"
    return df.rename(columns=ren)

def _parse_ns_ew(val: str):
    if not isinstance(val, str):
        return None, None
    tokens = re.findall(r'([+-]?\d+(?:\.\d+)?)([NnSsEeWw]?)', val.strip())
    nums = []
    for num, hemi in tokens:
        x = float(num)
        if hemi.upper() == 'S': x = -abs(x)
        elif hemi.upper() == 'N': x =  abs(x)
        elif hemi.upper() == 'W': x = -abs(x)
        elif hemi.upper() == 'E': x =  abs(x)
        nums.append(x)
    if len(nums) >= 2:
        return nums[0], nums[1]
    return None, None

def _extract_lat_lon_from_text(s: str):
    if not isinstance(s, str):
        return np.nan, np.nan
    lat, lon = _parse_ns_ew(s)
    if lat is not None and lon is not None:
        return lat, lon
    nums = re.findall(r'([+-]?\d+(?:\.\d+)?)', s)
    if len(nums) >= 2:
        return float(nums[0]), float(nums[1])
    return np.nan, np.nan

def _valid_lat_lon(lat, lon):
    try:
        if pd.isna(lat) or pd.isna(lon):
            return False
        lat = float(lat); lon = float(lon)
        return np.isfinite(lat) and np.isfinite(lon) and -90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0
    except Exception:
        return False

ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"

def enrich_coords_from_ena_by_runs(run_ids, chunk_size=200):
    """Devuelve DataFrame: Run, latitude, longitude, geographic_location (ENA read_run)."""
    records = []
    run_ids = [r for r in run_ids if isinstance(r, str) and r.strip()]
    for i in range(0, len(run_ids), chunk_size):
        chunk = run_ids[i:i+chunk_size]
        ors = " OR ".join([f'run_accession="{r}"' for r in chunk])
        params = {
            "result": "read_run",
            "query":  ors,
            "fields": ",".join(["run_accession", "lat", "lon", "location"]),
            "format": "tsv",
            "limit":  0,
        }
        url = f"{ENA_BASE}?{urlencode(params)}"
        try:
            r = requests.get(url, timeout=180, verify=CA_BUNDLE)
            r.raise_for_status()
            if r.content:
                dfc = pd.read_csv(io.StringIO(r.content.decode("utf-8")), sep="\t", dtype=str)
                if not dfc.empty:
                    dfc = dfc.rename(columns={
                        "run_accession": "Run",
                        "lat": "latitude",
                        "lon": "longitude",
                        "location": "geographic_location",
                    })
                    records.append(dfc)
        except Exception as e:
            print(f"[ENA][enrich run][WARN] chunk {i//chunk_size}: {e}")
    if not records:
        return pd.DataFrame(columns=["Run","latitude","longitude","geographic_location"])
    return pd.concat(records, ignore_index=True, sort=False)

def enrich_coords_from_ena_by_samples(sample_ids, chunk_size=200):
    """Devuelve DataFrame: Sample, latitude, longitude, geographic_location (ENA read_sample)."""
    records = []
    sample_ids = [s for s in sample_ids if isinstance(s, str) and s.strip()]
    for i in range(0, len(sample_ids), chunk_size):
        chunk = sample_ids[i:i+chunk_size]
        ors = " OR ".join([f'sample_accession="{s}"' for s in chunk])
        params = {
            "result": "sample",
            "query":  ors,
            "fields": ",".join(["sample_accession","lat","lon","location","country"]),
            "format": "tsv",
            "limit":  0,
        }
        url = f"{ENA_BASE}?{urlencode(params)}"
        try:
            r = requests.get(url, timeout=180, verify=CA_BUNDLE)
            r.raise_for_status()
            if r.content:
                dfc = pd.read_csv(io.StringIO(r.content.decode("utf-8")), sep="\t", dtype=str)
                if not dfc.empty:
                    dfc = dfc.rename(columns={
                        "sample_accession": "Sample",
                        "lat": "latitude",
                        "lon": "longitude",
                        "location": "geographic_location",
                    })
                    records.append(dfc)
        except Exception as e:
            print(f"[ENA][enrich sample][WARN] chunk {i//chunk_size}: {e}")
    if not records:
        return pd.DataFrame(columns=["Sample","latitude","longitude","geographic_location"])
    return pd.concat(records, ignore_index=True, sort=False)

def fill_and_filter_geo(all_df: pd.DataFrame) -> pd.DataFrame:
    # Asegura columnas mínimas
    for col in ["latitude","longitude","geographic_location","Sample","Run"]:
        if col not in all_df.columns:
            all_df[col] = np.nan

    # (1) ENA por Run
    need_geo = all_df["Run"].notna() & (all_df["latitude"].isna() | all_df["longitude"].isna())
    runs_to_enrich = all_df.loc[need_geo, "Run"].dropna().astype(str).unique().tolist()
    if runs_to_enrich:
        print(f"[GEO] ENA read_run por Run (n={len(runs_to_enrich)})…")
        enr = enrich_coords_from_ena_by_runs(runs_to_enrich, chunk_size=200)
        if not enr.empty:
            all_df = all_df.merge(enr[["Run","latitude","longitude","geographic_location"]],
                                  on="Run", how="left", suffixes=("", "_ena"))
            for col in ["latitude","longitude","geographic_location"]:
                all_df[col] = all_df[col].fillna(all_df.get(f"{col}_ena"))
                if f"{col}_ena" in all_df.columns:
                    all_df.drop(columns=[f"{col}_ena"], inplace=True)

    # (2) ENA por Sample
    still_missing = (all_df["latitude"].isna() | all_df["longitude"].isna()) & all_df["Sample"].notna()
    samples_to_enrich = all_df.loc[still_missing, "Sample"].dropna().astype(str).unique().tolist()
    if samples_to_enrich:
        print(f"[GEO] ENA sample por Sample (n={len(samples_to_enrich)})…")
        ens = enrich_coords_from_ena_by_samples(samples_to_enrich, chunk_size=200)
        if not ens.empty:
            all_df = all_df.merge(ens[["Sample","latitude","longitude","geographic_location"]],
                                  on="Sample", how="left", suffixes=("", "_sena"))
            for col in ["latitude","longitude","geographic_location"]:
                all_df[col] = all_df[col].fillna(all_df.get(f"{col}_sena"))
                if f"{col}_sena" in all_df.columns:
                    all_df.drop(columns=[f"{col}_sena"], inplace=True)

    # (3) Parseo texto geographic_location
    still_missing2 = all_df["latitude"].isna() | all_df["longitude"].isna()
    if still_missing2.any():
        print(f"[GEO] Parseando geographic_location para {int(still_missing2.sum())} filas…")
        latlon = all_df.loc[still_missing2, "geographic_location"].apply(_extract_lat_lon_from_text)
        lat_parsed = latlon.map(lambda t: t[0])
        lon_parsed = latlon.map(lambda t: t[1])
        all_df.loc[still_missing2, "latitude"]  = all_df.loc[still_missing2, "latitude"].fillna(lat_parsed)
        all_df.loc[still_missing2, "longitude"] = all_df.loc[still_missing2, "longitude"].fillna(lon_parsed)

    # Convertir a numérico y validar
    all_df["latitude"]  = pd.to_numeric(all_df["latitude"], errors="coerce")
    all_df["longitude"] = pd.to_numeric(all_df["longitude"], errors="coerce")
    valid_mask = all_df.apply(lambda r: _valid_lat_lon(r["latitude"], r["longitude"]), axis=1)
    before_n = len(all_df)
    all_df = all_df[valid_mask].copy()
    after_n  = len(all_df)
    print(f"[GEO] Filtrado por coordenadas válidas: {before_n} -> {after_n}")
    return all_df

# ---- Unificación ----
ena_df_norm = pd.DataFrame()
sra_df_norm = pd.DataFrame()

if 'ena_df' in globals() and not ena_df.empty:
    ena_df_norm = normalize_ena(ena_df.copy())
    ena_df_norm["source"] = "ENA"
    ena_df_norm["ecosystem"] = "mangrove"

if (META_DIR / "sra_mangrove_wgs_illumina.csv").exists():
    sra_df = pd.read_csv(META_DIR / "sra_mangrove_wgs_illumina.csv", dtype=str)
    if not sra_df.empty:
        sra_df_norm = normalize_sra(sra_df.copy())
        sra_df_norm["source"] = "SRA"
        sra_df_norm["ecosystem"] = "mangrove"

frames = [x for x in [ena_df_norm, sra_df_norm] if not x.empty]
if not frames:
    print("[WARN] No hay tablas para unir.")
    all_df = pd.DataFrame()
else:
    all_df = pd.concat(frames, ignore_index=True, sort=False)

# Columnas mínimas
for col in ["Run","Sample","Study","LibraryLayout","LibraryStrategy","Platform",
            "latitude","longitude","geographic_location","ecosystem","source"]:
    if col not in all_df.columns:
        all_df[col] = np.nan

# Filtro WGS + ILLUMINA y de-dup por Run
mask = (
    all_df["LibraryStrategy"].fillna("").str.upper().str.contains("WGS") &
    all_df["Platform"].fillna("").str.upper().str.contains("ILLUMINA")
)
all_df = all_df[mask].copy()
all_df["Run"] = all_df["Run"].astype(str)
all_df = all_df.drop_duplicates(subset=["Run"], keep="first")

# Enriquecer coordenadas y filtrar geo válido
if not all_df.empty:
    all_df = fill_and_filter_geo(all_df)

# Exportar
if all_df.empty:
    print("[WARN] Tras filtrar por coordenadas válidas no quedaron corridas.")
    out_eco = OUT_DIR / "mangrove_metadata.tsv"
    save_tsv(pd.DataFrame(), out_eco)
    unified_tsv = META_DIR / "metadatos_unificados.tsv"
    save_tsv(pd.DataFrame(), unified_tsv)
else:
    out_eco = OUT_DIR / "mangrove_metadata.tsv"
    save_tsv(all_df, out_eco)
    unified_tsv = META_DIR / "metadatos_unificados.tsv"
    save_tsv(all_df, unified_tsv)
    print(f"[OK] mangrove: {len(all_df)} filas -> {out_eco}")
    print(f"[OK] Metadatos unificados (manglar, geo-válidos) -> {unified_tsv}")

    summary = (
        all_df.groupby(["ecosystem","source"])["Run"]
        .nunique()
        .reset_index()
        .rename(columns={"Run":"n_runs"})
        .sort_values(["ecosystem","source"])
    )
    log_path = LOGS_DIR / "summary_mangrove_wgs_illumina.tsv"
    save_tsv(summary, log_path)
    print(f"[OK] Resumen -> {log_path}")
```

# Unificación y filtros técnicos de metadatos genómicos

Combina y depura metadatos de **ENA** y **NCBI/SRA** para generar un dataset unificado con coordenadas geográficas validadas.

1. **Normalización** de columnas (`normalize_ena`, `normalize_sra`).
2. **Unificación** de tablas de ENA y SRA.
3. **Filtrado técnico**:  
   - Mantiene solo registros con `LibraryStrategy` metagenómica.  
   - (Opcional) restringe además a `WGS` y `ILLUMINA`.  
4. **Eliminación de duplicados** por `Run`.
5. **Enriquecimiento geográfico**:  
   - Consulta coordenadas en ENA (`read_run`, `sample`).  
   - Parseo de campos de texto (`geographic_location`).  
   - Validación de lat/lon.
6. **Exportación** de resultados en formatos `TSV`.

---

##  Funciones principales

- `normalize_ena(df)`: renombra columnas clave de ENA.  
- `normalize_sra(df)`: renombra columnas clave de SRA.  
- `_parse_ns_ew(val)`: interpreta coordenadas con hemisferio N/S/E/W.  
- `_extract_lat_lon_from_text(s)`: extrae lat/lon desde cadenas libres.  
- `_valid_lat_lon(lat, lon)`: valida rango y formato de coordenadas.  
- `enrich_coords_from_ena_by_runs(run_ids)`: obtiene coordenadas desde ENA por *Run*.  
- `enrich_coords_from_ena_by_samples(sample_ids)`: obtiene coordenadas desde ENA por *Sample*.  
- `fill_and_filter_geo(all_df)`: completa coordenadas faltantes y filtra inválidas.

---

##  Entradas

- `metadata/ena_mangrove.tsv` → tabla ENA.  
- `metadata/sra_mangrove_wgs_illumina.csv` → tabla NCBI/SRA.  

---

##  Salidas

- `metadata/metadatos_unificados.tsv` → tabla consolidada.  
- `output/mangrove_metadata.tsv` → dataset final filtrado por ecosistema.  
- `logs/summary_mangrove_wgs_illumina.tsv` → resumen por `ecosystem` y `source`.  
- `metadata/metadatos_unificados.csv` → versión CSV con todos los metadatos asociados.  

---

##  Parámetros editables

- **Mapeos de columnas** (`rename`).  
- **Filtros aplicados**:  
  - `metagenomic`  
  - `WGS`  
  - `ILLUMINA`
