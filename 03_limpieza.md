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

## 3.3 WorldClimb 2.1 (BIO1–BIO19)
```python
INPUT_TSV  = OUT_DIR / "metadatos_con_altitud.tsv"
OUTPUT_TSV = OUT_DIR / "metadatos_con_clima.tsv"

df = pd.read_csv(INPUT_TSV, sep="\t")
for i in range(1, 20):
    df[f"bio{i}"] = np.nan

def extract_from_raster(src, lon, lat):
    try:
        row, col = src.index(lon, lat)
        return float(src.read(1)[row, col])
    except Exception:
        return np.nan

for i in range(1, 19+1):
    tif = CLIM_DIR / f"wc2.1_30s_bio_{i}.tif"
    if not tif.exists():
        print("Falta", tif, "omitiendo bio", i)
        continue
    with rasterio.open(tif) as src:
        for idx,(lon,lat) in enumerate(tqdm(zip(df["longitude"].astype(float), df["latitude"].astype(float)), total=len(df), desc=f"bio{i}")):
            v = extract_from_raster(src, float(lon), float(lat))
            if not np.isfinite(v):
                arr = src.read(1)
                r,c = src.index(float(lon), float(lat))
                r0, r1 = max(0, r-10), min(arr.shape[0], r+11)
                c0, c1 = max(0, c-10), min(arr.shape[1], c+11)
                window = arr[r0:r1, c0:c1]
                valid = window[np.isfinite(window)]
                v = float(valid[0]) if valid.size>0 else np.nan
            df.at[idx, f"bio{i}"] = v

save_tsv(df, OUTPUT_TSV)
print("WorldClim añadido ->", OUTPUT_TSV)

```
**Resumen del bloque — WorldClim**  
**Qué hace:** extrae BIO1–BIO19 por coordenada con fallback local.  
**Entradas:** `outputs/metadatos_con_altitud.tsv`, `wc2.1_30s_bio_*.tif`.  
**Salidas:** `outputs/metadatos_con_clima.tsv`.  
**Parámetros editables:** `CLIM_DIR`, vecindad.

## 3.4 Ecorregiones marinas (MEOW)
```python
INPUT_TSV  = OUT_DIR / "metadatos_con_clima.tsv"
OUTPUT_TSV = OUT_DIR / "metadatos_con_ecoregion_marina.tsv"

df = pd.read_csv(INPUT_TSV, sep="\t")
gdf_points = gpd.GeoDataFrame(
    df.copy(),
    geometry=[Point(xy) for xy in zip(df["longitude"].astype(float), df["latitude"].astype(float))],
    crs="EPSG:4326"
)

if not MEOW_SHP.exists():
    raise FileNotFoundError(f"No se encontró el shapefile MEOW en {MEOW_SHP}. Ajuste la ruta MEOW_SHP.")
gdf_meow = gpd.read_file(MEOW_SHP).to_crs("EPSG:4326")

joined = gpd.sjoin(
    gdf_points,
    gdf_meow[["ECOREGION", "PROVINCE", "REALM", "geometry"]],
    how="left",
    predicate="intersects"
)

out = joined.drop(columns="geometry")
save_tsv(out, OUTPUT_TSV)
print("MEOW añadido ->", OUTPUT_TSV)

```

**Resumen del bloque — MEOW**  
**Qué hace:** crea `GeoDataFrame` de puntos y une con MEOW para asignar ecorregión/provincia/reino.  
**Entradas:** `outputs/metadatos_con_clima.tsv`, shapefile `MEOW`.  
**Salidas:** `outputs/metadatos_con_ecoregion_marina.tsv`.  
**Parámetros editables:** `MEOW_SHP`, `predicate`.

## 3.5 Huella humana (HII 2020)
```python
INPUT_TSV = OUT_DIR / "metadatos_con_ecoregion_marina.tsv"
OUTPUT_TSV= OUT_DIR / "metadatos_con_hii.tsv"

if not HII_TIF.exists():
    raise FileNotFoundError(f"No se encontró el raster HII en {HII_TIF}. Ajuste la ruta HII_TIF.")

df = pd.read_csv(INPUT_TSV, sep="\t")
vals = []
with rasterio.open(HII_TIF) as src:
    arr = src.read(1)
    nodata = src.nodata
    for lon, lat in tqdm(zip(df["longitude"].astype(float), df["latitude"].astype(float)), total=len(df), desc="Extrayendo HII"):
        try:
            r, c = src.index(float(lon), float(lat))
            if 0 <= r < arr.shape[0] and 0 <= c < arr.shape[1]:
                v = arr[r, c]
            else:
                v = np.nan
            if (nodata is not None and v == nodata) or not np.isfinite(v):
                r0, r1 = max(0, r-10), min(arr.shape[0], r+11)
                c0, c1 = max(0, c-10), min(arr.shape[1], c+11)
                sub = arr[r0:r1, c0:c1]
                valid = sub[(sub != nodata) & np.isfinite(sub)]
                v = valid[0] if valid.size > 0 else np.nan
            vals.append(float(v) if np.isfinite(v) else np.nan)
        except Exception:
            vals.append(np.nan)

df["human_footprint_2020"] = vals
save_tsv(df, OUTPUT_TSV)
print("HII añadido ->", OUTPUT_TSV)
```
**Resumen del bloque — HII**  
**Qué hace:** extrae huella humana para cada coordenada con manejo de NoData y fallback local.  
**Entradas:** `outputs/metadatos_con_ecoregion_marina.tsv`, raster `hii_2020.tif`.  
**Salidas:** `outputs/metadatos_con_hii.tsv`.  
**Parámetros editables:** `HII_TIF`, tamaño de vecindad.

## 3.6 Suelos — SoilGrids v2.0 (0–30 cm) + WRB
```python
import time
import math
from typing import Dict, Tuple, Any, Optional, List
import requests
import numpy as np
import pandas as pd
from tqdm import tqdm

INPUT_TSV  = OUT_DIR / "metadatos_con_hii.tsv"
OUTPUT_TSV = OUT_DIR / "metadatos_con_suelos.tsv"

ALL_PROPS = ["bdod", "cec", "clay", "sand", "silt", "phh2o", "soc", "nitrogen", "cfvo", "ocs"]
DEPTHS    = ["0-5cm","5-15cm","15-30cm"]
VALUES    = ["mean","Q0.5"]

REQUESTS_PER_SEC = 6
SLEEP_BETWEEN = 1.0 / REQUESTS_PER_SEC

MAX_RETRIES = 4
RETRY_BACKOFF = [1.0, 2.0, 4.0, 8.0]

SOIL_ENDPOINTS = [
    "https://rest.isric.org/soilgrids/v2.0/properties/query",
    "https://api.soilgrids.org/soilgrids/v2.0/properties/query",
]

CONNECT_TIMEOUT_S = 8
READ_TIMEOUT_S = 20


RADIUSES_KM = [0.0, 0.5, 1.0, 1.5, 2.0]
DIRECTIONS_DEG = [0, 45, 90, 135, 180, 225, 270, 315]

ROUND_DECIMALS: Optional[int] = None  # p.ej., 4 => ~11 m

def move_point(lat: float, lon: float, distance_km: float, bearing_deg: float) -> Tuple[float, float]:
    R = 6371.0
    b = math.radians(bearing_deg)
    lat1 = math.radians(lat); lon1 = math.radians(lon)
    d = distance_km / R
    lat2 = math.asin(math.sin(lat1)*math.cos(d) + math.cos(lat1)*math.sin(d)*math.cos(b))
    lon2 = lon1 + math.atan2(math.sin(b)*math.sin(d)*math.cos(lat1), math.cos(d) - math.sin(lat1)*math.sin(lat2))
    return (math.degrees(lat2), (math.degrees(lon2) + 540) % 360 - 180)

def pick_endpoint() -> str:
    return SOIL_ENDPOINTS[int(time.time()) % len(SOIL_ENDPOINTS)]

_last_request_time = 0.0
def throttle():
    global _last_request_time
    now = time.time()
    wait = SLEEP_BETWEEN - (now - _last_request_time)
    if wait > 0:
        time.sleep(wait)
    _last_request_time = time.time()

def build_params(lat: float, lon: float, properties: List[str]) -> Dict[str, Any]:
    return {
        "lat": float(lat),
        "lon": float(lon),
        "property": ",".join(properties),
        "depth": ",".join(DEPTHS),
        "value": ",".join(VALUES),
        "prefix": "mean",
    }

def do_request(params: Dict[str, Any]) -> Optional[requests.Response]:
    for attempt in range(MAX_RETRIES):
        throttle()
        try:
            resp = requests.get(pick_endpoint(), params=params, timeout=(CONNECT_TIMEOUT_S, READ_TIMEOUT_S))
            if resp.status_code == 200:
                return resp
            if 500 <= resp.status_code < 600:
                time.sleep(RETRY_BACKOFF[min(attempt, len(RETRY_BACKOFF)-1)])
                continue
            return resp
        except requests.RequestException:
            time.sleep(RETRY_BACKOFF[min(attempt, len(RETRY_BACKOFF)-1)])
    return None

def _depth_variants(depth_label: str) -> List[str]:
    return [depth_label, depth_label.replace("-", "_")]

def any_numeric_value(props: Dict[str, Any], prop: str, depth_label: str) -> bool:
    for d in _depth_variants(depth_label):
        for stat in VALUES:
            v = props.get(f"{prop}_{d}_{stat}")
            if isinstance(v, (int, float)):
                return True
    return False

def extract_value(props: Dict[str, Any], prop: str, depth_label: str) -> Optional[float]:

    for d in _depth_variants(depth_label):
        for stat in VALUES:
            v = props.get(f"{prop}_{d}_{stat}")
            if isinstance(v, (int, float)):
                return float(v)
    return None

def probe_cell_has_numeric(lat: float, lon: float) -> bool:
    params = build_params(lat, lon, ["clay"])
    resp = do_request(params)
    if resp is None or resp.status_code != 200:
        return False
    try:
        data = resp.json()
        props = data.get("properties", {})
    except Exception:
        return False
    for d in DEPTHS:
        if any_numeric_value(props, "clay", d):
            return True
    return False

def fetch_all_properties(lat: float, lon: float, props_list: List[str]) -> Optional[Dict[str, Any]]:
    params = build_params(lat, lon, props_list)
    resp = do_request(params)
    if resp is None or resp.status_code != 200:
        return None
    try:
        return resp.json()
    except Exception:
        return None

def enrich_coord_with_fallback(lat: float, lon: float, props_list: List[str]) -> Optional[Dict[str, Any]]:
    for radius in RADIUSES_KM:
        candidates = [(lat, lon)] if radius == 0 else [move_point(lat, lon, radius, b) for b in DIRECTIONS_DEG]
        for lt, ln in candidates:
            if not (-90 <= lt <= 90 and -180 <= ln <= 180):
                continue
            if not probe_cell_has_numeric(lt, ln):
                continue
            data = fetch_all_properties(lt, ln, props_list)
            if data is not None:
                return data
    return None

df = pd.read_csv(INPUT_TSV, sep="\t").reset_index(drop=True)


for p in ALL_PROPS:
    for d in DEPTHS:
        col = f"{p}_{d}_mean"
        if col not in df.columns:
            df[col] = np.nan


def coord_key(lat: float, lon: float) -> Tuple[float, float]:
    if ROUND_DECIMALS is not None:
        return (round(lat, ROUND_DECIMALS), round(lon, ROUND_DECIMALS))
    return (lat, lon)

cache: Dict[Tuple[float, float], Optional[Dict[str, Any]]] = {}

lats = df["latitude"].astype(float).to_numpy()
lons = df["longitude"].astype(float).to_numpy()

unique_keys = []
seen = set()
for lt, ln in zip(lats, lons):
    k = coord_key(lt, ln)
    if k not in seen:
        seen.add(k)
        unique_keys.append(k)

for (lt, ln) in tqdm(unique_keys, desc="SoilGrids (props/dep/values actualizados)"):
    cache[(lt, ln)] = enrich_coord_with_fallback(lt, ln, ALL_PROPS)


for idx, (lt, ln) in enumerate(zip(lats, lons)):
    k = coord_key(lt, ln)
    data = cache.get(k)
    if not data:
        continue
    props_dict = data.get("properties", {}) if isinstance(data, dict) else {}
    for p in ALL_PROPS:
        for d in DEPTHS:
            col = f"{p}_{d}_mean"
            v = extract_value(props_dict, p, d)
            if v is not None:
                df.at[idx, col] = v

save_tsv(df, OUTPUT_TSV)
print("SoilGrids añadido ->", OUTPUT_TSV)
```

**Resumen del bloque — SoilGrids + WRB**  
**Qué hace:** consulta SoilGrids por punto con control de tasa y reintentos; agrega propiedades y WRB si disponible.  
**Entradas:** `outputs/metadatos_con_hii.tsv`, API SoilGrids.  
**Salidas:** `outputs/metadatos_con_suelos.tsv`.  
**Parámetros editables:** `PROPERTIES`, `DEPTHS`, `MAX_WORKERS`, `REQUESTS_PER_SEC`, `RETRIES`, `BACKOFF`.

## 3.7 GBIF — especies de manglar y biodiversidad de plantas
```python
from __future__ import annotations
import time, math, threading
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter

import numpy as np
import pandas as pd
import requests

INPUT_TSV  = OUT_DIR / "metadatos_con_suelos.tsv"
OUTPUT_TSV = OUT_DIR / "metadatos_con_gbif.tsv"

RADIUS_M = 2000
KINGDOM_KEY = None
YEAR_FROM = None
YEAR_TO   = None

GBIF_OCC_URL   = "https://api.gbif.org/v1/occurrence/search"
GBIF_TAXON_URL = "https://api.gbif.org/v1/species/"


PAGE_LIMIT = 300
MAX_PAGES_PER_POINT = 50
MAX_RECORDS_PER_POINT = PAGE_LIMIT * MAX_PAGES_PER_POINT


MAX_WORKERS = 8
GLOBAL_RPS = 3
MAX_BACKOFF_SEC = 60

_print_lock = threading.Lock()
def log(msg: str):
    with _print_lock:
        print(msg, flush=True)

class RateLimiter:
    def __init__(self, rate_per_sec: float, capacity: int | None = None):
        self.rate = float(rate_per_sec)
        self.capacity = int(capacity or max(1, int(rate_per_sec * 2)))
        self.tokens = self.capacity
        self.lock = threading.Lock()
        self.ts = time.time()
    def acquire(self):
        while True:
            with self.lock:
                now = time.time()
                refill = (now - self.ts) * self.rate
                if refill > 0:
                    self.tokens = min(self.capacity, self.tokens + refill)
                    self.ts = now
                if self.tokens >= 1:
                    self.tokens -= 1
                    return
            time.sleep(0.01)

RATE = RateLimiter(GLOBAL_RPS)

def make_session() -> requests.Session:
    s = requests.Session()
    ad = requests.adapters.HTTPAdapter(pool_connections=MAX_WORKERS*2, pool_maxsize=MAX_WORKERS*2, max_retries=0)
    s.mount("https://", ad); s.mount("http://", ad)
    s.headers.update({"User-Agent": "MAGENTA-GBIF-diversity/2km/1.0 (+local)"})
    return s

SESSION = make_session()

def _get(url: str, params: Dict[str, Any], max_attempts=8) -> Optional[requests.Response]:
    """GET con RPS global y backoff exponencial (+ Retry-After cuando esté presente)."""
    attempt = 0; backoff = 1.0
    while attempt < max_attempts:
        attempt += 1
        RATE.acquire()
        try:
            r = SESSION.get(url, params=params, timeout=30)
        except requests.RequestException:
            time.sleep(min(backoff, MAX_BACKOFF_SEC))
            backoff = min(MAX_BACKOFF_SEC, backoff * 2) + 0.1 * attempt
            continue
        if r.status_code == 200:
            return r
        if r.status_code in (429, 500, 502, 503, 504):
            ra = r.headers.get("Retry-After")
            wait = float(ra) if ra and str(ra).isdigit() else backoff
            time.sleep(min(wait, MAX_BACKOFF_SEC))
            backoff = min(MAX_BACKOFF_SEC, backoff * 2) + 0.1 * attempt
            continue

        return None
    return None

_taxon_cache: Dict[int, str | None] = {}

def gbif_species_name(species_key: int) -> Optional[str]:
    """Resuelve nombre científico desde speciesKey, con caché simple."""
    if species_key in _taxon_cache:
        return _taxon_cache[species_key]
    r = _get(f"{GBIF_TAXON_URL}{species_key}", {})
    if not r:
        _taxon_cache[species_key] = None
        return None
    try:
        js = r.json()
        name = js.get("scientificName") or js.get("canonicalName")
    except Exception:
        name = None
    _taxon_cache[species_key] = name
    return name

def _apply_common_filters(q: Dict[str, Any]):
    q["hasCoordinate"] = "true"
    q["hasGeospatialIssue"] = "false"
    q["occurrenceStatus"] = "PRESENT"
    if KINGDOM_KEY is not None:
        q["kingdomKey"] = KINGDOM_KEY
    if YEAR_FROM is not None:
        q["year"] = f"{YEAR_FROM},{YEAR_TO or ''}".strip(",")

def _bbox_wkt(lat: float, lon: float, radius_m: int) -> str:
    """BBox aproximado en grados a partir de radio en metros (depende de latitud)."""
    dlat = radius_m / 1000.0 / 111.32
    dlon = dlat / max(0.1, math.cos(math.radians(lat)))
    minx = lon - dlon; maxx = lon + dlon
    miny = lat - dlat; maxy = lat + dlat
    return f"POLYGON(({minx:.6f} {miny:.6f},{maxx:.6f} {miny:.6f},{maxx:.6f} {maxy:.6f},{minx:.6f} {maxy:.6f},{minx:.6f} {miny:.6f}))"

def _page_species_counts(params: Dict[str, Any]) -> Tuple[int, Counter]:
    """
    Pagina occurrence/search sumando speciesKey (o taxonKey) en cliente.
    Devuelve (total_reported, species_counter).
    """
    p0 = dict(params); p0["limit"] = 0
    total_reported = 0
    r0 = _get(GBIF_OCC_URL, p0)
    if r0:
        try:
            total_reported = int(r0.json().get("count", 0))
        except Exception:
            total_reported = 0

    species_counter: Counter = Counter()
    fetched = 0
    offset = 0
    pages = 0

    while fetched < min(total_reported, MAX_RECORDS_PER_POINT) and pages < MAX_PAGES_PER_POINT:
        q = dict(params)
        q["limit"] = PAGE_LIMIT
        q["offset"] = offset
        r = _get(GBIF_OCC_URL, q)
        if not r:
            break
        try:
            js = r.json()
        except Exception:
            break
        results = js.get("results", [])
        if not results:
            break
        for rec in results:
            sk = rec.get("speciesKey") or rec.get("taxonKey")
            if isinstance(sk, int):
                species_counter[sk] += 1
        n = len(results)
        fetched += n
        offset += n
        pages += 1
        if n < PAGE_LIMIT:
            break

    return total_reported, species_counter

def gbif_occurrence_counts_spatial(lat: float, lon: float, radius_m: int) -> Tuple[int, List[Tuple[int,int]], str]:
    """
    Devuelve (total_occurrences, [(speciesKey, count), ...], spatial_mode)
      A) geo_distance  (snake_case)
      B) geoDistance   (camelCase)
      C) geometry=BBOX (WKT)
    """
    baseA = {}
    _apply_common_filters(baseA)
    baseA["geo_distance"] = f"{radius_m}m,{lat:.6f},{lon:.6f}"
    total, counter = _page_species_counts(baseA)
    if total > 0 and counter:
        items = sorted(counter.items(), key=lambda x: x[1], reverse=True)
        return total, items, "geo_distance"

    baseB = {}
    _apply_common_filters(baseB)
    baseB["geoDistance"] = f"{radius_m}m,{lat:.6f},{lon:.6f}"
    total, counter = _page_species_counts(baseB)
    if total > 0 and counter:
        items = sorted(counter.items(), key=lambda x: x[1], reverse=True)
        return total, items, "geoDistance"

    baseC = {}
    _apply_common_filters(baseC)
    baseC["geometry"] = _bbox_wkt(lat, lon, radius_m)
    total, counter = _page_species_counts(baseC)
    items = sorted(counter.items(), key=lambda x: x[1], reverse=True)
    return total, items, "bbox_wkt" if total > 0 else "none"

def diversity_metrics(species_counts: List[int], label_suffix: str) -> Dict[str, float]:
    out = {f"gbif_richness_species_{label_suffix}": 0.0,
           f"gbif_shannon_H_{label_suffix}": np.nan,
           f"gbif_pielou_J_{label_suffix}": np.nan,
           f"gbif_simpson_1minD_{label_suffix}": np.nan,
           f"gbif_chao1_{label_suffix}": np.nan}
    S = int(sum(1 for c in species_counts if c > 0))
    out[f"gbif_richness_species_{label_suffix}"] = float(S)
    N = float(sum(species_counts))
    if S == 0 or N <= 0:
        return out
    p = [c / N for c in species_counts if c > 0]
    H = -sum(pi * math.log(pi) for pi in p)
    out[f"gbif_shannon_H_{label_suffix}"] = float(H)
    out[f"gbif_pielou_J_{label_suffix}"] = float(H / math.log(S)) if S > 1 else np.nan
    D = sum(pi * pi for pi in p)
    out[f"gbif_simpson_1minD_{label_suffix}"] = float(1.0 - D)
    F1 = sum(1 for c in species_counts if c == 1)
    F2 = sum(1 for c in species_counts if c == 2)
    chao1 = (S + (F1 * F1) / (2.0 * F2)) if F2 > 0 else (S + (F1 * (F1 - 1)) / 2.0)
    out[f"gbif_chao1_{label_suffix}"] = float(chao1)
    return out

def fetch_point_fixed_2km(lat: float, lon: float) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "gbif_occurrences_effective": 0,
        "gbif_top_species_effective": None,
        "gbif_scope": ("All" if KINGDOM_KEY is None else str(KINGDOM_KEY)),
        "gbif_radius_effective_m": RADIUS_M,
        "gbif_spatial_mode": None,
        "gbif_error": None
    }

    total, items, spatial_mode = gbif_occurrence_counts_spatial(lat, lon, RADIUS_M)
    out["gbif_spatial_mode"] = spatial_mode
    out["gbif_occurrences_effective"] = int(total) if total else 0

    counts = [c for _, c in items] if items else []
    out.update(diversity_metrics(counts, "2km"))

    top = (items or [])[:5]
    names = []
    for sk, c in top:
        nm = gbif_species_name(sk)
        if nm:
            names.append(f"{nm} ({c})")
    if names:
        out["gbif_top_species_effective"] = "; ".join(names)

    return out

def _timed_fetch(lat: float, lon: float):
    t0 = time.time()
    out = fetch_point_fixed_2km(lat, lon)
    return out, time.time() - t0

df = pd.read_csv(INPUT_TSV, sep="\t").reset_index(drop=True)
if not {"latitude","longitude"}.issubset(df.columns):
    raise ValueError("Se requieren columnas 'latitude' y 'longitude' en el TSV de entrada.")

idxs = [i for i in range(len(df)) if pd.notnull(df.at[i, "latitude"]) and pd.notnull(df.at[i, "longitude"])]

features_by_idx: Dict[int, Dict[str, Any]] = {}
completed = 0
total = len(idxs)
log(f"GBIF diversidad (ALL taxa) **FIJO 2 km** | filas={total} | workers={MAX_WORKERS} | rps={GLOBAL_RPS} | radius={RADIUS_M} m")

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
    futs = {ex.submit(_timed_fetch, float(df.at[i, "latitude"]), float(df.at[i, "longitude"])): i for i in idxs}
    for fut in as_completed(futs):
        i = futs[fut]
        try:
            row, elapsed = fut.result()
            ok = "OK"
        except Exception as e:
            row = {"gbif_error": str(e),
                   "gbif_occurrences_effective": 0,
                   "gbif_scope": ("All" if KINGDOM_KEY is None else str(KINGDOM_KEY)),
                   "gbif_radius_effective_m": RADIUS_M,
                   "gbif_spatial_mode": "error"}
            ok = f"ERR:{e}"
            elapsed = np.nan

        features_by_idx[i] = row
        lat = float(df.at[i, "latitude"]); lon = float(df.at[i, "longitude"])
        sfx = "2km"
        log(f"[{completed+1}/{total}] lat={lat:.5f} lon={lon:.5f} "
            f"S={row.get(f'gbif_richness_species_{sfx}', 0)} "
            f"H'={row.get(f'gbif_shannon_H_{sfx}', np.nan)} "
            f"occ={row.get('gbif_occurrences_effective', 0)} "
            f"r_eff={row.get('gbif_radius_effective_m')} mode={row.get('gbif_spatial_mode')} {ok} | {elapsed:.2f}s")
        completed += 1

feat_df = pd.DataFrame.from_dict(features_by_idx, orient="index").reindex(df.index)
out_df = pd.concat([df, feat_df], axis=1)
save_tsv(out_df, OUTPUT_TSV)
print("GBIF añadido ->", OUTPUT_TSV)
```
**Resumen del bloque — GBIF (manglar y plantas)**  
**Qué hace:** calcula riqueza de plantas y riqueza/top-5 de especies de manglar por coordenada.  
**Entradas:** `outputs/metadatos_con_suelos.tsv`, API GBIF.  
**Salidas:** `outputs/metadatos_con_gbif.tsv`.  
**Parámetros editables:** `MANGROVE_TAXA`, `RADIUS_KM`, `MAX_ROWS`, `REQUEST_PAUSE`.

