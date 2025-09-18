## 0) Estandarización de rutas y configuración

```python
from pathlib import Path
import os
import subprocess
import shlex
import shutil
import pandas as pd
import numpy as np
import requests
from urllib.parse import urlencode
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import Point
import rasterio
from rasterio.warp import transform as rio_transform

# MAGENTA_DIR: respeta env externo; si no, usa el por defecto
os.environ.setdefault("MAGENTA_DIR", "/nfs/testing/.jbalvino/MAGENTA/MAGENTA_DIR")
PROJECT_DIR = Path(os.environ.get("MAGENTA_DIR", Path.cwd())).resolve()

RAW_DIR  = PROJECT_DIR / "rawdata" / "fastq"
META_DIR = PROJECT_DIR / "metadata"
OUT_DIR  = PROJECT_DIR / "outputs"
LOGS_DIR = PROJECT_DIR / "logs"
AUX_DIR  = PROJECT_DIR / "aux"
for d in [RAW_DIR, META_DIR, OUT_DIR, LOGS_DIR, AUX_DIR]:
    d.mkdir(parents=True, exist_ok=True)

META_BASE = "metadatos_manglar"
META_FILE = OUT_DIR / f"{META_BASE}.csv"

def save_tsv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)

def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)

def ensure_float(df: pd.DataFrame, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def km_to_deg(km: float) -> float:
    return km / 111.32

print("PROJECT_DIR:", PROJECT_DIR)

# --- Entorno seguro para curl/EDirect y requests ---
CA_BUNDLE = certifi.where()
SAFE_ENV = dict(os.environ)
SAFE_ENV["CURL_CA_BUNDLE"] = CA_BUNDLE
SAFE_ENV["SSL_CERT_FILE"]  = CA_BUNDLE
SAFE_ENV.setdefault("TERM", "dumb")
SAFE_ENV.setdefault("EUTILS_TOOL", "magenta_edirect")
SAFE_ENV.setdefault("EUTILS_EMAIL", "jbalvino@masternew")
os.environ["REQUESTS_CA_BUNDLE"] = CA_BUNDLE

def run_cmd(cmd: str, retries: int = 3) -> int:
    print(">>", cmd)
    last = 1
    for i in range(1, retries+1):
        last = subprocess.call(cmd, shell=True, env=SAFE_ENV)
        if last == 0:
            return 0
        print(f"[run_cmd][retry {i}/{retries}] exit={last}")
    return last

def have_edirect() -> bool:
    for tool in ("esearch", "efetch"):
        code = subprocess.call(f"command -v {tool} >/dev/null 2>&1", shell=True, env=SAFE_ENV)
        if code != 0:
            return False
    return True

def postfilter_wgs_illumina_csv(csv_in: Path, csv_out: Path):
    """Filtra RunInfo por LibraryStrategy=WGS y Platform=ILLUMINA (case-insensitive)."""
    if not csv_in.exists() or csv_in.stat().st_size == 0:
        return
    with open(csv_in, newline='', encoding='utf-8', errors="ignore") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        if not rows:
            csv_out.write_text("") ; return
        fields = reader.fieldnames
    def ok(row):
        ls = (row.get("LibraryStrategy") or "").upper()
        pf = (row.get("Platform") or "").upper()
        return ("WGS" in ls) and ("ILLUMINA" in pf)
    kept = [r for r in rows if ok(r)]
    with open(csv_out, "w", newline='', encoding='utf-8') as g:
        w = csv.DictWriter(g, fieldnames=fields)
        w.writeheader()
        for r in kept:
            w.writerow(r)

```
### 1. Configuración y utilidades iniciales  

Importa librerías, define rutas estándar para el proyecto, crea las carpetas necesarias, configura un entorno seguro para ejecutar comandos de EDirect/requests, y declara funciones auxiliares (lectura/escritura de TSV, conversión de tipos, conversión km→grados, ejecución robusta de comandos y filtrado de metadatos WGS+ILLUMINA).  

**Entradas:**  
Ninguna (solo variables de entorno como `MAGENTA_DIR` si están definidas).  

**Salidas:**  
- Estructura de carpetas lista (`rawdata`, `metadata`, `outputs`, `logs`, `aux`)  
- Variables globales de rutas  
- Entorno seguro de red para `requests` y `curl/EDirect`  
- Funciones auxiliares disponibles (`save_tsv`, `load_tsv`, `ensure_float`, `km_to_deg`, `run_cmd`, `have_edirect`, `postfilter_wgs_illumina_csv`)  

**Parámetros editables:**  
- `MAGENTA_DIR`  
- `EUTILS_TOOL`  
- `EUTILS_EMAIL`  
- Número de reintentos en `run_cmd` 
