# 5) Preprocesamiento de lecturas: QC → trimming → ensamblaje

**Binarios requeridos en PATH:** `fastqc`, `multiqc`, `fastp`, `megahit`.

**Variables de entorno opcionales:**
- `MAG_PROJECT_DIR`: fuerza la base del proyecto (por defecto usa el `PROJECT_DIR`).
- `MAG_RAW_DIR`: ruta a FASTQ crudos. Por defecto: `PROJECT_DIR/mags/data/raw`.
- `MAG_QC_THREADS`, `MAG_FASTP_THREADS`, `MAG_FASTP_MAX_WORKERS`, `MAG_VALIDATE_READS`, `MAG_MEGAHIT_THREADS`.

```python
from __future__ import annotations
from pathlib import Path
import os, re, sys, shutil, subprocess, csv, json, gzip
from datetime import datetime
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# Raíz del proyecto: MAG_PROJECT_DIR tiene prioridad, si no MAGENTA_DIR, si no cwd.
PROJECT_DIR = Path(
    os.environ.get("MAG_PROJECT_DIR", os.environ.get("MAGENTA_DIR", Path.cwd()))
).resolve()

# Estructura base tipo MAGs dentro del proyecto
MAGS_DIR    = PROJECT_DIR / "mags"
DATA_DIR    = MAGS_DIR / "data"            # (no se usa como default de RAW)
RESULTS_DIR = MAGS_DIR / "results"
SCRIPTS_DIR = MAGS_DIR / "scripts"

# Subcarpetas de resultados
QC_DIR      = RESULTS_DIR / "01.qc"
TRIM_DIR    = RESULTS_DIR / "02.trimmed"
ASM_DIR     = RESULTS_DIR / "03.assembly"
ASM_LOG_DIR = ASM_DIR / "logs"

# Metadatos / reportes
PIPE_META          = RESULTS_DIR / "pipeline_meta"
PIPE_META.mkdir(parents=True, exist_ok=True)
PIPE_MANIFEST_CSV  = PIPE_META / "pipeline_manifest.csv"
TRIM_REPORT_CSV    = PIPE_META / "trim_report.csv"
ASM_RESUMEN_CSV    = ASM_DIR   / "resumen_ensamblaje.csv"
VALID_SAMPLES_CSV  = PIPE_META / "muestras_validas.csv"

# Donde 02_descargar_y_convertir_mangrove.py deja los FASTQ
DEFAULT_RAW_SRC = PROJECT_DIR / "rawdata" / "convertidos"
RAW_SRC = Path(os.environ.get("MAG_RAW_DIR", DEFAULT_RAW_SRC))

# Crear estructura necesaria
for d in [MAGS_DIR, DATA_DIR, RESULTS_DIR, SCRIPTS_DIR, QC_DIR, TRIM_DIR, ASM_DIR, ASM_LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Parámetros
N_THREADS_FASTQC   = int(os.environ.get("MAG_QC_THREADS",        "12"))
N_THREADS_FASTP    = int(os.environ.get("MAG_FASTP_THREADS",     "12"))
MAX_WORKERS_FASTP  = int(os.environ.get("MAG_FASTP_MAX_WORKERS", "4"))
MAX_READS_VALIDATE = int(os.environ.get("MAG_VALIDATE_READS",    "10000"))
MEGAHIT_THREADS    = int(os.environ.get("MAG_MEGAHIT_THREADS",   "40"))
SAMPLE_REGEX       = os.environ.get(
    "MAG_SAMPLE_REGEX",
    r"^([A-Za-z0-9_\-]+)_([12])\.fastq(?:\.gz)?$"
)

def have_bin(b: str) -> bool:
    return shutil.which(b) is not None

def run(cmd, **kwargs) -> int:
    printable = " ".join(map(str, cmd)) if isinstance(cmd, (list,tuple)) else str(cmd)
    print("[CMD]", printable)
    return subprocess.call(cmd, shell=isinstance(cmd, str), **kwargs)

print("[PORTABLE] PROJECT_DIR =", PROJECT_DIR)
print("[PORTABLE] RAW_SRC     =", RAW_SRC)
print("[PORTABLE] RESULTS_DIR  =", RESULTS_DIR)
```

## Rutas y estructura de directorios

**Base del proyecto (PROJECT_DIR):**  
Determinada en este orden:
1. `MAG_PROJECT_DIR`
2. `MAGENTA_DIR`
3. Directorio de trabajo actual (`cwd`)

**Estructura creada dentro de `PROJECT_DIR`:**
- `mags/`
- `data/` (no se usa como entrada por defecto de FASTQ)
- `results/`
  - `01.qc/` (FastQC crudo y post-trim, MultiQC)
  - `02.trimmed/` (salidas de fastp)
  - `03.assembly/` (ensamblajes MEGAHIT)
- `pipeline_meta/` (manifiestos y reportes)
- `scripts/`

**Entrada por defecto de FASTQ:**  
`rawdata/convertidos/`, generada por el script de descarga y conversión.  
Puede sobrescribirse con la variable `MAG_RAW_DIR`.

---

## Parámetros configurables

Los siguientes parámetros se leen de variables de entorno, con valores por defecto sensatos:

- `MAG_QC_THREADS` (12) → hilos para FastQC  
- `MAG_FASTP_THREADS` (12) → hilos para fastp  
- `MAG_FASTP_MAX_WORKERS` (4) → número máximo de muestras en paralelo para fastp  
- `MAG_VALIDATE_READS` (10000) → lecturas muestreadas para validación de FASTQ  
- `MAG_MEGAHIT_THREADS` (40) → hilos para MEGAHIT  
- `MAG_SAMPLE_REGEX` → expresión regular para identificar muestras a partir del nombre de archivo  

---

## Utilidades del bloque

- **`have_bin(b: str)`**  
  Comprueba si un binario requerido está disponible en el `PATH`.

- **`run(cmd, **kwargs)`**  
  Ejecuta un comando, imprime la línea para trazabilidad y devuelve el código de salida.

---

## Archivos de salida esperados

- `pipeline_manifest.csv` → listado de muestras y rutas de FASTQ detectadas  
- `trim_report.csv` → reporte de trimming  
- `resumen_ensamblaje.csv` → estado de los ensamblajes por muestra  
- `muestras_validas.csv` → resumen de pares de FASTQ validados  

---

## 🔍 5.1 Control de calidad inicial — FastQC + MultiQC
```python
import glob

# 1) Detectar pares (y single-end) a partir de SAMPLE_REGEX
def list_pairs(raw_dir: Path, sample_regex: str) -> dict[str, dict[str, Path]]:
    rx = re.compile(sample_regex)
    buckets: dict[str, dict[str, Path]] = {}
    for p in sorted(raw_dir.rglob("*.fastq*")):
        m = rx.match(p.name)
        if not m:
            continue
        sample, read = m.group(1), m.group(2)
        d = buckets.setdefault(sample, {"R1": None, "R2": None})
        if read == "1" and d["R1"] is None:
            d["R1"] = p
        elif read == "2" and d["R2"] is None:
            d["R2"] = p
    # Aceptar single-end con nombres <sample>.fastq(.gz)
    for p in sorted(raw_dir.rglob("*.fastq*")):
        if rx.match(p.name):
            continue
        base = p.name.replace(".gz", "")
        if base.endswith(".fastq"):
            sample = base[:-6]  # quita ".fastq"
            buckets.setdefault(sample, {"R1": None, "R2": None})
            if buckets[sample]["R1"] is None:
                buckets[sample]["R1"] = p
    return buckets

def write_manifest(pairs: dict[str, dict[str, Path]]) -> pd.DataFrame:
    rows = []
    for s, rr in sorted(pairs.items()):
        rows.append({"sample": s, "R1": str(rr.get("R1") or ""), "R2": str(rr.get("R2") or "")})
    df = pd.DataFrame(rows)
    df.to_csv(PIPE_MANIFEST_CSV, index=False)
    print(f"[META] Manifiesto muestras -> {PIPE_MANIFEST_CSV}")
    return df

# 2) FastQC crudo (paralelo por lotes para no saturar la CLI)
def run_fastqc(inputs: list[Path], outdir: Path, threads: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if not inputs:
        print("[FASTQC] Sin archivos de entrada.")
        return
    batch = []
    for p in inputs:
        batch.append(str(p))
        if len(batch) >= 64:
            run(["fastqc", "-t", str(threads), "-o", str(outdir)] + batch)
            batch = []
    if batch:
        run(["fastqc", "-t", str(threads), "-o", str(outdir)] + batch)

# Ejecutar manifiesto + FastQC
if not RAW_SRC.exists():
    raise SystemExit(f"[ERROR] No existe RAW_SRC: {RAW_SRC}")

pairs = list_pairs(RAW_SRC, SAMPLE_REGEX)
if not pairs:
    raise SystemExit(f"[ERROR] No se detectaron FASTQ en {RAW_SRC}")

manifest_df = write_manifest(pairs)

raw_fastqs = []
for rr in pairs.values():
    if rr.get("R1"): raw_fastqs.append(rr["R1"])
    if rr.get("R2"): raw_fastqs.append(rr["R2"])

if have_bin("fastqc"):
    run_fastqc(raw_fastqs, QC_DIR / "raw_fastqc", N_THREADS_FASTQC)
    if have_bin("multiqc"):
        run(["multiqc", str(QC_DIR / "raw_fastqc"), "-o", str(QC_DIR / "raw_fastqc")])
else:
    print("[AVISO] 'fastqc' no está en PATH. Se omite QC inicial.")

  ```
Este bloque se encarga de **detectar archivos FASTQ** en el directorio de entrada (`RAW_SRC`), generar un **manifiesto de muestras** y ejecutar **FastQC en paralelo** (con soporte opcional para MultiQC).

## Flujo principal

1. **Descubrimiento de muestras (`list_pairs`)**
   - Detecta pares de archivos `R1`/`R2` usando `SAMPLE_REGEX`.
   - Admite también **single-end** con nombres `<sample>.fastq(.gz)`.
   - Devuelve un diccionario con las rutas asociadas a cada muestra.

2. **Generación de manifiesto (`write_manifest`)**
   - Construye un `DataFrame` con columnas:
     - `sample`
     - `R1`
     - `R2`
   - Guarda el archivo `pipeline_manifest.csv`.
   - Imprime la ruta generada.

3. **Ejecución de FastQC (`run_fastqc`)**
   - Revisa que existan archivos de entrada.
   - Crea el directorio `QC_DIR/raw_fastqc` si no existe.
   - Procesa los archivos en **lotes de hasta 64** para no saturar la CLI.
   - Ejecuta:  
     ```bash
     fastqc -t <threads> -o <QC_DIR/raw_fastqc> <archivos>
     ```
     ### Parámetros en *paired-end*

Cuando hay **R1** y **R2**:

- `-I`: archivo de entrada **R2** (lecturas *reverse*).  
- `-O`: archivo de salida de lecturas **R2** procesadas.  
- `--detect_adapter_for_pe`: activa la detección automática de adaptadores para datos *paired-end*.  

4. **MultiQC opcional**
   - Si está disponible, ejecuta:  
     ```bash
     multiqc <QC_DIR/raw_fastqc> -o <QC_DIR/raw_fastqc>
     ```
---

## Condiciones y validaciones

- Si `RAW_SRC` no existe → termina con error.  
- Si no se detectan FASTQ → termina con error.  
- Si `fastqc` no está en el `PATH` → omite la etapa con aviso.  
- Si `multiqc` no está disponible → solo se ejecuta FastQC.  

## Archivos generados

- `pipeline_manifest.csv` → listado de muestras y archivos FASTQ detectados.  
- Resultados de **FastQC** en `QC_DIR/raw_fastqc/`.  
- Reporte combinado de **MultiQC** (si está instalado).

## 5.2 Trimming de lecturas — *fastp* en paralelo

```python
from __future__ import annotations
from pathlib import Path
import os, re, sys, shutil, subprocess, argparse, shlex
import pandas as pd
from typing import Dict, Optional, List, Tuple

# ---------------------------
# 0) Rutas y parámetros base
# ---------------------------
PROJECT_DIR = Path(os.environ.get("MAG_PROJECT_DIR", os.environ.get("MAGENTA_DIR", Path.cwd()))).resolve()
MAGS_DIR    = PROJECT_DIR / "mags"
RESULTS_DIR = MAGS_DIR / "results"
TRIM_DIR    = RESULTS_DIR / "02.trimmed"
PIPE_META   = RESULTS_DIR / "pipeline_meta"
PIPE_META.mkdir(parents=True, exist_ok=True)

PIPE_MANIFEST_CSV = PIPE_META / "pipeline_manifest.csv"
TRIM_REPORT_CSV   = PIPE_META / "trim_report.csv"

DEFAULT_RAW_SRC = PROJECT_DIR / "rawdata" / "convertidos"
RAW_SRC = Path(os.environ.get("MAG_RAW_DIR", DEFAULT_RAW_SRC))

SAMPLE_REGEX = os.environ.get("MAG_SAMPLE_REGEX", r"^([A-Za-z0-9_\-]+)_([12])\.fastq(?:\.gz)?$")

N_THREADS_FASTP = int(os.environ.get("MAG_FASTP_THREADS", "12"))

def have_bin(b: str) -> bool:
    from shutil import which
    return which(b) is not None

def run(cmd, **kwargs) -> int:
    printable = " ".join(map(str, cmd)) if isinstance(cmd, (list,tuple)) else str(cmd)
    print("[CMD]", printable, flush=True)
    return subprocess.call(cmd, shell=isinstance(cmd, str), **kwargs)

print(f"[TRIM] PROJECT_DIR = {PROJECT_DIR}")
print(f"[TRIM] RAW_SRC     = {RAW_SRC}")
print(f"[TRIM] TRIM_DIR    = {TRIM_DIR}")

# --------------------------------
# 1) Cargar o descubrir 'pairs'
# --------------------------------
def rebuild_pairs_from_manifest(manifest_csv: Path) -> Dict[str, Dict[str, Optional[Path]]]:
    df = pd.read_csv(manifest_csv)
    pairs: Dict[str, Dict[str, Optional[Path]]] = {}
    for _, row in df.iterrows():
        s = str(row["sample"])
        r1 = Path(row["R1"]) if isinstance(row["R1"], str) and row["R1"] else None
        r2 = Path(row["R2"]) if isinstance(row["R2"], str) and row["R2"] else None
        pairs[s] = {"R1": r1, "R2": r2}
    return pairs

def discover_pairs(raw_dir: Path, sample_regex: str) -> Dict[str, Dict[str, Optional[Path]]]:
    rx = re.compile(sample_regex)
    buckets: Dict[str, Dict[str, Optional[Path]]] = {}
    for p in sorted(raw_dir.rglob("*.fastq*")):
        m = rx.match(p.name)
        if not m: continue
        sample, read = m.group(1), m.group(2)
        d = buckets.setdefault(sample, {"R1": None, "R2": None})
        if read == "1" and d["R1"] is None: d["R1"] = p
        elif read == "2" and d["R2"] is None: d["R2"] = p
    # Single-end: <sample>.fastq(.gz)
    for p in sorted(raw_dir.rglob("*.fastq*")):
        if rx.match(p.name): continue
        base = p.name.replace(".gz", "")
        if base.endswith(".fastq"):
            sample = base[:-6]
            buckets.setdefault(sample, {"R1": None, "R2": None})
            if buckets[sample]["R1"] is None:
                buckets[sample]["R1"] = p
    return buckets

# --------------------------------
# 2) fastp (default, sin FastQC/MultiQC)
# --------------------------------
def build_fastp_cmd(R1: Path, R2: Optional[Path], r1_out: Path, r2_out: Optional[Path],
                    threads: int, extra_opts: List[str]) -> List[str]:
    cmd = ["fastp", "-w", str(threads), "-o", str(r1_out)]
    if R2:
        cmd += ["-O", str(r2_out), "-i", str(R1), "-I", str(R2), "--detect_adapter_for_pe"]
    else:
        cmd += ["-i", str(R1)]
    # NOTA: NO añadimos -j/-h para no generar reportes fastp explícitos (evitar ruido).
    # fastp aún podría crear fastp.json/html por defecto en algunos builds; se ignoran.
    cmd += extra_opts
    return cmd

def fastp_one(sample: str, R1: Optional[Path], R2: Optional[Path],
              threads: int, extra_opts: List[str]) -> dict:
    TRIM_DIR.mkdir(parents=True, exist_ok=True)
    if not have_bin("fastp"):
        return {"sample": sample, "ok": False, "reason": "fastp_not_found"}
    if not (R1 and Path(R1).exists()):
        return {"sample": sample, "ok": False, "reason": "no_R1"}
    if R2 and not Path(R2).exists():
        return {"sample": sample, "ok": False, "reason": "R2_missing"}

    r1_out = TRIM_DIR / f"{sample}_R1.fastq.gz"
    r2_out = TRIM_DIR / f"{sample}_R2.fastq.gz" if R2 else None

    code = run(build_fastp_cmd(R1, R2, r1_out, r2_out, threads, extra_opts))
    ok = (code == 0) and r1_out.exists() and (r2_out is None or r2_out.exists())
    return {
        "sample": sample,
        "ok": ok,
        "trimmer": "fastp",
        "R1_out": str(r1_out) if r1_out.exists() else "",
        "R2_out": str(r2_out) if (r2_out and r2_out.exists()) else "",
    }

# --------------------------------
# 3) Trim Galore (sin --fastqc)
# --------------------------------
def _find_tg_outputs_pe(R1: Path, R2: Path, outdir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    base1 = Path(R1).stem
    base2 = Path(R2).stem
    cand1 = list(outdir.glob(f"*{base1}*_val_1.fq.gz")) + list(outdir.glob(f"*{base1}*_val_1.fastq.gz"))
    cand2 = list(outdir.glob(f"*{base2}*_val_2.fq.gz")) + list(outdir.glob(f"*{base2}*_val_2.fastq.gz"))
    if not cand1:
        cand1 = list(outdir.glob("*_val_1.fq.gz")) + list(outdir.glob("*_val_1.fastq.gz"))
    if not cand2:
        cand2 = list(outdir.glob("*_val_2.fq.gz")) + list(outdir.glob("*_val_2.fastq.gz"))
    return (cand1[0] if cand1 else None, cand2[0] if cand2 else None)

def _find_tg_outputs_se(R1: Path, outdir: Path) -> Optional[Path]:
    base1 = Path(R1).stem
    cand = list(outdir.glob(f"*{base1}*trimmed.fq.gz")) + list(outdir.glob(f"*{base1}*trimmed.fastq.gz"))
    if not cand:
        cand = list(outdir.glob("*trimmed.fq.gz")) + list(outdir.glob("*trimmed.fastq.gz"))
    return cand[0] if cand else None

def trim_galore_one(sample: str, R1: Optional[Path], R2: Optional[Path],
                    tg_opts_pe: List[str], tg_opts_se: List[str]) -> dict:
    TRIM_DIR.mkdir(parents=True, exist_ok=True)
    if not have_bin("trim_galore"):
        return {"sample": sample, "ok": False, "reason": "trim_galore_not_found"}
    if not (R1 and Path(R1).exists()):
        return {"sample": sample, "ok": False, "reason": "no_R1"}
    if R2 and not Path(R2).exists():
        return {"sample": sample, "ok": False, "reason": "R2_missing"}

    # Importante: sin --fastqc y con --no_report_file para evitar QC redundante
    base_opts = ["--no_report_file"]
    if R2:
        cmd = ["trim_galore", "--paired"] + base_opts + tg_opts_pe + ["-o", str(TRIM_DIR), str(R1), str(R2)]
    else:
        cmd = ["trim_galore"] + base_opts + tg_opts_se + ["-o", str(TRIM_DIR), str(R1)]
    code = run(cmd)

    if R2:
        out1_src, out2_src = _find_tg_outputs_pe(R1, R2, TRIM_DIR)
        out1_dst = TRIM_DIR / f"{sample}_R1.fastq.gz"
        out2_dst = TRIM_DIR / f"{sample}_R2.fastq.gz"
        if out1_src and not out1_dst.exists(): shutil.move(str(out1_src), str(out1_dst))
        if out2_src and not out2_dst.exists(): shutil.move(str(out2_src), str(out2_dst))
        ok = (code == 0) and out1_dst.exists() and out2_dst.exists()
        return {"sample": sample, "ok": ok, "trimmer": "trim_galore",
                "R1_out": str(out1_dst) if out1_dst.exists() else "",
                "R2_out": str(out2_dst) if out2_dst.exists() else ""}
    else:
        out_src = _find_tg_outputs_se(R1, TRIM_DIR)
        out_dst = TRIM_DIR / f"{sample}_R1.fastq.gz"
        if out_src and not out_dst.exists(): shutil.move(str(out_src), str(out_dst))
        ok = (code == 0) and out_dst.exists()
        return {"sample": sample, "ok": ok, "trimmer": "trim_galore",
                "R1_out": str(out_dst) if out_dst.exists() else "",
                "R2_out": ""}

# --------------------------------
# 4) Main
# --------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Trimming con fastp (default) o Trim Galore, sin QC post-trimming.")
    parser.add_argument("--trimmer", choices=["fastp","trim_galore"], default=os.environ.get("MAG_TRIM_TOOL","fastp"),
                        help="Herramienta de trimming (default: fastp).")
    parser.add_argument("--threads-fastp", type=int, default=N_THREADS_FASTP, help="Hilos para fastp.")
    parser.add_argument("--fastp-opts", type=str, default=os.environ.get("MAG_FASTP_OPTS",""),
                        help="Flags extra para fastp (string estilo shell).")
    parser.add_argument("--tg-pe-opts", type=str,
                        default=os.environ.get("MAG_TRIMGALORE_PE","--cores 20 --quality 15 --length 50"),
                        help="Flags Trim Galore para paired-end (sin --fastqc).")
    parser.add_argument("--tg-se-opts", type=str,
                        default=os.environ.get("MAG_TRIMGALORE_SE","--cores 20 --quality 15 --length 50"),
                        help="Flags Trim Galore para single-end (sin --fastqc).")
    args = parser.parse_args(argv)

    TRIM_DIR.mkdir(parents=True, exist_ok=True)

    # Pairs desde manifiesto o descubrimiento
    if PIPE_MANIFEST_CSV.exists():
        pairs = rebuild_pairs_from_manifest(PIPE_MANIFEST_CSV)
        print(f"[TRIM] Cargado pairs desde {PIPE_MANIFEST_CSV} ({len(pairs)} muestras).")
    else:
        if not RAW_SRC.exists():
            print(f"[ERROR] No existe RAW_SRC: {RAW_SRC}", file=sys.stderr)
            return 2
        pairs = discover_pairs(RAW_SRC, SAMPLE_REGEX)
        print(f"[TRIM] Descubierto pairs escaneando {RAW_SRC} ({len(pairs)} muestras).")
    if not pairs:
        print("[ERROR] No hay muestras detectadas.", file=sys.stderr)
        return 3

    reports: List[dict] = []
    for sample, rr in sorted(pairs.items()):
        R1, R2 = rr.get("R1"), rr.get("R2")
        if args.trimmer == "fastp":
            extra = shlex.split(args.fastp_opts) if args.fastp_opts else []
            rep = fastp_one(sample, R1, R2, threads=args.threads_fastp, extra_opts=extra)
        else:
            tg_pe = shlex.split(args.tg_pe_opts) if args.tg_pe_opts else []
            tg_se = shlex.split(args.tg_se_opts) if args.tg_se_opts else []
            rep = trim_galore_one(sample, R1, R2, tg_pe, tg_se)
        reports.append(rep)

    pd.DataFrame(reports).sort_values("sample").to_csv(TRIM_REPORT_CSV, index=False)
    print(f"[TRIM] Reporte -> {TRIM_REPORT_CSV}")
    print("[TRIM] Finalizado (sin QC post-trim: usar 04_fastqc_parallel.py).")
    return 0

if __name__ == "__main__":
    sys.exit(main())

 ```
Implementa un módulo de **trimming de lecturas** para un pipeline de ensamblaje de MAGs.  
Es compatible con dos herramientas de recorte de secuencias:

- **fastp** (por defecto)  
- **Trim Galore**

Su propósito es limpiar las lecturas crudas, generar archivos de salida filtrados y producir un reporte de trimming.

1. Configura rutas de proyecto, resultados y directorios de calidad (`QC_DIR`, `TRIM_DIR`, `PIPE_META`).  
2. Carga un **manifiesto de muestras** (`pipeline_manifest.csv`) con la información de las lecturas.  
3. Ejecuta el trimming con la herramienta seleccionada:
   - **fastp**: usando hilos definidos (`N_THREADS_FASTP`) y soporte para *single-end* o *paired-end*.  
   - **Trim Galore**: ajustando automáticamente según el tipo de datos (single o paired).  
4. Guarda un reporte (`trim_report.csv`) con el detalle de los archivos generados por cada muestra.  

## Inputs
- **Archivo de manifiesto** (`pipeline_manifest.csv`), con columnas:  
  - `sample`: nombre de la muestra  
  - `R1`: archivo de lecturas forward  
  - `R2`: archivo de lecturas reverse (opcional, si es paired-end)  
- Variables de entorno (opcionales):  
  - `MAG_PROJECT_DIR` o `MAGENTA_DIR`: define el directorio base del proyecto.  
  - `MAG_FASTP_THREADS`: número de hilos para `fastp` (default: 12).  
  - `MAG_FASTP_MAX_WORKERS`: número máximo de procesos paralelos (default: 4).  
  - `MAG_TRIM_TOOL`: herramienta de trimming (`fastp` o `trim_galore`).  

## Outputs
- Archivos **FASTQ recortados**, ubicados en `02.trimmed/<sample>/`:  
  - `<sample>_R1.trimmed.fastq.gz`  
  - `<sample>_R2.trimmed.fastq.gz` (si aplica)  
- **Reporte de trimming** (`trim_report.csv`) con columnas:  
  - `sample` → nombre de muestra  
  - `tool` → herramienta usada (`fastp` o `trim_galore`)  
  - `R1` → archivo trimmed forward  
  - `R2` → archivo trimmed reverse (si aplica)  

## Parámetros usados

### Trim Galore
Parámetros mínimos aplicados en el flujo:

- `--paired` → para lecturas pareadas  
- `--no_report_file` → desactiva reportes por muestra (se usa MultiQC después)  
- `--cores 20` → número de hilos  
- `--quality 15` → trimming de bases con calidad < Q15  
- `--length 50` → descarta lecturas post-trimming con longitud < 50 nt  
- `--fastqc` → corre FastQC por muestra  

### fastp (por defecto)

- `-q 15` / `--qualified_quality_phred 15` → calidad mínima por base (Q ≥ 15)  
- `-u 40` / `--unqualified_percent_limit 40` → descarta la lectura si >40% de bases < Q15  
- `-l 15` / `--length_required 15` → descarta lecturas más cortas que 15 nt  
- `--detect_adapter_for_pe` → detección y recorte automático de adaptadores (paired-end)  
- **Poly-G trimming** → activado automáticamente en datos de NovaSeq/NextSeq  

- **Trim Galore** en el pipeline está configurado para ser más estricto con la longitud (≥50 nt) y siempre corre FastQC.  
- **fastp**, por defecto, aplica filtros de calidad equivalentes (Q15) y longitud mínima (≥15 nt), además de recorte automático de adaptadores.

## 5.3 Ensamblaje metagenómico (MEGAHIT por defecto, MetaSPAdes opcional)

```python
from __future__ import annotations
import os, sys, csv, shutil, subprocess, argparse
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Optional

# ==========
# Rutas base
# ==========
PROJECT_DIR = Path(os.environ.get("MAG_PROJECT_DIR", os.environ.get("MAGENTA_DIR", Path.cwd()))).resolve()
MAGS_DIR    = PROJECT_DIR / "mags"
RESULTS_DIR = MAGS_DIR / "results"

TRIM_DIR    = RESULTS_DIR / "02.trimmed"
ASM_DIR     = RESULTS_DIR / "03.assembly"
ASM_LOG_DIR = ASM_DIR / "logs"
RESUMEN_CSV = ASM_DIR / "resumen_ensamblaje.csv"

ASM_LOG_DIR.mkdir(parents=True, exist_ok=True)

# ==========
# Utilidades
# ==========
def have_bin(b: str) -> bool:
    from shutil import which
    return which(b) is not None

def run(cmd, log_file: Path | None = None) -> int:
    printable = " ".join(map(str, cmd))
    print("[CMD]", printable, flush=True)
    if log_file:
        with open(log_file, "a") as log:
            log.write(f"[CMD] {printable}\n")
            return subprocess.call(cmd, stdout=log, stderr=log)
    else:
        return subprocess.call(cmd)

def size_sum(f1: Path, f2: Path) -> int:
    try:
        return f1.stat().st_size + f2.stat().st_size
    except Exception:
        return 0

# =========================
# Descubrimiento de pares
# =========================
R1_GLOBS_DEFAULT = [
    "*_R1*.fastq.gz", "*_1*.fastq.gz", "*R1*.fastq.gz", "*1*.fastq.gz",
    "*_R1*.fq.gz",    "*_1*.fq.gz",    "*R1*.fq.gz",    "*1*.fq.gz",
]
TOKEN_PAIRS = [("_R1", "_R2"), ("_1", "_2"), ("R1", "R2"), ("1", "2")]

def _try_match_r2_by_tokens(r1: Path) -> Optional[Path]:
    """Intenta encontrar R2 reemplazando tokens comunes en el nombre de R1, en el mismo directorio."""
    name = r1.name
    for t1, t2 in TOKEN_PAIRS:
        if t1 in name:
            cand = r1.with_name(name.replace(t1, t2, 1))
            if cand.exists():
                return cand
    return None

def _best_effort_r2_in_dir(r1: Path) -> Optional[Path]:
    """
    Si el reemplazo directo falla, intenta localizar un R2 plausible en el mismo directorio:
    - Algún archivo que contenga un token R2 y termine en .f*q.gz
    - Preferir el que comparte mayor prefijo con R1.
    """
    dirp = r1.parent
    candidates: List[Path] = []
    for patt in ["*_R2*.fastq.gz", "*_2*.fastq.gz", "*R2*.fastq.gz", "*2*.fastq.gz",
                 "*_R2*.fq.gz",   "*_2*.fq.gz",   "*R2*.fq.gz",   "*2*.fq.gz"]:
        candidates.extend(dirp.glob(patt))
    if not candidates:
        return None

    def common_prefix_len(a: str, b: str) -> int:
        n = min(len(a), len(b)); i = 0
        while i < n and a[i] == b[i]:
            i += 1
        return i

    best = max(candidates, key=lambda p: common_prefix_len(p.name, r1.name))
    return best

def _collect_sample_dirs(root: Path) -> List[Path]:
    """Devuelve subdirectorios de primer nivel si existen; si no, devuelve [root] para fallback plano."""
    subdirs = [d for d in sorted(root.iterdir()) if d.is_dir()]
    return subdirs if subdirs else [root]

def find_pairs(trim_root: Path, debug: bool=False) -> List[Tuple[str, Path, Path]]:
    """
    Busca pares R1/R2 en estructura:
      trim_root/<sample>/**/{R1,R2}*.f*q.gz
    - Explora recursivamente debajo de cada subcarpeta de muestra.
    - Acepta múltiples patrones de nombres.
    Devuelve: lista de (sample, R1, R2) con sample = nombre del subdirectorio de primer nivel
              o, en fallback plano, el prefijo derivado del archivo R1.
    """
    pairs: List[Tuple[str, Path, Path]] = []
    sample_dirs = _collect_sample_dirs(trim_root)

    for sample_dir in sample_dirs:
        if sample_dir == trim_root and sample_dir.is_dir() and not any(sample_dir.iterdir()):
            continue

        sample_name = sample_dir.name if sample_dir != trim_root else None

        r1_candidates: List[Path] = []
        for patt in R1_GLOBS_DEFAULT:
            r1_candidates.extend(sample_dir.rglob(patt))
        r1_candidates = sorted(set(r1_candidates))

        if debug:
            print(f"[DEBUG] Dir muestra: {sample_dir}")
            for r1 in r1_candidates:
                print(f"[DEBUG]   R1 cand: {r1}")

        used: set[Path] = set()
        for r1 in r1_candidates:
            if r1 in used:
                continue
            r2 = _try_match_r2_by_tokens(r1)
            if r2 is None or not r2.exists():
                r2 = _best_effort_r2_in_dir(r1)
            if r2 and r2.exists() and r2 not in used:
                used.add(r1); used.add(r2)
                if sample_dir == trim_root:
                    sname = r1.name
                    for suf in ["_R1.fastq.gz","_1.fastq.gz","R1.fastq.gz","1.fastq.gz",
                                "_R1.fq.gz","_1.fq.gz","R1.fq.gz","1.fq.gz"]:
                        if sname.endswith(suf):
                            sname = sname[:-len(suf)]
                            break
                else:
                    sname = sample_name
                pairs.append((sname, r1, r2))

        if debug and not pairs:
            print(f"[DEBUG]   No se encontraron pares en {sample_dir}")

    return pairs

# =========================
# Normalización de salidas
# =========================
def normalize_and_prefix_outputs_megahit(sample: str, out_dir: Path, log: Path):
    """
    MEGAHIT produce final.contigs.fa.
    Crear contigs.fasta y renombrar a <sample>_contigs.fasta para uniformidad.
    """
    src = out_dir / "final.contigs.fa"
    mid = out_dir / "contigs.fasta"
    if src.exists():
        try:
            shutil.copyfile(src, mid)
        except Exception as e:
            with open(log, "a") as lg:
                lg.write(f"[WARN] No se pudo copiar a contigs.fasta: {e}\n")
    for fname in ["contigs.fasta", "scaffolds.fasta"]:
        source = out_dir / fname
        if source.exists():
            renamed = out_dir / f"{sample}_{fname}"
            if not renamed.exists():
                source.rename(renamed)
            print(f"[RENAME] {renamed.name}")

def prefix_outputs_metaspades(sample: str, out_dir: Path):
    """
    MetaSPAdes produce contigs.fasta y scaffolds.fasta.
    Renombrar con prefijo de muestra.
    """
    for fname in ["contigs.fasta", "scaffolds.fasta"]:
        source = out_dir / fname
        if source.exists():
            renamed = out_dir / f"{sample}_{fname}"
            if not renamed.exists():
                source.rename(renamed)
            print(f"[RENAME] {renamed.name}")

# ==========
# Ensamblaje
# ==========
def assemble_one(sample: str, R1: Path, R2: Path, assembler: str, threads: int, mem_gb: int,
                 overwrite: bool=False) -> tuple[str, str, str]:
    sample_outdir = ASM_DIR / f"{sample}_assembled"
    log_file = ASM_LOG_DIR / f"{sample}.log"

    # Manejo del directorio de salida
    if sample_outdir.exists():
        if overwrite:
            try:
                shutil.rmtree(sample_outdir)
            except Exception as e:
                return (sample, "ERROR", f"No se pudo borrar salida previa: {e}")
        else:
            return (sample, "SKIPPED", "Ya existe (use --overwrite para rehacer)")

    try:
        if assembler == "megahit":
            if not have_bin("megahit"):
                return (sample, "ERROR", "megahit no encontrado en PATH")
            cmd = [
                "megahit",
                "-1", str(R1),
                "-2", str(R2),
                "-o", str(sample_outdir),   # dejar que megahit cree la carpeta
                "-t", str(threads),
                "--presets", "meta-sensitive",
                # Correcciones añadidas (flags extra de ejemplo):
                "--min-contig-len", "1000",
                "--k-list", "21,29,39,59,79,99,119",
            ]
            code = run(cmd, log_file)
            if code != 0:
                try:
                    print(f"[LOG TAIL] {log_file}")
                    os.system(f"tail -n 200 {log_file}")
                except Exception:
                    pass
                return (sample, "ERROR", f"megahit exit code {code}")
            normalize_and_prefix_outputs_megahit(sample, sample_outdir, log_file)

        elif assembler == "metaspades":
            exe = "metaspades.py" if have_bin("metaspades.py") else ("spades.py" if have_bin("spades.py") else None)
            if exe is None:
                return (sample, "ERROR", "metaspades.py no encontrado en PATH")
            cmd = [exe, "-1", str(R1), "-2", str(R2), "-o", str(sample_outdir), "-t", str(threads)]
            if mem_gb and int(mem_gb) > 0:
                cmd += ["-m", str(int(mem_gb))]
            # Correcciones añadidas (flag extra de ejemplo):
            cmd += ["--only-assembler"]  # si ya tienes corrección hecha
            # cmd += ["--meta"]  # (metaspades.py ya implica modo meta)
            code = run(cmd, log_file)
            if code != 0:
                try:
                    print(f"[LOG TAIL] {log_file}")
                    os.system(f"tail -n 200 {log_file}")
                except Exception:
                    pass
                return (sample, "ERROR", f"metaspades exit code {code}")
            prefix_outputs_metaspades(sample, sample_outdir)

        else:
            return (sample, "ERROR", f"Assembler no soportado: {assembler}")

        return (sample, "OK", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    except Exception as e:
        try:
            print(f"[LOG TAIL] {log_file}")
            os.system(f"tail -n 200 {log_file}")
        except Exception:
            pass
        return (sample, "ERROR", str(e))

def obtener_muestras_pendientes(pairs: List[Tuple[str, Path, Path]]) -> List[str]:
    """
    Filtra muestras sin carpeta de ensamblaje y las ordena por tamaño total ascendente.
    """
    info = []
    for sample, r1, r2 in pairs:
        outfolder = ASM_DIR / f"{sample}_assembled"
        if not outfolder.exists():
            sz = size_sum(r1, r2)
            info.append((sample, sz))
    return [m for m, _ in sorted(info, key=lambda x: x[1])]

# ==========
# Main
# ==========
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Ensamblaje metagenómico en serie (MEGAHIT/MetaSPAdes) con 02.trimmed/<sample>/**/R1,R2.")
    parser.add_argument("--assembler", choices=["megahit","metaspades"],
                        default=os.environ.get("MAG_ASSEMBLER", "megahit"),
                        help="Algoritmo de ensamblaje (default: megahit).")
    parser.add_argument("--threads", type=int,
                        default=int(os.environ.get("MAG_MEGAHIT_THREADS", os.environ.get("MAG_SPADES_THREADS", "40"))),
                        help="Hilos por muestra (default: env o 40).")
    parser.add_argument("--mem-gb", type=int,
                        default=int(os.environ.get("MAG_SPADES_MEM_GB", "0")),
                        help="Memoria (GB) para MetaSPAdes; 0=auto.")
    parser.add_argument("--debug", action="store_true",
                        help="Imprime información de depuración del descubrimiento de pares.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Si existe el directorio de salida de la muestra, lo elimina y vuelve a ensamblar.")
    args = parser.parse_args(argv)

    if not TRIM_DIR.exists():
        print(f"[ERROR] No existe carpeta de trimmed: {TRIM_DIR}", file=sys.stderr)
        return 2

    pairs = find_pairs(TRIM_DIR, debug=args.debug)
    if not pairs:
        print(f"[ERROR] No se detectaron pares R1/R2 dentro de subdirectorios en {TRIM_DIR}", file=sys.stderr)
        if not args.debug:
            print("Sugerencia: ejecute nuevamente con --debug para imprimir candidatos encontrados.", file=sys.stderr)
        return 3

    pendientes = obtener_muestras_pendientes(pairs)
    print(f" Muestras detectadas: {len(pairs)}")
    print(f" Muestras a ensamblar (ordenadas por tamaño): {len(pendientes)}")
    print(f" Iniciando ensamblaje en serie con {args.threads} hilos por muestra usando '{args.assembler}'...\n")

    idx = {s: (r1, r2) for s, r1, r2 in pairs}

    resultados = []
    for sample in pendientes:
        R1, R2 = idx[sample]
        print(f" Ensamblando: {sample}")
        res = assemble_one(sample, R1, R2, args.assembler, args.threads, args.mem_gb, overwrite=args.overwrite)
        resultados.append(res)
        print(f" Resultado: {res[1]} - {res[2]}\n")

    ASM_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESUMEN_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Muestra", "Estado", "Detalle"])
        writer.writerows(resultados)

    print(f"\n Logs por muestra en: {ASM_LOG_DIR}")
    print(f" Contigs generados en: {ASM_DIR}")
    print(f" Resumen CSV: {RESUMEN_CSV}")
    print("\n Ensamblajes completados.")
    return 0

if __name__ == "__main__":
    sys.exit(main())


```
# Ensamblaje metagenómico (MEGAHIT / MetaSPAdes)

Ensamblar muestras **paired-end** a partir de lecturas *trimmed*, procesando **una muestra a la vez**, puedes paralelizar este proceso.  
Usa **MEGAHIT** por defecto y **MetaSPAdes** como alternativa.

1. **Localiza pares R1/R2** dentro de `mags/results/02.trimmed/<sample>/**/{R1,R2}*.f*q.gz` con patrones y reemplazo de tokens.
2. **Ordena** las muestras pendientes por **tamaño total** (R1+R2) ascendente.
3. **Ensamblaje** por muestra:
   - **MEGAHIT** (default) con `--presets meta-sensitive` **y** flags adicionales sugeridos.
   - **MetaSPAdes** (`metaspades.py` o `spades.py`) con `--only-assembler` como ejemplo extra.
4. **Normaliza salidas**:
   - MEGAHIT: `final.contigs.fa → contigs.fasta → <sample>_contigs.fasta`.
   - MetaSPAdes: renombra `contigs.fasta`/`scaffolds.fasta` a `<sample>_*`.
5. **Registra logs** por muestra y compone un **resumen CSV**.

## Entradas (inputs)
- Directorio con lecturas **recortadas**: `mags/results/02.trimmed/`
- Archivos esperados:
  - `<sample>/**/R1*.f*q.gz` y `<sample>/**/R2*.f*q.gz`
- **Variables de entorno** (opcionales):
  - `MAG_PROJECT_DIR` o `MAGENTA_DIR`: raíz del proyecto (default: `cwd`).
  - `MAG_ASSEMBLER`: `megahit` o `metaspades` (si no se usa `--assembler`).
  - `MAG_MEGAHIT_THREADS` / `MAG_SPADES_THREADS`: hilos por muestra (fallback: 40).
  - `MAG_SPADES_MEM_GB`: memoria (GB) para MetaSPAdes (0 = auto).

## Salidas (outputs)
- Ensamblajes por muestra: `mags/results/03.assembly/<sample>_assembled/`
  - **MEGAHIT**: `<sample>_contigs.fasta` (a partir de `final.contigs.fa`)
  - **MetaSPAdes**: `<sample>_contigs.fasta` y opcional `<sample>_scaffolds.fasta`
- **Logs**: `mags/results/03.assembly/logs/<sample>.log`
- **Resumen**: `mags/results/03.assembly/resumen_ensamblaje.csv`
  - Columnas: `Muestra`, `Estado` (`OK`, `ERROR`, `SKIPPED`), `Detalle` (timestamp o mensaje)

## Descubrimiento de pares (R1/R2)
- Patrones R1 (recursivo): `*_R1*.fastq.gz`, `*_1*.fq.gz`, `*R1*.fq.gz`, etc.
- R2 por **reemplazo de tokens** (`_R1↔_R2`, `_1↔_2`, `R1↔R2`, `1↔2`) o **búsqueda mejor-esfuerzo** en el mismo directorio (prefiere mayor prefijo compartido).
- Fallback plano: si no hay subcarpetas en `02.trimmed`, deriva el nombre de la muestra del prefijo del archivo.

---

## Ensambladores y **parámetros usados por el script**

### Parámetros — MEGAHIT

- `-1`, `-2` → FASTQ pareados.  
- `-o <sample_outdir>` → carpeta de salida (`<sample>_assembled`).  
- `-t <threads>` → hilos por muestra.  
- `--presets meta-sensitive` → preset sensible para metagenomas.  
- `--min-contig-len 1000` → **nuevo**: contigs mínimos a reportar (ejemplo didáctico).  
- `--k-list 21,29,39,59,79,99,119` → **nuevo**: lista de *k-mers* (ejemplo didáctico).  

Puedes ajustar `--min-contig-len` y `--k-list` según cobertura/heterogeneidad.  

---

### Parámetros — MetaSPAdes

- `-1`, `-2` → FASTQ pareados.  
- `-o <sample_outdir>` → carpeta de salida (`<sample>_assembled`).  
- `-t <threads>` → hilos por muestra.  
- `-m <mem_gb>` → **opcional** si `--mem-gb > 0`; si `0`, deja el autoajuste.  
- `--only-assembler` → **nuevo**: omite la fase de corrección si ya hiciste QC/corrección aguas arriba.  

## Ejemplos prácticos

### MEGAHIT con contigs ≥1 kb y k-list amplio
```bash
python 06_assembly_serial.py --assembler megahit --threads 48
(El script ya incluye --min-contig-len 1000 y --k-list 21,29,39,59,79,99,119 como ejemplo didáctico.)

MetaSPAdes con 64 hilos, 512 GB y solo ensamblado
python 06_assembly_serial.py --assembler metaspades --threads 64 --mem-gb 512


