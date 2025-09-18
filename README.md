# MAGENTA 🦠


<div align="justify">
Los manglares son una conocida reserva de diversidad biológica y un ecosistema altamente productivo. Diversos estudios metagenómicos en diferentes partes del mundo han reconocido a la comunidad microbiana del manglar como un agente importante dentro de los ciclos biogeoquímicos, en los cuales se llevan a cabo procesos tales como la transformación del carbono, la fotosíntesis, la fijación de nitrógeno y la reducción de azufre. En la actualidad, no contamos con una herramienta informática que nos permita entender estos procesos y relaciones a una  <b>escala global </b>.
  <br><br>
 <b>MAGENTA </b> (o Global MAngrove GENe CaTAlogue) actúa como un catálogo global de genes únicos y no redundantes a nivel de especie (agrupados al 95% de identidad de nucleótidos) que, a partir de los datos de libre acceso disponibles en bases de datos especializadas (WGS, metagenomas de acceso público - ENA) de cinco de los principales hábitats de la comunidad microbiana del manglar (rizosfera, agua de mar, sedimento, suelo, humedal), busca formular nuevas hipótesis sobre la abundancia, distribución y funciones metabólicas de los microorganismos en el ecosistema del manglar, con miras a atender esta necesidad.

</div>

## Descarga y enriquecimiento de metadatos

<div align="justify">
Al integrar metagenomas públicos, incorporamos sesgos de curación heterogénea (especialmente en NCBI). Usando coordenadas como eje común, incorporamos metadatos de repositorios globales y de acceso publico. Su contexto ambiental explicito altitud, clima (BIO1–BIO19), ecorregión marina, huella humana, tipo de suelos, pH, secuestro de carbono, cobertura vegetal (canopy), y diversidad de especies, sitúan a cada metadato en su nicho ambiental.  Esta estandarización de metadatos permite modelos comparativos y reproducibles, con el potencial de generar información útil en el desarrollo de estrategias de conservación, restauración, y la generación de hipótesis funcionales
</div>



## Reproducibilidad y entorno

- Requisitos (sugeridos): `python>=3.10`, `pandas`, `requests`, `tqdm`, `rasterio`, `geopandas`, `numpy`, `elevation`, `biopython`, `entrez-direct`, `sra-tools`, `aria2c`/`wget`, opcionalmente `aspera`.
- Sistema de archivos: el flujo se diseñó para coexistir con dos estilos de rutas previamente usadas en MAGENTA.

  
