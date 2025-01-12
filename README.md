# genomics_pipeline_TFM
Pipeline Nextflow para análisis genómico en el TFM

Este repositorio contiene un pipeline de análisis bioinformático desarrollado como parte del Trabajo Final de Máster (TFM). Utiliza herramientas bioinformáticas modernas como **Nextflow**, **FastQC**, **Cutadapt**, y otras para realizar un análisis completo, desde el preprocesamiento de lecturas hasta la anotación y priorización de variantes genéticas.

## **Características principales**
- **Preprocesamiento de lecturas**: Incluye control de calidad, eliminación de adaptadores, deduplicación y normalización.
- **Alineamiento**: Lecturas alineadas al genoma de referencia (GRCh38).
- **Llamado de variantes**: Identificación de variantes genómicas con herramientas como GATK y BCFtools.
- **Anotación funcional**: Predicción del impacto funcional utilizando herramientas como ANNOVAR.

---

## **Requisitos previos**
Para ejecutar este pipeline, necesitas las siguientes herramientas instaladas:
- **[Nextflow](https://www.nextflow.io/)**: Para orquestar y ejecutar el pipeline.
- **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**: Para realizar el control de calidad de las lecturas.
- **[Cutadapt](https://cutadapt.readthedocs.io/)**: Para la eliminación de adaptadores.
- **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**: Para mejorar la calidad de las lecturas.
- **[SAMtools](http://www.htslib.org/)**: Para manipular archivos SAM y BAM.
- **[GATK](https://gatk.broadinstitute.org/)**: Para el llamado de variantes.
- **[ANNOVAR](https://annovar.openbioinformatics.org/)**: Para la anotación funcional de variantes.
- **[BBMap](https://sourceforge.net/projects/bbmap/)**: Herramienta de normalización y deduplicación.
- **[BCFtools](http://www.htslib.org/doc/bcftools.html)**: Análisis de variantes.

### Descarga del genoma humano:
Para realizar el alineamiento y el llamado de variantes, necesitas el genoma humano de referencia. Se recomienda descargarlo desde **ENSEMBL**:

1. Descarga el genoma de referencia (GRCh38):
   - URL: https://www.ensembl.org/info/data/ftp/index.html
   - Archivo necesario: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`

2. Asegúrate de indexar el genoma con las herramientas adecuadas:
   ```bash
   samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
   bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
---

## **Instalación**
1. Clona este repositorio:
   ```bash
   git clone https://github.com/julianv-12/genomics_pipeline_TFM.git
   cd genomics_pipeline_TFM

## Estructura del Proyecto

El proyecto está organizado de la siguiente manera:

genomics_pipeline_TFM/
├── aligned_2/             # Archivos alineados
├── normalized_2/          # Lecturas normalizadas
├── fastqc_results_2/      # Resultados de FastQC
├── variants_2/            # Variantes detectadas
├── trimmomatic_results_2/ # Resultados de Trimmomatic
├── fastq_2/               # Lecturas preprocesadas
├── genomics_pipeline.nf   # Pipeline principal
├── tools.sh               # Script auxiliar de instalación de herramientas
└── README.md              # Documentación del proyecto

## Ejecucion del pipelline
nextflow run genomics_pipeline.nf


### Notas sobre la estructura

- Cada carpeta contiene resultados generados por diferentes herramientas.
- El archivo `genomics_pipeline.nf` es el pipeline principal que puedes ejecutar con Nextflow.
- `tools.sh` incluye comandos útiles para instalar las herramientas necesarias.

