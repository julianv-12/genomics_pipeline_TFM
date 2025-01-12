#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "/TFM/genomics_pipeline/data/SRR*_{1,2}.fastq"
params.outdir = "/TFM/genomics_pipeline/results/fastqc_results"

// Define el canal de lecturas, asegurando que se encuentren archivos fastq
Channel
    .fromFilePairs(params.reads, size: 2)
    .ifEmpty { error "No se encontraron archivos fastq." }
    .set { reads_ch }


//Se define el flujo de trabajo principal

workflow {
// Proceso para realizar FastQC antes de cualquier recorte
runFastQC(reads_ch)

// Proceso para recortar adaptadores usando Cutadapt
cutadapt_results = CUTADAPT(reads_ch)

// Proceso para mejorar la calidad de las lecturas con Trimmomatic
trimmomatic_results = TRIMMOMATIC(cutadapt_results.collect())

// Proceso para realizar FastQC después del recorte
fastqc_trimmomatic = FastQC_PostTrim(trimmomatic_results)

// Proceso para eliminar duplicados despues de aplicar trimmomatic
dedup_results = DeduplicateReads(trimmomatic_results)

// Proceso para normalizar
normalized_results = NormalizeReads(dedup_results) 

// Proceso para realizar FastQC en lecturas normalizadas
fastqc_normalized = FastQC_Normalized(normalized_results)
 
// Alineamiento de lecturas normalizadas 
aligned_results = AlignReads(normalized_results)

// Ordenar las lecturas 
sorted_bam = SortBam(aligned_results)

// Proceso para agregar los readGroups
read_group_bam = AddReadGroup(sorted_bam)

// Proceso para el llamado de variantes
variants = CallVariants(read_group_bam)

//Generar estadisticas 
variant_stats = GenerateVariantStats(variants)

//Filtrar variantes de baja calidad 
filtered_variants = FilterLowQualityVariants(variants)

//Generar estadisticas de las variables filtradas
filtered_variant_stats = GenerateFilteredVariantStats(filtered_variants)

// Anotar variantes
annotated_variants = AnnotateVariants(filtered_variants)


// Preparar archivo para visualización en IGV
    igv_index = PrepareForIGV(filtered_variants_gz)
}

process runFastQC {
    tag "${sample_id}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    
    """
    echo "Ejecutando FastQC para ${sample_id}..."
    fastqc -o . ${reads.join(' ')}
    """
}


process CUTADAPT {
     tag "${pair_id}"
     publishDir "/TFM/fastq_2", mode: 'copy'  // Especifica dónde se guardarán los archivos de salida. 

     input:
     tuple val(pair_id), path(reads)

     output:
     tuple val(pair_id), path("${pair_id}_1_trimmed.fastq"), path("${pair_id}_2_trimmed.fastq"), emit: trimmed_files

     script:
    """
    echo "Eliminando adaptadores Poly-A con Cutadapt..."
    cutadapt -a "A{10}" -A "A{10}" \
      -o ${pair_id}_1_trimmed.fastq -p ${pair_id}_2_trimmed.fastq \
      ${reads[0]} ${reads[1]}

    echo "Verificación de archivos creados:"
    ls -l ${pair_id}_1_trimmed.fastq ${pair_id}_2_trimmed.fastq
    """
}

process TRIMMOMATIC {

     tag "${sample_id}"
     publishDir "/TFM/trimmomatic_results_2", mode: 'copy'

     input:
     tuple val(sample_id), path(read1), path(read2)

     output:
     tuple val(sample_id), path("${sample_id}_1_paired.fastq"), path("${sample_id}_2_paired.fastq"), emit: trimmed_paired

     script:
    """
     echo "Ejecutando Trimmomatic para mejorar la calidad de las lecturas..."
     java -jar /usr/share/java/trimmomatic-0.39.jar PE -phred33 \\
        $read1 $read2 \\
        ${sample_id}_1_paired.fastq ${sample_id}_1_unpaired.fastq \\
        ${sample_id}_2_paired.fastq ${sample_id}_2_unpaired.fastq \\
        ILLUMINACLIP:/usr/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \\
        LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

    """     
}

process FastQC_PostTrim {
    tag "${sample_id}"
    publishDir "/TFM/fastqc_results_2", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*.html", emit: qc_html
    path "*.zip", emit: qc_zip

    script:
    """
    fastqc -o . ${read1} ${read2}
    """
}

process DeduplicateReads {
    tag "${sample_id}"
    publishDir "/TFM/dedupe_results_2", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:    
    tuple val(sample_id), path ("${sample_id}_interleaved_dedup.fastq"), emit: interleaved_dedup
    
    script:
    """
    if [[ ! -f $read1 ]] || [[ ! -f $read2 ]]; then
        echo "Error: Archivos de entrada faltantes para ${sample_id}."
        exit 1
    fi
    echo "Archivos de entrada verificados para ${sample_id}. Procediendo con la deduplicación..."
    /TFM/BBMap_39.11/bbmap/dedupe.sh in1=$read1 in2=$read2 out=${sample_id}_interleaved_dedup.fastq
    """
    
}     


process NormalizeReads {
    tag "${sample_id}"
    publishDir "/TFM/normalized_2", mode: 'copy'

    input:
    tuple val(sample_id), path(interleaved_dedup) 
     

    output:
    tuple val(sample_id), path ("${sample_id}_norm.fastq.gz"), emit: normalized_fastq
    
    script:
    """
    echo "Iniciando la normalización de lecturas para ${sample_id}..."
    /TFM/BBMap_39.11/bbnorm.sh in=$interleaved_dedup  out=${sample_id}_norm.fastq.gz target=150 min=5
        echo "Normalización completada para ${sample_id}."
    """
}

process FastQC_Normalized {
    tag "${sample_id}"
    publishDir "/TFM/fastqc_results_2", mode: 'copy'

    input:
    tuple val(sample_id), path(norm_reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    echo "Ejecutando FastQC en lecturas normalizadas para ${sample_id}..."
    fastqc -o . ${norm_reads}
    """
}

process AlignReads {
    tag "${sample_id}"
    publishDir "/TFM/aligned_2", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ("${sample_id}_aligned.sam"), emit: aligned_sam

    script:
    """
    echo "Alineando lecturas al genoma de referencia para ${sample_id}..."
    bwa mem /TFM/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa $reads > ${sample_id}_aligned.sam
    """
}

process SortBam {
    tag "${sample_id}"
    publishDir "/TFM/aligned_2", mode: 'copy'

    input:
    tuple val(sample_id), path(aligned_sam)

    output:
    tuple val(sample_id), path ("${sample_id}_sorted.bam"), emit: sorted_bam

    script:
    """
    echo "Convirtiendo SAM a BAM y ordenando para ${sample_id}..."
    samtools view -S -b $aligned_sam > ${sample_id}_aligned.bam
    samtools sort ${sample_id}_aligned.bam -o ${sample_id}_sorted.bam
    """
}

process AddReadGroup {
    tag "${sample_id}"
    publishDir "/TFM/aligned_2", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path ("${sample_id}_sorted_rg.bam"), emit: sorted_rg_bam

    script:
    """
    echo "Agregando Read Group ID a ${sample_id}..."
    samtools addreplacerg -r 'ID:${sample_id}' -r 'LB:lib1' -r 'PL:ILLUMINA' -r 'PU:unit1' -r 'SM:${sample_id}' -o ${sample_id}_sorted_rg.bam $sorted_bam
    """
}

process CallVariants {
    tag "${sample_id}"
    publishDir "/TFM/variants_2", mode: 'copy'

    input:
    tuple val(sample_id), path(rg_bam)

    output:
    tuple val(sample_id), path ("${sample_id}_variants.vcf.gz"), emit: vcf

    script:
    """
    echo "Llamando variantes con bcftools para ${sample_id}..."
    bcftools mpileup -Ou -f /TFM/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa $rg_bam | \
    bcftools call -mv -Oz -o ${sample_id}_variants.vcf.gz
    """
}

process GenerateVariantStats {
    tag "${sample_id}"
    publishDir "/TFM/variants_2", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    tuple val(sample_id), path ("${sample_id}_variants_stats.txt"), emit: vcf_stats

    script:
    """
    echo "Generando estadísticas del archivo de variantes para ${sample_id}..."
    bcftools stats $vcf_file > ${sample_id}_variants_stats.txt
    """
}

process FilterLowQualityVariants {
    tag "${sample_id}"
    publishDir "/TFM/variants_2", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    tuple val(sample_id), path ("${sample_id}_variants_filtered.vcf.gz"), emit: filtered_vcf

    script:
    """
    echo "Filtrando variantes de baja calidad para ${sample_id}..."
    bcftools filter -e 'QUAL<30' -Oz -o ${sample_id}_variants_filtered.vcf.gz $vcf_file
    """
}
process GenerateFilteredVariantStats {
    tag "${sample_id}"
    publishDir "/TFM/variants_2", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf)

    output:
    tuple val(sample_id), path ("${sample_id}_variants_filtered_stats.txt"), emit: filtered_stats

    script:
    """
    echo "Generando estadísticas del archivo de variantes filtrado para ${sample_id}..."
    bcftools stats $filtered_vcf > ${sample_id}_variants_filtered_stats.txt
    """
}


process AnnotateVariants {
    tag "${sample_id}"
    publishDir "/TFM/variants_2", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf)

    output:
    path "${sample_id}_annotated_variants.hg38_multianno.csv", emit: annotated_csv

    script:
    """
    echo "Iniciando la anotación de variantes para ${sample_id}..."
    bcftools convert $filtered_vcf -O u -o ${sample_id}_filtered.avinput
    perl /TFM/annovar/table_annovar.pl ${sample_id}_filtered.avinput /TFM/annovar/humandb/ \
        -buildver hg38 \
        -out ${sample_id}_annotated_variants \
        -remove \
        -protocol refGene,exac03,dbnsfp33a \
        -operation g,f,f \
        -nastring . \
        -csvout
    """
}

process PrepareForIGV {
    tag "${sample_id}"
    publishDir "/TFM/variants_2", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf_gz)

    output:
    path "${sample_id}_variants_filtered.vcf.gz.tbi", emit: igv_index

    script:
    """
    echo "Indexando archivo para visualización en IGV..."
    tabix -p vcf ${filtered_vcf_gz}
    """
}

