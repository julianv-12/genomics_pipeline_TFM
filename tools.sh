#!/bin/bash

# Actualizar paquetes
apt-get update

# Instalar herramientas necesarias
echo "Instalando Samtools..."
apt-get install -y samtools

echo "Instalando BWA para alineamiento..."
apt-get install -y bwa

echo "Instalando FastQC..."
if ! command -v fastqc &> /dev/null; then
    apt-get install -y fastqc
else
    echo "FastQC ya está instalado."
fi

echo "Instalando SRA Toolkit..."
apt-get install -y sra-toolkit

echo "Instalando Trimmomatic..."
if ! command -v trimmomatic &> /dev/null; then
    apt-get install -y trimmomatic
else
    echo "Trimmomatic ya está instalado."
fi

# Descargar adaptadores para Trimmomatic
echo "Verificando adaptadores para Trimmomatic..."
if [ ! -f /usr/share/trimmomatic/adapters/TruSeq3-PE.fa ]; then
    wget -O /usr/share/trimmomatic/adapters/TruSeq3-PE.fa https://github.com/timflutre/trimmomatic/raw/master/adapters/TruSeq3-PE.fa
else
    echo "Adaptador TruSeq3-PE.fa ya existe."
fi

# Instalar Cutadapt
echo "Instalando Cutadapt..."
if ! command -v cutadapt &> /dev/null; then
    apt-get install -y pipx
    pipx install cutadapt
else
    echo "Cutadapt ya está instalado."
fi

# Configurar el PATH para Cutadapt
echo 'export PATH=$PATH:/root/.local/bin' >> ~/.bashrc
source ~/.bashrc

# Descargar e instalar BBMap
echo "Instalando BBMap..."
if [ ! -d "/TFM/bbmap" ]; then
    wget -O /TFM/BBMap_39.11.tar.gz https://sourceforge.net/projects/bbmap/files/latest/download
    tar -xzvf /TFM/BBMap_39.11.tar.gz -C /TFM/
    mv /TFM/bbmap /TFM/BBMap_39.11
else
    echo "BBMap ya está instalado."
fi

# Descargar Picard Tools
echo "Instalando Picard Tools..."
if [ ! -f /TFM/picard.jar ]; then
    wget -O /TFM/picard.jar https://github.com/broadinstitute/picard/releases/download/2.27.1/picard.jar
else
    echo "Picard Tools ya está descargado."
fi

echo "Instalación completa."


