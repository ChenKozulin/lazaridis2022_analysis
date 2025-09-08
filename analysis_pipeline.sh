#install the full bioinformatics toolkit
# Create and activate environment
conda create -n adna_analysis python=3.9 -y
conda activate adna_analysis

# Install everything at once
conda install -y -c bioconda -c conda-forge \
    fastqc adapterremoval bwa samtools bcftools picard gatk4 \
    mapdamage2 angsd admixture eigensoft r-base r-ggplot2 bc

# Create working directory structure
mkdir lazaridis2022_analysis
cd lazaridis2022_analysis
mkdir raw_data processed_data results

# Move your downloaded fastq.gz files to raw_data
mv /path/to/your/downloaded/*fastq.gz raw_data/

# List your samples
ls raw_data/*.fastq.gz > sample_list.txt

# Quick quality check with FastQC on a subset
fastqc *.fastq.gz -o ../results/ -t 4

# Check read counts and basic statistics
echo "Sample Read_Counts File_Type" > ../sample_read_counts.txt

for file in *.fastq.gz; do
    sample=$(basename $file .fastq.gz)  # get filename without extension
    reads=$(gzcat $file | wc -l | awk '{print $1/4}')  # use gzcat on macOS
    
    if [[ $file == *"_R1_"* ]] || [[ $file == *"_R1."* ]] || [[ $file == *"_1."* ]]; then
        echo "$sample $reads R1"
    elif [[ $file == *"_R2_"* ]] || [[ $file == *"_R2."* ]] || [[ $file == *"_2."* ]]; then
        echo "$sample $reads R2"
    else
        echo "$sample $reads Single-end"
    fi
done >> ../sample_read_counts.txt

echo "Sample read count summary:"
cat ../sample_read_counts.txt

#Process Samples with AdapterRemoval
# Create processing script for all samples
ls ERR*.fastq.gz | sed 's/\.fastq\.gz//' > sample_ids.txt
cat sample_ids.txt

#Get read counts for all samples:
echo "Sample Read_Count" > ../sample_read_counts.txt
for file in ERR*.fastq.gz; do
    sample=$(echo $file | sed 's/\.fastq\.gz//')
    reads=$(gunzip -c $file | wc -l | awk '{print $1/4}')
    echo "$sample $reads" >> ../sample_read_counts.txt
done

cat ../sample_read_counts.txt

mkdir -p ../processed_data

while read sample; do
    echo "Processing sample: $sample"
    
    # Single-end processing for your ERR files
    AdapterRemoval --file1 ${sample}.fastq.gz \
      --basename ../processed_data/${sample} \
      --minlength 25 \
      --maxns 10 \
      --minquality 20 \
      --trimns \
      --threads 4
    
done < sample_ids.txt

# Compile preprocessing statistics
echo -e "Sample\tRaw_Reads\tTrimmed_Reads\tRetention_Rate" > preprocessing_stats.txt

while read sample; do
    json_file="processed_data/${sample}.json"
    if [ -f "$json_file" ]; then
        raw_reads=$(grep -A 2 '"input":' $json_file | grep '"reads":' | grep -o '[0-9]*')
        retained_reads=$(grep -A 2 '"output":' $json_file | grep '"reads":' | grep -o '[0-9]*')
        
        if [ "$raw_reads" -gt 0 ] 2>/dev/null; then
            retention_rate=$(echo "scale=2; $retained_reads * 100 / $raw_reads" | bc)
            echo -e "$sample\t$raw_reads\t$retained_reads\t${retention_rate}%"
        fi
    else
        echo "JSON file not found for $sample"
    fi
done < sample_ids.txt >> preprocessing_stats.txt

cat preprocessing_stats.txt

# Download human reference genome (GRCh38/hg38)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz

# Index the reference
bwa index GCA_000001405.15_GRCh38_genomic.fna
samtools faidx GCA_000001405.15_GRCh38_genomic.fna

#3.2 Alignment with Ancient DNA Parameters
# Create results directory if it doesn't exist
mkdir -p results

# Run alignment with corrected paths
while read sample; do
    echo "Aligning sample: $sample"
    
    # Use the trimmed single-end file from AdapterRemoval
    input_file="processed_data/${sample}.r1.fastq.gz"
    
    if [ -f "$input_file" ]; then
        echo "Processing $input_file"
        
        # BWA alignment with correct reference path
        bwa aln -l 1024 -n 0.01 -o 2 -t 4 ../reference/GCA_000001405.15_GRCh38_genomic.fna $input_file > results/${sample}.sai
        
        # Convert to BAM
        bwa samse ../reference/GCA_000001405.15_GRCh38_genomic.fna results/${sample}.sai $input_file | \
          samtools view -Sb -F 4 - > results/${sample}_raw.bam
        
        # Check if files were created
        if [ -s "results/${sample}.sai" ]; then
            echo "SAI file created successfully for $sample"
        else
            echo "ERROR: SAI file is empty for $sample"
        fi
        
        if [ -s "results/${sample}_raw.bam" ]; then
            echo "BAM file created successfully for $sample"
        else
            echo "ERROR: BAM file is empty for $sample"
        fi
    else
        echo "Input file not found: $input_file"
    fi
    
done < sample_ids.txt

# Alignment 
# Check if BWA index exists
ref_genome="reference/GCA_000001405.15_GRCh38_genomic.fna"

echo "Starting alignment pipeline..."
if [ ! -f "${ref_genome}.bwt" ]; then
    echo "ERROR: BWA index not found. Creating index..."
    echo "This will take 30-60 minutes..."
    bwa index "$ref_genome"
fi

# Process each sample
while read sample; do
    sample=$(echo "$sample" | tr -d '\r\n ')
    echo "=== Aligning sample: $sample ==="
    
    input_file="processed_data/${sample}.r1.fastq.gz"
    
    if [ -f "$input_file" ]; then
        echo "Processing $input_file"
        
        echo "Running BWA aln..."
        bwa aln -l 1024 -n 0.01 -o 2 -t 4 "$ref_genome" "$input_file" > "results/${sample}.sai"
        
        echo "Converting to BAM..."
        bwa samse "$ref_genome" "results/${sample}.sai" "$input_file" | \
          samtools view -Sb -F 4 - > "results/${sample}_raw.bam"
        
        echo "Checking results..."
        if [ -s "results/${sample}.sai" ]; then
            echo "SUCCESS: SAI file created for $sample"
            echo "SAI size: $(ls -lh results/${sample}.sai | awk '{print $5}')"
        else
            echo "ERROR: SAI file empty for $sample"
        fi
        
        if [ -s "results/${sample}_raw.bam" ]; then
            echo "SUCCESS: BAM file created for $sample"
            echo "BAM size: $(ls -lh results/${sample}_raw.bam | awk '{print $5}')"
        else
            echo "ERROR: BAM file empty for $sample"
        fi
    else
        echo "ERROR: Input file not found: $input_file"
    fi
    
    echo "=== Finished sample: $sample ==="
    echo ""
    
done < sample_ids.txt

echo "All samples processed!"



