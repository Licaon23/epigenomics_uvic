# EN-TEx ChIP-seq data: navigating ENCODE and running `chip-nf`

This document summarizes my work for tutorial **3.1. EN-TEx ChIP-seq data: how to navigate the portal and run the chipnf pipeline** from the `epigenomics_uvic` repository.

The goal of this exercise was to learn how to:

1. Work inside the course Docker environment.
2. Retrieve ENCODE metadata for EN-TEx ChIP-seq experiments.
3. Parse the metadata to identify the FASTQ files needed for a specific experiment.
4. Manually retrieve control FASTQ identifiers when they are not included in the metadata.
5. Prepare the `chip-nf` pipeline index file.
6. Run the `chip-nf` ChIP-seq pipeline using a reduced example dataset.

---

## 1. Starting the Docker environment

The tutorial was run using the Docker image provided for the course. I mounted the current GitHub repository inside the container and set the working directory to the repository path:

```bash
docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

Explanation of the options:

- `-v $PWD:$PWD`: bind-mounts the current local directory into the same path inside the container.
- `-w $PWD`: sets the working directory inside the container.
- `--rm`: removes the container after exiting.
- `-it`: opens the container in interactive mode.

This setup allowed the files created inside the container, such as downloaded metadata, FASTQ files, and pipeline index files, to remain available in the local repository.

---

## 2. Downloading ENCODE metadata

The first step was to download metadata for DNA-binding experiments from the ENCODE portal for a specific EN-TEx donor, restricting the search to:

- donor UUID: `d370683e-81e7-473f-8475-7716d027849b`
- assembly: `GRCh38`
- biosamples: `sigmoid colon` and `stomach`
- assay category: `DNA binding`
- released experiments only

The metadata was downloaded with:

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assembly=GRCh38&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+binding"
```

The ENCODE search returned **28 DNA-binding experiments**, including transcription factor ChIP-seq, histone ChIP-seq, and control experiments.

For this exercise, I focused only on one experiment:

> **H3K4me3 ChIP-seq in sigmoid colon**

This is a histone modification ChIP-seq experiment.

---

## 3. Inspecting the metadata file

Since ChIP-seq pipelines usually start from sequencing output files, the next task was to retrieve the FASTQ file identifiers for the selected experiment.

To understand the structure of the metadata file, I printed the column names together with their positions:

```bash
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```

The relevant parts of the command are:

- `FS=OFS="\t"`: sets the input and output field separators to tab.
- `for (i=1;i<=NF;i++)`: loops over all fields in the header line.
- `print $i, i`: prints each column name and its corresponding position.

Some important columns were:

```text
File_accession          1
File_format             2
File_type               3
File_format_type        4
Output_type             5
File_assembly           6
Experiment_accession    7
Assay                   8
```

For this exercise, the most relevant fields were:

- file accession
- file format / file type
- assay or experiment target
- biosample term name
- biological replicate

---

## 4. Counting FASTQ files for H3K4me3 in sigmoid colon

To count the FASTQ files corresponding to the selected experiment, I filtered the metadata for:

1. rows containing `H3K4me3`
2. rows containing `sigmoid_colon`
3. rows where the file format was `fastq`

```bash
grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
awk 'BEGIN{FS="\t"}$2=="fastq"{n++}END{print n}'
```

An equivalent approach would be:

```bash
grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
awk -F "\t" '$2=="fastq"' | \
wc -l
```

The result was:

```text
4
```

So, there were **4 FASTQ files** available for the H3K4me3 sigmoid colon experiment.

When checking the ENCODE experiment page manually, the experiment also contained other file types, such as:

- FASTQ files
- BAM files
- bigWig files
- BED / bigBed peak files

However, for running the pipeline from raw sequencing data, the FASTQ files were the key input.

---

## 5. Retrieving FASTQ accessions for the ChIP sample

To retrieve the FASTQ accessions, I filtered the metadata and printed column 1, which contains the ENCODE file accession:

```bash
grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
awk 'BEGIN{FS=OFS="\t"}$2=="fastq"{print $1}'
```

The FASTQ accessions were:

```text
ENCFF166KFA
ENCFF963MZL
ENCFF880IAN
ENCFF215QVS
```

These accessions can be used to download the corresponding FASTQ files from ENCODE.

The general download command is:

```bash
grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
awk 'BEGIN{FS=OFS="\t"}$2=="fastq"{print $1}' | \
while read filename; do
  wget -P data/fastq.files "https://www.encodeproject.org/files/$filename/@@download/$filename.fastq.gz"
done
```

The logic is:

- `grep` selects the experiment of interest.
- `awk` extracts only FASTQ file accessions.
- `while read filename` loops over each accession.
- `wget` downloads the corresponding FASTQ file into `data/fastq.files`.

In the practical exercise, subsampled FASTQ files were used instead of the full files to reduce time and storage requirements.

---

## 6. Retrieving control FASTQ files

When checking the ENCODE experiment page, there was an audit warning indicating that the control samples were not reported in the metadata file.

This is important because ChIP-seq analysis requires control/input samples to correctly distinguish true enrichment from background signal.

The control experiment had to be accessed manually through the ENCODE portal. From the control experiment page, I retrieved the following FASTQ accessions:

```text
ENCFF102SFU
ENCFF187FWF
ENCFF549PVM
ENCFF599MFK
ENCFF876EUR
ENCFF950GED
```

These control files can be downloaded with:

```bash
echo -e "ENCFF102SFU\nENCFF599MFK\nENCFF187FWF\nENCFF549PVM\nENCFF876EUR\nENCFF950GED" | \
while read filename; do
  wget -P data/fastq.files/ "https://www.encodeproject.org/files/$filename/@@download/$filename.fastq.gz"
done
```

Again, in the exercise, only a reduced set of subsampled FASTQ files was used.

---

## 7. Preparing the `chip-nf` pipeline input file

The `chip-nf` pipeline requires an index file describing the samples and input FASTQ files.

According to the pipeline format, the index file contains five tab-separated columns:

```text
sample_id    run_id    fastq_path    control_id    mark
```

Conceptually, these columns represent:

| Column | Meaning |
|---|---|
| `sample_id` | Identifier used to merge BAM files belonging to the same biological sample |
| `run_id` | Identifier for a specific sequencing run |
| `fastq_path` | Path to the FASTQ file |
| `control_id` | Identifier of the input/control sample, or `-` if no control is used |
| `mark` | Histone mark, transcription factor, or `input` for control samples |

An example structure is:

```text
sample1    sample1_run1    /path/to/sample1_run1.fastq.gz    control1    H3K4me3
control1   control1_run1   /path/to/control1_run1.fastq.gz   control1    input
```

---

## 8. Creating index entries for the H3K4me3 sample

For the H3K4me3 sigmoid colon sample, I used the metadata fields:

- column 1: file accession
- column 2: file format
- column 11: biosample term name
- column 23: experiment target
- column 35: biological replicate

The following command builds the sample lines for the pipeline index file:

```bash
mypath=$(pwd)

grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
cut -f1,2,11,23,35 | \
awk -v mypath="$mypath" 'BEGIN{FS=OFS="\t";n=0}$2=="fastq"{n++; split($4, a,"-"); print $3"_"a[1]"_"$5, $3"_"a[1]"_"$5"_run"n, mypath"/data/fastq.files/"$1".sub.fastq", "sigmoid_colon_control_1", a[1]}' | \
head -2 > chip-nf/pipeline.index.tsv
```

Explanation of the generated fields:

- `sample_id`: built as `sigmoid_colon_H3K4me3_<replicate_number>`
- `run_id`: same as `sample_id`, with `_run<n>` appended
- `fastq_path`: absolute path to the subsampled FASTQ file
- `control_id`: `sigmoid_colon_control_1`
- `mark`: `H3K4me3`

Only the first two FASTQ files were used with `head -2` to keep the example lightweight.

---

## 9. Adding control entries to the pipeline index

The control FASTQ entries were appended to the same index file:

```bash
echo -e "ENCFF102SFU\nENCFF599MFK\nENCFF187FWF\nENCFF549PVM\nENCFF876EUR\nENCFF950GED" | \
head -2 | \
awk -v mypath="$mypath" 'BEGIN{FS=OFS="\t";n=0}{n++; print "sigmoid_colon_control_1", "sigmoid_colon_control_1_run"n, mypath"/data/fastq.files/"$1".sub.fastq", "sigmoid_colon_control_1", "input"}' >> chip-nf/pipeline.index.tsv
```

Here:

- `>>` appends the control lines to the existing index file.
- `head -2` keeps only two control FASTQ files for the reduced exercise.
- the mark column is set to `input`, because these are control samples.

The final `pipeline.index.tsv` therefore contains:

- two H3K4me3 ChIP FASTQ files
- two control/input FASTQ files

---

## 10. Running the `chip-nf` pipeline

The pipeline was run from the `chip-nf` directory:

```bash
cd ChIP-seq/chip-nf
```

The command used was:

```bash
NXF_VER=22.04.0 nextflow run guigolab/chip-nf \
  -dsl1 \
  -bg \
  -r v0.2.3 \
  --index pipeline.index.tsv \
  --genome refs/GRCh38.primary_assembly.genome.chr19.fa \
  --genomeSize hs \
  -with-docker > pipeline.log.txt
```

Explanation of the main arguments:

- `NXF_VER=22.04.0`: runs the pipeline with a specific Nextflow version.
- `nextflow run guigolab/chip-nf`: runs the `chip-nf` pipeline from the Guigó Lab repository.
- `-dsl1`: uses Nextflow DSL1 syntax, required by this pipeline version.
- `-bg`: runs the pipeline in the background.
- `-r v0.2.3`: uses release `v0.2.3` of the pipeline.
- `--index pipeline.index.tsv`: provides the sample index file.
- `--genome refs/GRCh38.primary_assembly.genome.chr19.fa`: uses the chromosome 19 reference FASTA.
- `--genomeSize hs`: indicates the human genome size setting.
- `-with-docker`: runs the pipeline processes using Docker.
- `> pipeline.log.txt`: saves the pipeline log to a text file.

For time and computational reasons, the exercise used only chromosome 19 instead of the full human genome.

After launching the pipeline, I returned to the parent directory with:

```bash
cd ../../..
```

---

## 11. Pipeline outputs

The pipeline generated several types of output files. The output types listed in my notes include:

```text
Alignments
narrowPeak
broadPeak
gappedPeak
fcSignal
pileupSignal
pvalueSignal
zeroneMatrix
zeroneBed
zeroneMergedBed
```

Broadly, these outputs correspond to:

| Output type | Interpretation |
|---|---|
| `Alignments` | mapped reads, usually BAM-related outputs |
| `narrowPeak` | narrow peak calls |
| `broadPeak` | broad peak calls |
| `gappedPeak` | peak calls allowing gapped regions |
| `fcSignal` | fold-change signal tracks |
| `pileupSignal` | read pileup signal tracks |
| `pvalueSignal` | statistical signal tracks |
| `zeroneMatrix` | matrix used by Zerone |
| `zeroneBed` | Zerone peak output in BED format |
| `zeroneMergedBed` | merged Zerone peak regions |

These files summarize different stages of the ChIP-seq analysis, including mapping, signal generation, and peak calling.

---

## Useful commands summary

### Run the course Docker container

```bash
docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

### Download metadata

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&assembly=GRCh38&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assay_slims=DNA+binding"
```

### Inspect metadata columns

```bash
head -1 metadata.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```

### Count FASTQ files

```bash
grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
awk 'BEGIN{FS="\t"}$2=="fastq"{n++}END{print n}'
```

### Download H3K4me3 FASTQ files

```bash
grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
awk 'BEGIN{FS=OFS="\t"}$2=="fastq"{print $1}' | \
while read filename; do
  wget -P data/fastq.files "https://www.encodeproject.org/files/$filename/@@download/$filename.fastq.gz"
done
```

### Download control FASTQ files

```bash
echo -e "ENCFF102SFU\nENCFF599MFK\nENCFF187FWF\nENCFF549PVM\nENCFF876EUR\nENCFF950GED" | \
while read filename; do
  wget -P data/fastq.files/ "https://www.encodeproject.org/files/$filename/@@download/$filename.fastq.gz"
done
```

### Create sample entries for the pipeline index

```bash
mypath=$(pwd)

grep -F H3K4me3 metadata.tsv | \
grep -F sigmoid_colon | \
cut -f1,2,11,23,35 | \
awk -v mypath="$mypath" 'BEGIN{FS=OFS="\t";n=0}$2=="fastq"{n++; split($4, a,"-"); print $3"_"a[1]"_"$5, $3"_"a[1]"_"$5"_run"n, mypath"/data/fastq.files/"$1".sub.fastq", "sigmoid_colon_control_1", a[1]}' | \
head -2 > chip-nf/pipeline.index.tsv
```

### Add control entries to the pipeline index

```bash
echo -e "ENCFF102SFU\nENCFF599MFK\nENCFF187FWF\nENCFF549PVM\nENCFF876EUR\nENCFF950GED" | \
head -2 | \
awk -v mypath="$mypath" 'BEGIN{FS=OFS="\t";n=0}{n++; print "sigmoid_colon_control_1", "sigmoid_colon_control_1_run"n, mypath"/data/fastq.files/"$1".sub.fastq", "sigmoid_colon_control_1", "input"}' >> chip-nf/pipeline.index.tsv
```

### Run `chip-nf`

```bash
cd ChIP-seq/chip-nf

NXF_VER=22.04.0 nextflow run guigolab/chip-nf \
  -dsl1 \
  -bg \
  -r v0.2.3 \
  --index pipeline.index.tsv \
  --genome refs/GRCh38.primary_assembly.genome.chr19.fa \
  --genomeSize hs \
  -with-docker > pipeline.log.txt

cd ../../..
```
