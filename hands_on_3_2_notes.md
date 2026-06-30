# ChIP-seq data: downstream analyses

In the previous section, we saw how to navigate the ENCODE portal and retrieve our data of interest. We also learned how to run the `chip-nf` pipeline to process raw ChIP-seq data.

Here, we will work with processed ChIP-seq data available on the ENCODE portal for the same experiment: H3K4me3 ChIP-seq in sigmoid colon. We will compare these results with H3K4me3 ChIP-seq experiments performed in stomach tissue and learn some commonly performed downstream analyses in epigenomics projects.

## 1. Download processed ChIP-seq files

First:

1. Create an `analyses` folder in the `ChIP-seq` directory.
2. Download the data for the peak calling files (`bigBed`) and the FC signal files (`bigWig`) from the repository, using the metadata file.

### 1.1 Download `bigBed` peak files

```bash
 grep -F H3K4me3 metadata.tsv |\
 grep -F "bigBed_narrowPeak" |\
 grep -F "pseudoreplicated_peaks" |\
 grep -F "GRCh38" |\
 awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
 sort -k2,2 -k1,1r |\
 sort -k2,2 -u > analyses/bigBed.peaks.ids.txt
```

```bash
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

Notes:

- `pseudoreplicated_peaks` are the peak calling files after QC, used for reproducibility.
- There are two files for each H3K4me3 experiment in stomach and sigmoid colon. We use `sort` to sort them by tissue type, and then by file ID in descending order, to have the most recent file of each tissue on top of the list (`r` flag). Then, we use `sort -u` on column 2 to get a unique row, keeping the first one, which corresponds to the most recent file.
- After that, we get and save the file ID, and use the `wget` command with the corresponding URL to download the files. We use the `while` structure previously used.

### 1.2 Download `bigWig` FC signal files

The same process is used for the `bigWig` files containing the FC signal.

```bash
grep -F H3K4me3 metadata.tsv |\
grep -F "bigWig" |\
grep -F "fold_change_over_control" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigWig.FC.ids.txt
```

```bash
cut -f1 analyses/bigWig.FC.ids.txt |\
while read filename; do
  wget -P data/bigWig.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigWig"
done
```

### 1.3 Check file integrity using MD5 checksums

```bash
for file_type in bigBed bigWig; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```

## 2. Create an aggregation plot

H3K4me3 typically marks promoters of actively expressed genes. This plot will allow us to visualize the FC signal, which represents how enriched a DNA region is for this specific protein compared to controls, and compare the signal between highly and lowly expressed genes.

We will restrict the analysis to protein-coding genes. For each tissue, we will split them into two groups: highly expressed and lowly expressed genes.

To do so, we need:

1. The `bigWig` FC file downloaded before.
2. A GENCODE annotation file to have the list of protein-coding genes and genomic coordinates over which we will compute the aggregated signal.
3. Expression matrices of stomach and sigmoid colon tissues. From these matrices, we will select the top 1,000 most and least expressed genes.

### 2.1 Download the GENCODE annotation

We need the same GENCODE version used by ENCODE to generate the expression matrices from En-TEx RNA-seq experiments for stomach and sigmoid colon:

- Assay type: transcription.
- Assay title: total RNA-seq.
- Biosample names: stomach and sigmoid colon.

We will use assembly GRCh38 and GENCODE annotation v24.

Download the annotation file: `gencode.v24.primary_assembly.annotation.gtf.gz`.

```bash
mkdir annotation
wget -P annotation "https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz"
```

Then, uncompress the GTF file and prepare it in BED format.

```bash
gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz
```

### 2.2 Explore the GTF file format

The GTF file contains the following fields:

```text
1   chr1
2   HAVANA
3   gene
4   11869
5   14409
6   .
7   +
8   .
9   gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
```

These correspond to chromosome name, source, feature type, genomic positions, score, strand, frame, and attributes. The third column reports the type of region: `CDS`, `exon`, `gene`, `UTR`, `stop_codon`, `transcript`, etc.

### 2.3 Create a BED file with protein-coding gene body coordinates

We need to retrieve gene body coordinates of protein-coding genes:

- chromosome: column 1
- start: column 4
- end: column 5
- strand: column 7

We also remove mitochondrial genes located on `chrM` and convert the coordinate system from 1-based to 0-based. This creates a BED file for protein-coding genes.

```bash
awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\     
grep -F "protein_coding" |\    
cut -d ";" -f1 |\  
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\  
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed
```

Explanation:

- The first command checks whether the third column contains `gene`.
- Then, we filter rows containing `protein_coding`.
- The attributes in the last column are separated by `;`. The `cut` instruction establishes `;` as a delimiter and keeps the first field, as if it was all in the same column with tab delimiters.
- Then, we take the columns of interest and print them.
- `sed` is a stream editor, a text transformer. Here, the substitution command `s/\"//g` removes the literal double quotes from the annotation.
- Then, we select rows that do not contain `chrM` in the chromosome column and transform the coordinate system from 1-based GTF to 0-based BED. A 1-based coordinate puts the coordinate exactly on the base pair, whereas 0-based coordinates are in between. We subtract 1 from the start coordinate.
- Finally, we pipe the output to a new file in BED format.

Example output:

```text
chr1    69090   70008   ENSG00000186092.4  0   +   ENSG00000186092.4
chr1    182392  184158  ENSG00000279928.1  0   +   ENSG00000279928.1
chr1    184922  200322  ENSG00000279457.3  0   -   ENSG00000279457.3
chr1    450739  451678  ENSG00000278566.1  0   -   ENSG00000278566.1
chr1    685715  686654  ENSG00000273547.1  0   -   ENSG00000273547.1
chr1    924879  944581  ENSG00000187634.10 0   +   ENSG00000187634.10
chr1    944203  959309  ENSG00000188976.10 0   -   ENSG00000188976.10
chr1    960586  965715  ENSG00000187961.13 0   +   ENSG00000187961.13
chr1    966496  975865  ENSG00000187583.10 0   +   ENSG00000187583.10
chr1    975203  982093  ENSG00000187642.9  0   -   ENSG00000187642.9
```

The BED file for the annotated protein-coding genes is saved in the `annotation` directory.

## 3. Retrieve highly and lowly expressed protein-coding genes for each tissue

To know the expression levels of genes, we need to retrieve an expression matrix from RNA-seq experiments.

This information is not in the metadata file, since these are from DNA binding experiments, not transcription assays. We retrieve the files manually from:

```text
https://www.encodeproject.org/search/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&limit=all&assay_slims=Transcription&assay_title=total+RNA-seq&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon
```

```bash
echo -e "ENCFF268RWA\tsigmoid_colon\nENCFF918KPC\tstomach" > analyses/tsv.totalRNASeq.ids.txt
```

### 3.1 Download the expression matrices

Download the expression matrices just like before.

```bash
mkdir data/tsv.files
```

```bash
cut -f1 analyses/tsv.totalRNASeq.ids.txt |\
while read filename; do
  wget -P data/tsv.files "https://www.encodeproject.org/files/$filename/@@download/$filename.tsv"
done
```

### 3.2 Explore one total RNA-seq gene expression matrix

```bash
head data/tsv.files/ENCFF268RWA.tsv | column -t
```

```text
gene_id transcript_id(s) length effective_length expected_count TPM FPKM posterior_mean_count posterior_standard_deviation_of_count pme_TPM pme_FPKM TPM_ci_lower_bound TPM_ci_upper_bound FPKM_ci_lower_bound FPKM_ci_upper_bound
10904   10904            93.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12954   12954            94.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12956   12956            72.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12958   12958            82.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12960   12960            73.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12962   12962            72.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12964   12964            74.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12965   12965            82.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
12967   12967            73.00  0.00             0.00           0.00 0.00 0.00                 0.00                                  0.00    0.00     0                  0                  0                   0
```

Which column contains the TPM value?

```bash
head -1 data/tsv.files/ENCFF268RWA.tsv | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++){print $i, i}}'
```

Column 6.

This matrix contains the whole set of genes, not only protein-coding genes. How many are Ensembl genes?

```bash
grep "^ENSG" data/tsv.files/ENCFF268RWA.tsv | wc -l
```

```text
60725
```

### 3.3 Subset expression matrices to protein-coding genes

We will work only with protein-coding genes. We subset the expression matrices to this set of genes, and retrieve the gene ID and TPM information.

```bash
cut -f1 analyses/tsv.totalRNASeq.ids.txt |\
while read filename; do 
  ../bin/selectRows.sh <(cut -f4 annotation/gencode.v24.protein.coding.gene.body.bed) <(cut -f1,6 data/tsv.files/"$filename".tsv) > tmp
  mv tmp data/tsv.files/"$filename".tsv
done
```

Steps:

- Get the accession number of each total RNA-seq downloaded data file, for stomach and sigmoid colon.
- For each file:
  - Get the gene IDs from the GENCODE BED file that are protein-coding.
  - Get the gene ID and TPM from the total RNA-seq expression experiment.
  - Select rows from the second file that are present in the first file, obtaining TPM expression values only for genes that are protein-coding.
  - Save this into a temporary file.
  - Finally, move the temporary file to rewrite the original `.tsv` RNA expression file using `mv`.

Example output:

```text
ENSG00000000003.14  0.21
ENSG00000000005.5   0.01
ENSG00000000419.12  1.46
ENSG00000000457.13  0.37
ENSG00000000460.16  0.20
ENSG00000000938.12  0.33
ENSG00000000971.15  1.04
ENSG00000001036.13  0.68
ENSG00000001084.10  1.47
ENSG00000001167.14  0.35
```

### 3.4 Select the 1,000 most expressed genes in each tissue

```bash
cat analyses/tsv.totalRNASeq.ids.txt |\
while read filename tissue; do
  sort -k2,2gr data/tsv.files/"$filename".tsv |\
  head -1000 |\
  cut -f1 > analyses/"$tissue".1000.most.expressed.genes.txt
done
```

### 3.5 Select the 1,000 least expressed genes in each tissue

```bash
cat analyses/tsv.totalRNASeq.ids.txt |\
while read filename tissue; do 
  sort -k2,2gr data/tsv.files/"$filename".tsv |\
  tail -1000 |\
  cut -f1 > analyses/"$tissue".1000.least.expressed.genes.txt
done
```

This reads the file IDs and, for each row with columns `filename` and `tissue`, sorts according to column 2, which contains TPM, in reverse order from high to low. It then gets the top or bottom 1,000 rows, cuts the first column containing the gene ID, and saves it to a file.

### 3.6 Prepare BED files for these gene sets

For these genes, we need to prepare the BED file just like before: chromosome, start, end, gene ID, strand, coordinate type, and gene ID.

```bash
for tissue in stomach sigmoid_colon; do
  ../bin/selectRows.sh analyses/"$tissue".1000.least.expressed.genes.txt <(awk 'BEGIN{FS=OFS="\t"}{print $4, $0}' annotation/gencode.v24.protein.coding.gene.body.bed) |\
  cut -f2- > annotation/"$tissue".1000.least.expressed.genes.bed
done
```

Here, we take the GENCODE BED file, place the gene ID in the first column, and select those IDs that are present in the top 1,000 least expressed genes. Then, we select all columns from column 2 onward and save the file.

In other words, IDs in the second file that are present in the first file are returned as a subset.

Do the same for the 1,000 most expressed genes:

```bash
for tissue in stomach sigmoid_colon; do
  ../bin/selectRows.sh analyses/"$tissue".1000.most.expressed.genes.txt <(awk 'BEGIN{FS=OFS="\t"}{print $4, $0}' annotation/gencode.v24.protein.coding.gene.body.bed) |\
  cut -f2- > annotation/"$tissue".1000.most.expressed.genes.bed
done
```

## 4. Compute the aggregated signal over the TSS of the selected genes

Create an `aggregation.plot` directory inside `analyses` to store the following results.

```bash
mkdir analyses/aggregation.plot
```

Then, use the `aggregate` function from the `bwtool` package. This needs to be installed correctly on the computer.

Use it on the 1,000 highest expressed genes in stomach, as well as the 1,000 least expressed genes:

```bash
bwtool aggregate 2000:2000 -starts -keep-bed annotation/stomach.1000.most.expressed.genes.bed data/bigWig.files/ENCFF391KDD.bigWig analyses/aggregation.plot/stomach.1000.most.expressed.genes.aggregate.tsv
```

```bash
bwtool aggregate 2000:2000 -starts -keep-bed annotation/stomach.1000.least.expressed.genes.bed data/bigWig.files/ENCFF391KDD.bigWig analyses/aggregation.plot/stomach.1000.least.expressed.genes.aggregate.tsv
```

And the same for sigmoid colon:

```bash
bwtool aggregate 2000:2000 -starts -keep-bed annotation/sigmoid_colon.1000.most.expressed.genes.bed data/bigWig.files/ENCFF886LUE.bigWig analyses/aggregation.plot/sigmoid_colon.1000.most.expressed.genes.aggregate.tsv
```

```bash
bwtool aggregate 2000:2000 -starts -keep-bed annotation/sigmoid_colon.1000.least.expressed.genes.bed data/bigWig.files/ENCFF886LUE.bigWig analyses/aggregation.plot/sigmoid_colon.1000.least.expressed.genes.aggregate.tsv
```

Then, we will use the R script to use the aggregated H3K4me3 signal computed over the set of most and least expressed genes.

### 4.1 Summary of the aggregation plot workflow

First, we downloaded total RNA-seq data for both tissues, and also the GENCODE annotation in the correct assembly and version. From the GTF file, we created a BED file with the correct format and selected the genes that are protein-coding.

Second, for each tissue, we filtered the total RNA-seq matrices by genes that are protein-coding, and selected the gene ID and TPM. Then, we selected the top 1,000 most expressed and least expressed genes, and created a `.txt` file with their IDs. From this, we created a BED file for their genomic locations using the GENCODE file and filtering the matching IDs.

Finally, we used the `aggregate` function from `bwtool` to aggregate the FC signal around the TSS. It works like this: it takes the start coordinate of all genes in the BED file for the 1,000 most expressed genes, looks across 2,000 bp before and after the start coordinate of all genes, and averages the signal at each position relative to the TSS across all genes.

## 5. Compute the correlation coefficient between expression and H3K4me3 levels

From this link, we download the matrix of H3K4me3 FC signals over promoter regions:

```text
https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/H3K4me3.matrix.tsv
```

We will use `wget` to save this into the `analyses` directory.

```bash
wget --no-check-certificate -P analyses "https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/H3K4me3.matrix.tsv"
```

Unlike gene expression, there is no clear consensus on how to assign the value of a histone mark to a gene. In this case, we are provided with a matrix built for our set of protein-coding genes. Specifically, the matrix reports, for each gene, the mean H3K4me3 FC signal around the TSS, ±2 kb. For genes having multiple TSSs, the TSS with the highest value was selected.

Steps:

1. Obtain an expression matrix with the same structure as the one previously used for TPM.
2. Calculate Pearson's and Spearman's correlation coefficients between expression and H3K4me3 levels measured around the TSS of each gene.
3. Visualize the results with a scatterplot.

### 5.1 Prepare the expression matrix

First, verify that the order of the gene IDs is the same in the expression matrices of the two tissues.

```bash
diff <(cut -f1 data/tsv.files/ENCFF268RWA.tsv) <(cut -f1 data/tsv.files/ENCFF918KPC.tsv)
```

Result: OK.

The H3K4me3 matrix has an ID column and then the expression values for sigmoid colon and stomach. We need to extract the TPM column for each tissue and then concatenate them with the same structure.

Take the sigmoid colon matrix `ENCFF268RWA` containing ID and TPM, and paste column 2 from the stomach matrix `ENCFF918KPC`. Then, print the first row as headers.

```bash
paste data/tsv.files/ENCFF268RWA.tsv <(cut -f2 data/tsv.files/ENCFF918KPC.tsv) |\
awk 'BEGIN{FS=OFS="\t"; print "sigmoid_colon", "stomach"}{print}' > analyses/expression.matrix.tsv
```

Now, verify that the row order of the H3K4me3 matrix is the same as the RNA-seq expression matrix.

```bash
diff <(cut -f1 analyses/H3K4me3.matrix.tsv) <(cut -f1 analyses/expression.matrix.tsv)
```

### 5.2 Run the correlation scatterplot script

For this, we will run an R script stored in the `bin` directory: `scatterplot.correlation.R`.

Arguments:

- `--expression`: the RNA expression matrix.
- `--mark`: the histone mark matrix.
- `--tissue`: the column name for the tissue, identifying which column to get from the expression matrix.

```bash
for tissue in sigmoid_colon stomach; do
  Rscript ../bin/scatterplot.correlation.R --expression analyses/expression.matrix.tsv --mark analyses/H3K4me3.matrix.tsv --tissue "$tissue" --output analyses/scatterplot.correlation/scatterplot.correlation."$tissue".pdf
done
```

## 6. Retrieve genes with tissue-specific marking

Here, we will explore the extent of H3K4me3 tissue-specific marking: which genes are marked in only one of the tissues.

We will define four groups of genes:

1. Genes marked by H3K4me3 in both tissues.
2. Genes with sigmoid colon-specific marking.
3. Genes with stomach-specific marking.
4. Genes not marked in any of the tissues.

To do so, we will use the `intersect` command from `bedtools`.

Then, we will perform a functional analysis on these sets of genes by looking at their GO term enrichment using Metascape. Finally, we will compare the distribution of expression values between the four sets of genes.

We will store this analysis in the `peaks.analysis` directory. For peaks, we will work with `bigBed` files.

### 6.1 Genes with H3K4me3 peaks in each tissue

First, convert the `bigBed` files for each tissue into BED files, using the function `bigBedToBed`.

Careful: I had some environment conflicts with the UCSC tools, so I had to create a specific environment to work with this function. It uses the `bigBed` accession IDs for both tissues and saves them in BED files.

```bash
mkdir data/bed.files
```

```bash
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

Then, download the list of promoters and store it in the `annotation` directory.

```bash
wget --no-check-certificate -P annotation "https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/gencode.v24.protein.coding.non.redundant.TSS.bed"
```

We then use this list of genes in BED format as a base to select those for which there is a peak in each ChIP-seq experiment.

```bash
cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u |\
  cut -f7 |\
  sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".H3K4me3.txt
done
```

Explanation:

- From the `bigBed` index file, we select the first two columns: ID and tissue.
- Then, using the `while` clause for the two columns, `filename` and `tissue`, we run `bedtools intersect` on the list of promoters.
- `bedtools intersect` reports overlaps between two feature files, specified with `-a` and `-b`, in BED, GFF, VCF, or BAM format.
- The first file is the GENCODE BED file for non-redundant TSSs of protein-coding genes, containing start-end coordinates.
- The second file is the BED file containing ChIP-seq peaks for the histone mark in different genomic regions.
- It checks, for each TSS region in `-a`, whether this interval overlaps any peak in `-b`.
- If yes, it reports the particular genomic region from `-a`.
- The `-u` flag indicates that this match is reported only once, in case the same region overlaps one or more peaks. For example, this could happen if the TSS region is wide enough to have several peaks reported.
- We get the gene ID from the seventh column of the BED file, with the format: chromosome, start, end, transcript ID, score, strand, gene ID.
- The obtained list has the coordinates of TSS regions where there is at least one H3K4me3 peak.
- These relate to transcripts. The same gene can have several transcripts, each with its own TSS. In the end, we can end up with a list with more than one TSS with a histone mark peak.
- The final instruction, `sort -u`, keeps only one gene ID.

How many genes have H3K4me3 peaks in each tissue?

```bash
ls analyses/peaks.analysis/ |\
while read filename; do
  wc -l analyses/peaks.analysis/"$filename"
done
```

```text
15032 analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt
15557 analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt
```

Questions about `bedtools intersect`:

- Why are we using the flag `-u`?
- Why are we using the command `sort -u` then?

These are already answered in the explanation above.

### 6.2 Genes marked in both tissues

For this task, we will use the `selectRows.sh` function to select which gene IDs from the gene list with peaks in sigmoid colon are also present in stomach.

```bash
../bin/selectRows.sh analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt |\
cut -d "." -f1 > analyses/peaks.analysis/gene.marked.both.tissues.H3K4me3.txt
```

### 6.3 Sigmoid colon-specific marks

For this, we will use another custom `.sh` function, `discardRows.sh`. It selects rows from file 2 that are not present in file 1.

```bash
../bin/discardRows.sh analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt |\
cut -d "." -f1 > analyses/peaks.analysis/genes.with.sigmoid_colon.specific.peaks.H3K4me3.txt
```

### 6.4 Stomach-specific marks

```bash
../bin/discardRows.sh analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt |\
cut -d "." -f1 > analyses/peaks.analysis/genes.with.stomach.specific.peaks.H3K4me3.txt
```

### 6.5 Genes not marked in any tissue

Use `discardRows.sh` to select those genes in GENCODE that do not appear in any of the three lists: both tissues, stomach-specific, or sigmoid colon-specific. We use `cat` to concatenate all three lists.

```bash
../bin/discardRows.sh <(cat analyses/peaks.analysis/genes.marked.both.tissues.H3K4me3.txt analyses/peaks.analysis/genes.with.stomach.specific.peaks.H3K4me3.txt analyses/peaks.analysis/genes.with.sigmoid_colon.specific.peaks.H3K4me3.txt) <(cut -f7 annotation/gencode.v24.protein.coding.gene.body.bed |\
cut -d "." -f1) > analyses/peaks.analysis/genes.not.marked.H3K4me3.txt
```

### 6.6 Represent expression values between the four gene sets

We will use the R script `boxplot.expression.R`, which needs the expression matrix from both tissues and the lists of genes from each group.

```bash
Rscript ../bin/boxplot.expression.R --expression analyses/expression.matrix.tsv --marked_both_tissues analyses/peaks.analysis/genes.marked.both.tissues.H3K4me3.txt --stomach_specific analyses/peaks.analysis/genes.with.stomach.specific.peaks.H3K4me3.txt --sigmoid_colon_specific analyses/peaks.analysis/genes.with.sigmoid_colon.specific.peaks.H3K4me3.txt --not_marked analyses/peaks.analysis/genes.not.marked.H3K4me3.txt --output analyses/peaks.analysis/boxplot.expression.pdf
```

### 6.7 GO enrichment analysis

Use the set of protein-coding genes as the universe for the GO enrichment analysis.

Extract the gene ID from the seventh column of the GENCODE protein-coding gene body coordinates and remove everything after the `.`.

```bash
cut -f7 annotation/gencode.v24.protein.coding.gene.body.bed |\
cut -d "." -f1 > analyses/peaks.analysis/universe.genes.txt
```

## 7. Compute the percentage of genes with peaks of H3K4me3 and POLR2A

We will retrieve peak calls targeting POLR2A for the same donor and tissues. We will then compute the percentage of genes with peaks of H3K4me3 and POLR2A in the two tissues.

### 7.1 Retrieve peak calls in BED format for POLR2A

```bash
grep -F POLR2A-human metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_IDR_thresholded_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u >> analyses/bigBed.peaks.ids.txt
```

We are appending the new searches to the existing file.

Download the `bigBed` files for the corresponding IDs, only the new ones from POLR2A.

```bash
awk '$3=="POLR2A-human"{print $1}' analyses/bigBed.peaks.ids.txt |\
while read filename; do 
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

### 7.2 Check MD5 integrity

Retrieve MD5 hashes of the files from the metadata.

```bash
../bin/selectRows.sh <(awk '$3=="POLR2A-human"{print $1}' analyses/bigBed.peaks.ids.txt) metadata.tsv | cut -f1,46 > data/bigBed.files/tmp
```

Compute MD5 hashes on the downloaded files.

```bash
cat data/bigBed.files/tmp |\
while read filename original_md5sum; do 
  md5sum data/bigBed.files/"$filename".bigBed |\
  awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}'
done >> data/bigBed.files/md5sum.txt
rm data/bigBed.files/tmp
```

Make sure there are no files for which the original and computed MD5 hashes differ.

```bash
awk '$2!=$3' data/bigBed.files/md5sum.txt
```

Convert the `bigBed` files to BED files using the UCSC tools.

```bash
awk '$3=="POLR2A-human"{print $1}' analyses/bigBed.peaks.ids.txt |\
while read filename; do 
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

### 7.3 Genes with POLR2A peaks in each tissue

Now we have a BED file with regions with peaks associated with POLR2A. Just like we did before, we want to know how many of these peaks are located in TSS positions of protein-coding genes, stored in the `gencode.v24.protein.coding.non.redundant.TSS.bed` file.

```bash
grep -F POLR2A analyses/bigBed.peaks.ids.txt |\
cut -f-2 |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u |\
  cut -f7 |\
  sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".POLR2A.txt
done
```

### 7.4 Venn diagram

Use a custom R function, `VennDiagram.4groups.R`, which takes the list of genes with peaks of POLR2A and H3K4me3 in stomach or sigmoid colon, and calculates the overlap.

```bash
Rscript ../bin/VennDiagram.4groups.R --setA analyses/peaks.analysis/genes.with.peaks.stomach.H3K4me3.txt --setB analyses/peaks.analysis/genes.with.peaks.stomach.POLR2A.txt --setC analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.H3K4me3.txt --setD analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.POLR2A.txt --output analyses/peaks.analysis/Venn.Diagram.H3K4me3.POLR2A.png
```
