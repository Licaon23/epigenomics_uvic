# Distal regulatory activity

The ATAC-seq peaks lying outside gene coordinates in stomach and sigmoid colon,
obtained in the previous section, are used here as a starting point to build a
catalogue of distal regulatory regions. For each tissue, the candidate distal
regulatory elements are defined, their distance to the closest gene is computed,
and the mean and median of these distances are reported.

All commands are run inside the course Docker container:

```bash
docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

---

## Task 1 - Folder structure

A folder `regulatory_elements` is created inside `epigenomics_uvic` to store all
subsequent results, together with the directory structure for the bigBed and BED
data files, the analyses and the annotation.

```bash
mkdir regulatory_elements
cd regulatory_elements
mkdir -p data/bigBed.files data/bed.files
mkdir analyses
mkdir annotation
```

---

## Tasks 2 and 3 - Candidate distal regulatory elements and their start coordinates

Distal regulatory regions are usually found to be flanked by both H3K27ac and
H3K4me1. From the starting catalogue of open regions in each tissue, those that
overlap peaks of H3K27ac **AND** H3K4me1 in the corresponding tissue are
selected, producing a list of candidate distal regulatory elements per tissue
(Task 2). Then, the regulatory elements located on chromosome 1 are kept, and a
file with the name of the regulatory region (the name of the original ATAC-seq
peak) and its start (5') coordinate is generated (Task 3).

### Retrieve the ChIP-seq metadata

First, the metadata for the ChIP-seq experiments of the histone markers H3K27ac
and H3K4me1 is retrieved for the same subject in ENCODE, in stomach and sigmoid
colon. The URL was obtained by filtering through these criteria.

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_slims=DNA+binding&target.label=H3K27ac&target.label=H3K4me1&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&type=Experiment"
```

### Select and download the peak files

The ChIP-seq peaks for these two markers (bigBed narrow, pseudoreplicated peaks,
GRCh38 assembly) are retrieved from the metadata for stomach and sigmoid colon.
The file accession (first column) is selected and the list is looped over to
download each file from ENCODE with `wget`, storing it in `data/bigBed.files`.

```bash
grep -E "(H3K27ac|H3K4me1)" metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23, $46, $55}' |\
sort -k2,2 -k3,3 -k5,5r |\
sort -k2,2 -k3,3 -u > analyses/bigBed.ChIPseq.peaks.ids.txt
```

```bash
cut -f1 analyses/bigBed.ChIPseq.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

### Verify the md5sum of the downloaded files

The integrity of the files is checked by computing the md5sum of each bigBed file
and comparing it with the original one stored in the id file. Any row printed by
the last command indicates a mismatch (empty output means all files are correct).

```bash
awk 'BEGIN{FS=OFS="\t"} {print $1, $4}' analyses/bigBed.ChIPseq.peaks.ids.txt | \
while read filename original_md5sum; do
  computed_md5sum=$(md5sum data/bigBed.files/"$filename".bigBed | awk '{print $1}')
  printf "%s\t%s\t%s\n" "$filename" "$original_md5sum" "$computed_md5sum"
done > data/bigBed.files/md5sum.txt

awk '$2 != $3' data/bigBed.files/md5sum.txt
```

### Convert bigBed files to BED

```bash
cut -f1 analyses/bigBed.ChIPseq.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

### Intersect the ATAC-seq peaks with the two ChIP-seq markers

It is determined which ATAC-seq peaks lying outside gene body coordinates overlap
ChIP-seq peaks of H3K27ac **AND** H3K4me1 in each tissue. For each tissue, the
ATAC-seq peaks that overlap each ChIP-seq marker are first obtained separately,
keeping each overlapping ATAC-seq peak.

```bash
cut -f-3 analyses/bigBed.ChIPseq.peaks.ids.txt |\
while read filename tissue experiment; do

  bedtools intersect \
    -a ../ATAC-seq/analyses/peaks.outside.geneBodyCoordinates."$tissue".bed \
    -b data/bed.files/"$filename".bed \
    -u > analyses/peaks.outside.geneBodyCoordinates."$tissue"."$experiment".bed
done
```

Then, the ATAC-seq peaks present in both intersections are kept: the two tissue
names are obtained and the previous peak files are intersected. The number of
rows in each final list is counted, giving the number of regulatory elements
obtained per tissue. In addition, for Task 3, a subset of the peaks located on
chromosome 1 is made, retrieving the peak name (column 4) and the start 5'
coordinate (column 2).

```bash
echo -e "Element\ttissue\tCount" > analyses/summaryRegElements.tsv
cut -f2 analyses/bigBed.ChIPseq.peaks.ids.txt |\
sort -u |\
while read tissue; do

  bedtools intersect \
    -a analyses/peaks.outside.geneBodyCoordinates."$tissue".H3K27ac-human.bed \
    -b analyses/peaks.outside.geneBodyCoordinates."$tissue".H3K4me1-human.bed \
    -u > analyses/peaks.outside.geneBodyCoordinates."$tissue".distalRegulatoryRegions.bed

  # Count number of regulatory elements per tissue
  count=$(cat analyses/peaks.outside.geneBodyCoordinates."$tissue".distalRegulatoryRegions.bed | wc -l)
  echo -e "Regulatory elements\t$tissue\t$count" >> analyses/summaryRegElements.tsv

  # Focus on chr1 elements and subset the name of the peak and the 5' start coordinate
  cat analyses/peaks.outside.geneBodyCoordinates."$tissue".distalRegulatoryRegions.bed |\
  awk 'BEGIN{FS=OFS="\t"}$1=="chr1"{print $4, $2}' >> analyses/"$tissue".regulatory.elements.starts.tsv

done

column -t -s $'\t' analyses/summaryRegElements.tsv
```

The number of candidate distal regulatory elements per tissue:

| Element             | tissue        | Count |
|---------------------|---------------|------:|
| Regulatory elements | sigmoid_colon | 14215 |
| Regulatory elements | stomach       |  8454 |

---

## Task 4 - Gene start coordinates on chromosome 1

The protein-coding genes located on chromosome 1 are considered. From the BED
file of gene body coordinates generated during the ChIP-seq hands-on, a
tab-separated file `gene.starts.tsv` is prepared, storing the name of the gene in
the first column and the start coordinate of the gene in the second column. For
genes located on the minus strand, the start coordinate is taken at the 3'.

```bash
cat ../ChIP-seq/annotation/gencode.v24.protein.coding.gene.body.bed |\
awk 'BEGIN{FS=OFS="\t"} $1=="chr1" {if ($6=="+"){start=$2} else {start=$3}; print $4, start}' > annotation/gene.starts.tsv
```

---

## Task 5 - Complete the distance script

The `get.distance.py` script is copied into the `epigenomics_uvic/bin` folder.
For a given coordinate `--start`, the script is completed so that it returns the
closest gene, the start of the gene and the distance of the regulatory element.
The core of the completion iterates over every gene start, keeps the smallest
absolute distance to the regulatory element, and records the corresponding gene.

```python
for line in open_input.readlines():
    gene, y = line.strip().split('\t')
    position = int(y)
    distance = abs(position - enhancer_start)

    if distance < x:
        x = distance
        selectedGene = gene
        selectedGeneStart = position
```

The script is checked with the reference example:

```bash
../bin/get.distance.py --input annotation/gene.starts.tsv --start 980000
```

which returns the expected result:

```
ENSG00000187642.9	982093	2093
```

---

## Task 6 - Closest gene and distance for each regulatory element

For each regulatory element contained in the file
`regulatory.elements.starts.tsv`, the closest gene and the distance to it are
retrieved using the python script completed above. For each tissue, an output
file is initialised with a header, and the regulatory-element list is looped
over: the distance to the nearest gene is computed by the script and the
regulatory element name is pasted together with the script output.

```bash
cut -f2 analyses/bigBed.ChIPseq.peaks.ids.txt |\
sort -u |\
while read tissue; do
    rm analyses/"$tissue".regulatoryElements.genes.distances.tsv
    touch analyses/"$tissue".regulatoryElements.genes.distances.tsv
    echo -e "RegulatoryElementName\tclosestGeneName\tStartCoordGene\tDistance" >> analyses/"$tissue".regulatoryElements.genes.distances.tsv

    cat analyses/"$tissue".regulatory.elements.starts.tsv |\
    while read element start; do
        # Calculate distance from nearest gene in list and save results in out
        out=$(python ../bin/get.distance.py --input annotation/gene.starts.tsv --start "$start")
        # Paste regulatory element name with output
        echo -e "$element\t$out" >> analyses/"$tissue".regulatoryElements.genes.distances.tsv

    done

done
```

---

## Task 7 - Mean and median distance per tissue

R is used to compute the mean and the median of the distances stored in
`regulatoryElements.genes.distances.tsv` for each tissue.

```bash
Rscript - << 'end'
    files <- grep('regulatoryElements.genes.distance',list.files('analyses/'),value=TRUE)
    lr <- t(sapply(c('stomach','sigmoid_colon'),function(tissue){
        d <- read.table(paste0('analyses/',files[grep(tissue,files)]),sep='\t',header=TRUE)
        c('Tissue'=tissue,'Median.distance'=median(d[,'Distance'],na.rm=TRUE),'Mean.distance'=round(mean(d[,'Distance'],na.rm=TRUE),2))
    }))
    write.table(lr,paste0('analyses/distance.summary.txt'),sep='\t',row.names=FALSE)
end

column -t -s $'\t' analyses/distance.summary.txt
```

These are the results:

| Tissue        | Median.distance | Mean.distance |
|---------------|----------------:|--------------:|
| stomach       |           27728 |      45370.43 |
| sigmoid_colon |           35802 |      73635.89 |
