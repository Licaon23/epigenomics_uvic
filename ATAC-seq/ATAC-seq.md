# EN‑TEx ATAC‑seq data: downstream analyses

Downstream analysis of ATAC‑seq peaks for **stomach** and **sigmoid colon** from
the same EN‑TEx donor. For each tissue we retrieve the pseudoreplicated peak
calls from ENCODE, convert them to BED, and use `bedtools intersect` to count
how many peaks fall in promoter regions and how many fall outside gene bodies.

---

## 1. Set up

Enter the Docker container and move to the ATAC‑seq directory.

```bash
# Enter the docker
docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course

# And go to the ATAC-seq directory
```

Create the directory structure to store data, analyses and annotation files.

```bash
mkdir -p data/bigBed.files data/bed.files
mkdir analyses
mkdir annotation
```

---

## 2. Download the metadata

The metadata was obtained from **ENCODE -  functional genomics experiments** for
the individual donor `ENCDO451RUA`, filtering by **Assay type: ATAC-seq** and
**Biosample: stomach & sigmoid colon**, and downloading the selected files.

Reference report URL (ENCODE portal):

```
https://www.encodeproject.org/report/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&limit=all&assay_title=ATAC-seq&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon
```

Download the metadata file:

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon&type=Experiment"
```

---

## 3. Retrieve the ATAC-seq peak files

From the metadata the ATAC-seq peaks (bigBed narrow, pseudoreplicated
peaks) were retrieved for assembly **GRCh38**, for stomach and sigmoid colon from the same
donor, making sure the md5sum values coincide with the ones provided by ENCODE.

First, inspect the header to identify the relevant fields:

```bash
head -n 1 metadata.tsv | awk -F'\t' '{for (i=1; i<=NF; i++) print i, $i}'
```

Fields 1, 8, 11 and 4 provide the information on file accession, assay,
biosample, file format and output type. Parse the **bigBed narrowPeak**
files for **pseudoreplicated peaks** in the **GRCh38** assembly.

```bash
grep -F ATAC-seq metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $4, $6, $46}' > analyses/bigBed.peaks.ids.txt
```

Then the file accession (first column) is selected and loop over the list to
download each file from ENCODE using `wget`, storing it in `data/bigBed.files`.

```bash
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

---

## 4. Verify the md5sum of the downloaded files

The file accession of each bigBed file is used to compute the md5sum and compare
it with the original one, which was already retrieved previously in the same id
file. Then it was check whether the original and computed values differ.

```bash
awk 'BEGIN{FS=OFS="\t"} {print $1, $6}' analyses/bigBed.peaks.ids.txt | \
while read filename original_md5sum; do
  computed_md5sum=$(md5sum data/bigBed.files/"$filename".bigBed | awk '{print $1}')
  printf "%s\t%s\t%s\n" "$filename" "$original_md5sum" "$computed_md5sum"
done > data/bigBed.files/md5sum.txt

# Print any file whose original and computed md5sum differ (empty output = all OK)
awk '$2 != $3' data/bigBed.files/md5sum.txt
```

An empty output from the last command means every file matches its expected
checksum.

---

## 5. Intersection analysis: promoters and gene bodies

For each tissue an intersection analysis was run with **BEDTools** and report:

1. the number of peaks that intersect the **promoter regions**, and
2. the number of peaks that fall **outside gene coordinates** (the whole gene
   body, not just the promoter regions).

### 5.1 Convert bigBed to BED

```bash
cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

### 5.2 Download the annotation files

Download the annotations needed: (1) the **GRCh38 v24 GENCODE** assembly, to
get the coordinates of the gene bodies, and (2) a list of **promoters (�-2 kb to
+2 kb)** around the TSS of protein-coding genes.

```bash
wget -P annotation "https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz"
wget --no-check-certificate -P annotation "https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/gencode.v24.protein.coding.non.redundant.TSS.bed"
```

### 5.3 Prepare the gene-body BED file for protein-coding genes

```bash
gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz

awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed
```

Step by step, this pipeline:

- gets the third column and checks whether it contains `gene`;
- filters rows containing `protein_coding` (fixed string);
- the last column of attributes is separated by `;`.The 'cut' instruction establishes `;` as
  the delimiter and keeps the first field, as if it were all in the same column  with tab delimiters;
- then takes these fields and prints the columns of interest;
- `sed` uses the substitution to remove quotes from the annotations;
- selects the rows that do not contain `chrM` in the chromosome column (1)
  and transforms the coordinate system from 1-based (GTF) to 0-based (BED):
  1-based puts the coordinate exactly on the bp, whereas 0-based sits in between,
  so we subtract 1 from the start. This returns every row not containing `chrM`,
  mutating the start column by subtracting 1 and returning the whole line;
- pipes the output to a new file in BED format.

### 5.4 Run the intersection analysis

This list of genes (BED file) was used to compare how many peaks overlap the
promoter regions, or fall outside the gene body.

- From the bigBed index file, the two first columns are selected (`-f-2`): ID and
  tissue.
- Then, using the `while` clause for the two columns (filename, tissue):
  - run the `bedtools intersect` command;
  - the first file (`-a`) is the BED file for ATAC-seq peaks; the second (`-b`)
    is the GENCODE BED file for non-redundant TSS of protein-coding genes
    (start-end coordinates), or the gene-body coordinates. It checks, for each
    peak delimited by start and end coordinates (`-a`), whether this interval
    overlaps a promoter region or a gene body;
  - if yes, it reports that particular peak from `-a`;
  - the `-u` flag indicates that we report each match only once, in case the same
    peak overlaps one or more promoter regions, or matches more than one gene
    body;
  - in the "outside gene bodies" analysis we use the `-v` flag to report those
    peaks (`-a`) that have **no** overlap with any gene‑body interval (`-b`).

```bash
echo -e "tissue\tanalysis\tn_peaks" > analyses/peaks.intersections.summary.tsv

cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do

  promoter_out="analyses/peaks.intersect.promoter.${tissue}.bed"
  outside_gene_out="analyses/peaks.outside.geneBodyCoordinates.${tissue}.bed"

  bedtools intersect \
    -a data/bed.files/"$filename".bed \
    -b annotation/gencode.v24.protein.coding.non.redundant.TSS.bed \
    -u > "$promoter_out"

  bedtools intersect \
    -a data/bed.files/"$filename".bed \
    -b annotation/gencode.v24.protein.coding.gene.body.bed \
    -v > "$outside_gene_out"

  n_promoter=$(wc -l < "$promoter_out")
  n_outside_gene=$(wc -l < "$outside_gene_out")

  echo -e "${tissue}\tpeaks_intersecting_promoters\t${n_promoter}" >> analyses/peaks.intersections.summary.tsv
  echo -e "${tissue}\tpeaks_outside_gene_body_coordinates\t${n_outside_gene}" >> analyses/peaks.intersections.summary.tsv

done
```

Report the counts for each intersection analysis:

```bash
column -t -s $'\t' analyses/peaks.intersections.summary.tsv
```

---

## 6. Results

| tissue        | analysis                            | n_peaks |
|---------------|-------------------------------------|--------:|
| sigmoid_colon | peaks_intersecting_promoters        |   47871 |
| sigmoid_colon | peaks_outside_gene_body_coordinates |   37035 |
| stomach       | peaks_intersecting_promoters        |   44749 |
| stomach       | peaks_outside_gene_body_coordinates |   34537 |


