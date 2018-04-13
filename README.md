# MAKER2_PM_genome_annotation
# Using MAKER2 FOR ANNOTATION OF FOUR DICOT PM genomes
# Reference to  https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2

Software & Data
Software prerequisites:
RepeatModeler(1.0.4) and RepeatMasker (4.0.5) with all dependencies (I used NCBI BLAST) and RepBase (version used was 20150807).
MAKER version 2.31.9 (though any other version 2 releases should be okay).
Augustus version 3.3.
BUSCO version 3.
SNAP https://github.com/KorfLab/SNAP
BEDtools version 2.24.0

Raw data/resources:
1: Genome scaffolds:
U*.fa: The de novo assembled reference genome using CLCbio. 

2: EST data:
Trinigy.fasta from de novo transcriptome assembly using Trinity
Reads: mycelia RNAseq + haustoria reads (which can be mapped to scaffolds)
It may make sense to do some post-processing of this assembly.

3: Full protein set
Complete UniProtKB/Swiss-Prot data set in FASTA format: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
Bgh proteins:Blumeria graminis f. sp. hordei DH14 (6495)
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/151/065/GCA_000151065.3_ASM15106v3/GCA_000151065.3_ASM15106v3_protein.faa.gz
Bgt Proteins: Blumeria graminis f. sp. tritici 96224 . (6525)
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/418/435/GCA_000418435.1_Bgt_454_newbler_assembly/GCA_000418435.1_Bgt_454_newbler_assembly_protein.faa.gz
Ene Proteins: Erysiphe necator C strain (6484)
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/798/715/GCA_000798715.1_ASM79871v1/GCA_000798715.1_ASM79871v1_protein.faa.gz

#4: Repeat lib
#A curated snake repeat library derived from 14 snake species (from internal efforts).

Repeat Annotation
1. De Novo Repeat Identification
For genome annotation, it is very important to identify repetitive content. Sometimes, we can download existing libraries from Repbase or from other efforts. However,it is also important to identify repeats from de novo assembly using RepeatModeler. 

BuildDatabase -name UCSC1 -engine ncbi ../../UCSC1_CLC_de_novo_rmhost_mod.fa
RepeatModeler -engine ncbi -pa 8 -database UCSC1 1>UCSC1_repeatmodeler.o 2>UCSC1_repeatmodeler.e



