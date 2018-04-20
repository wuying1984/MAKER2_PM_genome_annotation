# MAKER2_PM_genome_annotation
# Using MAKER2 FOR ANNOTATION OF FOUR DICOT PM genomes
# Reference to  https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
#             & https://github.com/CompSynBioLab-KoreaUniv/FunGAP#step1

#Software & Data
#Software prerequisites:
RepeatModeler(1.0.4) and RepeatMasker (4.0.5) with all dependencies (I used NCBI BLAST) and RepBase (version used was 20150807).
MAKER version 2.31.9 (though any other version 2 releases should be okay).
Augustus version 3.3.
BUSCO version 3.
SNAP https://github.com/KorfLab/SNAP
BEDtools version 2.24.0

#Raw data/resources:
#1: Genome scaffolds:
U*.fa: The de novo assembled reference genome using CLCbio. 
          Ambiguous trim = Yes
          Ambiguous limit = 2
          Quality trim = Yes
          Quality limit = 0.05
          Use colorspace = No
          Create report = Yes
          Also search on reversed sequence = Yes
          Save discarded sequences = Yes
          Remove 5' terminal nucleotides = Yes
          Number of 5' terminal nucleotides = 9
          Minimum number of nucleotides in reads = 83
          Discard short reads = Yes
          Remove 3' terminal nucleotides = No
          Trim adapter list = Illumina Trim Adapter List
          Discard long reads = No
          Save broken pairs = No

#2: EST data:
Transcriptome.fasta from de novo transcriptome assembly using Trinity
Reads: mycelia + haustoria RNA-seq reads 
First, map the haustoria RNAseq reads to the contig using Tophat
Secondly, do transcriptome assembly
Use --jaccard_clip for Trinity because high gene density leads to UTR overlap in the assembly. This option helps avoid fusion of neighbor transcripts. 
#reference guided :A BAM-format file for genome-guided assembly is generated by Tophat and Samtools format converter (SAM file to sorted BAM file). 
Trinity --genome_guided_bam <BAM_FILE> --genome_guided_max_intron 2000 --max_memory 10G --CPU <NUMBER_OF_CORES> --output <OUTPUT_DIR>  --jaccard_clip

It may make sense to do some post-processing of this assembly. However, I did not.

#3: Full protein set
Complete UniProtKB/Swiss-Prot data set in FASTA format: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
Bgh proteins:Blumeria graminis f. sp. hordei DH14 (6495)
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/151/065/GCA_000151065.3_ASM15106v3/GCA_000151065.3_ASM15106v3_protein.faa.gz
Bgt Proteins: Blumeria graminis f. sp. tritici 96224 . (6525)
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/418/435/GCA_000418435.1_Bgt_454_newbler_assembly/GCA_000418435.1_Bgt_454_newbler_assembly_protein.faa.gz
Ene Proteins: Erysiphe necator C strain (6484)
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/798/715/GCA_000798715.1_ASM79871v1/GCA_000798715.1_ASM79871v1_protein.faa.gz

#4. Repeat Annotation
De Novo Repeat Identification
For genome annotation, it is very important to identify repetitive content. Sometimes, we can download existing libraries from Repbase or from other efforts. However,it is also important to identify repeats from de novo assembly using RepeatModeler. 
  BuildDatabase -name UCSC1 -engine ncbi ../../UCSC1_CLC_de_novo_rmhost_mod.fa
  RepeatModeler -engine ncbi -pa 8 -database UCSC1 1>UCSC1_repeatmodeler.o 2>UCSC1_repeatmodeler.e
  
  
#5. Run BUSCO
run_BUSCO.py -i UCSC1_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UCSC1_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea --long >UCSC1_BUSCO_so_long.out&
run_BUSCO.py -i UMSG1_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UMSG1_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea --long >UMSG1_BUSCO_so_long.out&
run_BUSCO.py -i UMSG2_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UMSG2_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea --long >UMSG2_BUSCO_so_long.out&
run_BUSCO.py -i UMSG3_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UMSG3_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea --long >UMSG3_BUSCO_so_long.out&

