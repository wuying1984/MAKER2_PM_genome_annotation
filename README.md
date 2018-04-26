# MAKER2_PM_genome_annotation
## Using MAKER2 FOR ANNOTATION OF FOUR DICOT PM genomes
### Reference to  https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
###             & https://github.com/CompSynBioLab-KoreaUniv/FunGAP#step1

### Software & Data
### Software prerequisites:
###### RepeatModeler(1.0.4) and RepeatMasker (4.0.5) with all dependencies (I used NCBI BLAST) and RepBase (version used was 20150807).
###### MAKER version 2.31.9 (though any other version 2 releases should be okay).
###### Augustus version 3.3.
###### BUSCO version 3.
###### SNAP https://github.com/KorfLab/SNAP
###### BEDtools version 2.24.0

### Raw data/resources:
#### *1: Genome scaffolds*:
##### >Genome.fa: The de novo assembled reference genome using CLCbio. 
```
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
```

#### *2: EST data*:
##### Transcriptome.fasta from de novo transcriptome assembly using Trinity
##### Reads: mycelia + haustoria RNA-seq reads 
##### *First*, map the haustoria RNAseq reads to the contig using Tophat
##### *Secondly*, do transcriptome assembly
```
Trinity --seqType fq --left R1-common.fastq.gz --right R2-common.fastq.gz --jaccard_clip --max_memory 100G --CPU 24 --output trinity_out
```
###### Use --jaccard_clip for Trinity because high gene density leads to UTR overlap in the assembly. This option helps avoid fusion of neighbor transcripts. 
###### It may make sense to do some post-processing of this assembly. However, I did not. Because powdery mildew transcripts should be quite different from its host.

#### *3: Full protein set*:
##### Complete UniProtKB/Swiss-Prot data set in FASTA format: `ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`
##### Bgh proteins:Blumeria graminis f. sp. hordei DH14 (6495)
`ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/151/065/GCA_000151065.3_ASM15106v3/GCA_000151065.3_ASM15106v3_protein.faa.gz`
##### Bgt Proteins: Blumeria graminis f. sp. tritici 96224 . (6525)
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/418/435/GCA_000418435.1_Bgt_454_newbler_assembly/GCA_000418435.1_Bgt_454_newbler_assembly_protein.faa.gz
##### Ene Proteins: Erysiphe necator C strain (6484)
`ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/798/715/GCA_000798715.1_ASM79871v1/GCA_000798715.1_ASM79871v1_protein.faa.gz`
##### And also Gor CSEP

#### *4.Repeat Annotation*:
##### De Novo Repeat Identification
###### For genome annotation, it is very important to identify repetitive content. Sometimes, we can download existing libraries from Repbase or from other efforts. However,it is also important to identify repeats from de novo assembly using RepeatModeler. 
```
  BuildDatabase -name UCSC1 -engine ncbi ../../UCSC1_CLC_de_novo_rmhost_mod.fa
  RepeatModeler -engine ncbi -pa 8 -database UCSC1 1>UCSC1_repeatmodeler.o 2>UCSC1_repeatmodeler.e
``` 
  
#### *5. Run BUSCO*:
```
run_BUSCO.py -i UCSC1_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UCSC1_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea  >UCSC1_BUSCO_so_long.out&
run_BUSCO.py -i UMSG1_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UMSG1_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea  >UMSG1_BUSCO_so_long.out&
run_BUSCO.py -i UMSG2_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UMSG2_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea  >UMSG2_BUSCO_so_long.out&
run_BUSCO.py -i UMSG3_CLC_de_novo_rmhost_mod.fa -l ~/program/BUSCO/sordariomyceta_odb9 -o UMSG3_BUSCO_so_long -m geno -c 1 -sp botrytis_cinerea  >UMSG3_BUSCO_so_long.out&
```
##### I added --long option to generate the augustus training models.
##### However, if --long option is used, the genome completeness will goes down from 84% to 85%. I do not know why.
##### But if I used trained models to perform estimation, the estimated completeness went back to 84%-85% again.

#### *6. Do three iterative round of Maker*
#### *First_round*:
#### maker_opts.ctl files:
```
est=Trinity.fasta
protein=unipro_sport_add4genome.fasta
rmlib=repeat.consensi.fa #############very important
repeat_protein=repeat.consensi.fa
softmask=1
snaphmm=snap.hmm ##########training from CEGMA gff file
gmhmm=gmhmm.mod ###########training from genome sequence
augustus=species ##########model derived from BUSCO analysis (using --long option)
est2genome=1
protein2genome=1
max_dna_len=100000
min_contig=500
AED_threshold=1
min_protein=30
always_complete=1
split_hit=3000 #######intron size limitation
single_exon=1 ########turn it on for fungi genome annotation
single_length=250 ####single exon length 
correct_est_fusion=0 #Did not turn one in the first round (refer to https://groups.google.com/forum/#!topic/maker-devel/J_ZLTFQ3xN4)
```

#### *Second round*:
#### retraining snap and augustus
#### 1)SNAP
`gff3_merge -d *index.log`
##### export 'confident' gene models from MAKER and rename to something meaningful
`maker2zff -d ../../Bcon_rnd1.maker.output/Bcon_rnd1_master_datastore_index.log`
##### gather some stats and validate
`fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna  -validate > validate.log 2>&1`
##### collect the training sequences and annotations, plus 1000 surrounding bp for training
`fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1`
##### create the training parameters
```
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl Bcon_rnd1.zff.length50_aed0.25 params > Bcon_rnd1.zff.length50_aed0.25.hmm
```

#### 2)maker_opts.ctl files:
```
est=Trinity.fasta
protein=unipro_sport_add4genome.fasta
rmlib=repeat.consensi.fa #############very important
repeat_protein=repeat.consensi.fa
softmask=1
snaphmm=snap.hmm ##########training from CEGMA gff file
gmhmm=gmhmm.mod ###########training from genome sequence
augustus=species ##########model derived from BUSCO analysis (using --long option)
est2genome=0 ##############set to 0 from second round
protein2genome=0 ##############set to 0 from second round
max_dna_len=100000
min_contig=500
AED_threshold=1
min_protein=30
always_complete=1
split_hit=5000 #######intron size limitation
single_exon=1 ########turn it on for fungi genome annotation
single_length=250 ####single exon length 
correct_est_fusion=1 #turn it on from the second round
```
#### *Third round*:
```
est=Trinity.fasta
protein=unipro_sport_add4genome.fasta
rmlib=repeat.consensi.fa #############very important
repeat_protein=repeat.consensi.fa
softmask=1
snaphmm=snap.hmm ##########training from CEGMA gff file
gmhmm=gmhmm.mod ###########training from genome sequence
augustus=species ##########model derived from BUSCO analysis (using --long option)
est2genome=0 ##############set to 0 from second round
protein2genome=0 ##############set to 0 from second round
max_dna_len=100000
min_contig=500
AED_threshold=1
min_protein=30
always_complete=1
split_hit=5000 #######intron size limitation
single_exon=1 ########turn it on for fungi genome annotation
single_length=250 ####single exon length 
correct_est_fusion=1 #turn it on from the second round
```




