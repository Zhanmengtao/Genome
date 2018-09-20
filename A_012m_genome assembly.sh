############################################### Prepare data ###################################################################
# Merge subreads of 2 cell
samtools merge 012m.subreads.bam run1/m54171_180501_055307.subreads.bam run2/m54171_180502_184922.subreads.bam -@ 32
samtools fastq -0 012m.subreads.fq -@ 32 subreads.bam

########################################### Data quality control and error correction ##################################
# quality control
mkdir Trimmomatic
cd Trimmomatic
java -jar /opt/biosoft/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 32 -phred33 ../012m_L1_1.fq ../012m_L1_2.fq 012m_Trimmomatic.1.fq 012m_Trimmomatic.unpaired.1.fq 012m_Trimmomatic.2.fq 012m_Trimmomatic.unpaired.2.fq \
ILLUMINACLIP:/opt/biosoft/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75 TOPHRED33
cd ../
# correct Illumina reads | Assessment of genome size and heterozygosity
mkdir Finderrors
source ~/.bash.pacbio
ErrorCorrectReads.pl PHRED_ENCODING=33 READS_OUT=012m FILL_FRAGMENTS=0 KEEP_KMER_SPECTRA=1 PAIRED_READS_A_IN=../fastuniq/012m.fastuniq.1.fastq PAIRED_READS_B_IN=../fastuniq/012m.fastuniq.1.fastq PLOIDY=2 PAIRED_SEP=251 PAIRED_STDEV=48
#Wed Sep 05 22:01:16 2018 (KSC): Genome size estimate        =     82,542,732 bases
#Wed Sep 05 22:01:16 2018 (KSC): Genome size estimate CN = 1 =     65,412,360 bases (  79.2 % )
#Wed Sep 05 22:01:16 2018 (KSC): Genome size estimate CN > 1 =     17,130,372 bases (  20.8 % )
#Wed Sep 05 22:01:16 2018 (KSC): Coverage estimate           =             97 x
#Wed Sep 05 22:01:16 2018 (KSC): Bias stddev at scale > K    =           0.20
#Wed Sep 05 22:01:16 2018 (KSC): Base error rate estimate    =         0.0006 (Q = 32.3)
#Wed Sep 05 22:01:16 2018 (KSC): SNP rate: always verify with kmer spectrum plot.
#Wed Sep 05 22:01:16 2018 (KSC): Ploidy                      =              2
#Wed Sep 05 22:01:16 2018 (KSC): SNP rate                   ~=           1/48
#Wed Sep 05 22:01:16 2018 (KSC): SNPs closer than K         ~=             65 %
#Wed Sep 05 22:01:16 2018 (KSC): --------------------------------------------------------------
#PERFSTAT: estimated genome size in bases [genome_size_est] = 82542732
#PERFSTAT: % genome estimated to be repetitive (at K=25 scale) [genome_repetitiveness_est] = 20.0
#PERFSTAT: estimated genome coverage by fragment reads [genome_cov_est] = 97
#PERFSTAT: estimated standard deviation of sequencing bias (at K=25 scale) [bias_stddev_est] = 0.20

# kmer plot
cd 012m.fastq.kspec
perl -p -i -e 's/\@fns = \(\"frag_reads_filt.25mer.kspec\", \"frag_reads_edit.24mer.kspec\", \"frag_reads_corr.25mer.kspec\"\);\n//' /opt/biosoft/ALLPATHS-LG/bin/KmerSpectrumPlot.pl
KmerSpectrumPlot.pl SPECTRA=1 FREQ_MAX=255
perl -p -i -e 's/\@fns = \(\"frag_reads_filt.25mer.kspec\", \"frag_reads_edit.24mer.kspec\", \"frag_reads_corr.25mer.kspec\"\);\n//' /opt/biosoft/ALLPATHS-LG/bin/KmerSpectrumPlot.pl
convert kmer_spectrum.cumulative_frac.log.lin.eps kmer_spectrum.cumulative_frac.log.lin.png
convert kmer_spectrum.distinct.lin.lin.eps kmer_spectrum.distinct.lin.lin.png
convert kmer_spectrum.distinct.log.log.eps kmer_spectrum.distinct.log.log.png
cd ../
# Using LoRDEC to modify PacBio Reads
mkdir LoRDEC
cd LoRDEC
lordec-correct -2 ../Finderrors/012m.paired.A.fastq ../Finderrors/012m.paired.B.fastq -i ../012m.subreads.fq -k 19 -o pacbio.LoRDEC.corrected.fasta -s 3 -T 32 &> lordec-correct.log
seqkit seq -u pacbio.LoRDEC.corrected.fasta > pacbio.corrected.fasta
cd ../
########################################################### Genome assembly | polish #################################################
# Genome assembly using canu
mkdir canu
cd canu
source ~/.bashrc.pacbio
canu genomeSize=83000000 useGrid=false corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" rawErrorRate=0.150 correctedErrorRate=0.035 -p 012m -d ./ -pacbio-raw ../LoRDEC/pacbio.corrected.fasta minReadLength=1000 minOverlapLength=500
mkdir genome_feature_analysis
cd genome_feature_analysis
# Genome repeat analysis using RepeatMasker
mkdir repeat_analysis
cd repeat_analysis
fasta_no_blank.pl ../../012m.contigs.fasta > genome.fasta
mkdir repeatMasker
/opt/biosoft/RepeatMasker/RepeatMasker -pa 32 -e ncbi -species "Basidiomycota" -dir repeatMasker/ -gff genome.fasta
# Using repeatModeler to search for repetitive sequences of this species
mkdir repeatModeler
cd repeatModeler
/opt/biosoft/RepeatModeler-open-1.0.11/BuildDatabase -name Armillaria_012m -engine ncbi ../genome.fasta
/opt/biosoft/RepeatModeler-open-1.0.11/RepeatModeler -engine ncbi -pa 32 -database Armillaria_012m
/opt/biosoft/RepeatMasker/RepeatMasker -pa 32 -e ncbi -lib Armillaria_012m-families.fa -dir ./ -gff ../genome.fasta
cd ../
# merge the results of RepeatMasker and repeatModeler
merge_repeatMasker_out.pl repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats
# masked genome
maskedByGff.pl genome.repeat.gff3 genome.fasta --mask_type softmask > 012m_genome.softmask.fasta
# Using HaploMerger2 to obtain haploid genome
mkdir -p /opt/biosoft/HaploMerger2_20180603/mydata/ 
mv 012m_genome.softmask.fasta /opt/biosoft/HaploMerger2_20180603/mydata/
cd /opt/biosoft/HaploMerger2_20180603/mydata/
cp /opt/biosoft/HaploMerger2_20180603/project_template/all_lastz.ctl ./
cp /opt/biosoft/HaploMerger2_20180603/project_template/hm.* ./
cp /opt/biosoft/HaploMerger2_20180603/project_template/scoreMatrix.q ./
cp /opt/biosoft/HaploMerger2_20180603/project_template/run_all.batch ./
mkdir libraries
cd libraries
cp ~/data/012m/Finderrors/012m.paired.?.fastq ./
echo "# need no insert std

#maximal read length
max_rd_len=150

[LIB]
#average insert size
avg_ins=251
#if sequence needs to be reversed 
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#in which order the reads are used while scaffolding
rank=1
#fastq file for read 1 
q1=/opt/biosoft/HaploMerger2_20180603/mydata/libraries/012m.paired.A.fastq
q2=/opt/biosoft/HaploMerger2_20180603/mydata/libraries/012m.paired.B.fastq " >  gapCloser.cfg

echo 'Lib1 bwa/opt/biosoft/HaploMerger2_20180603/mydata/libraries/012m.paired.A.fastq /opt/biosoft/HaploMerger2_20180603/mydata/libraries/012m.paired.B.fastq 251 0.25 FR' > sspace_libraries.list
cd ../
mv 012m_genome.softmask.fasta genome.fa
gzip genome.fa
sh ./run_all.batch >run_all.log 2>&1
cd ~/data/012m
# Genome upgrade using FinisherSC
mkdir FinisherSC
cd FinisherSC
ln -s ../LoRDEC/pacbio.corrected.fasta raw_reads.fasta
cp /opt/biosoft/HaploMerger2_20180603/mydata/genome_A_ref_C_D_E.fa.gz ./
gunzip genome_A_ref_C_D_E.fa.gz
mv genome_A_ref_C_D_E.fa contigs.fasta
python /opt/biosoft/finishingTool/finisherSC.py -par 32 -o contigs.fasta_improved3.fasta ./ /opt/biosoft/finishingTool/MUMmer3.23/
cd ../
# plishing with pilon
mkdir pilon
cd pilon
seqkit seq -u ../FinisherSC/improved3.fasta > 012m_FinisherSC.fasta
bowtie2-build 012m_FinisherSC.fasta 012m
bowtie2 -p 32 -x 012m -1 ../Finderrors/012m.paired.A.fastq -2 ../Finderrors/012m.paired.B.fastq -S 012m.sam 2> 012m.bowtie2.log
java -jar /opt/biosoft/picard-tools-2.17.3/picard.jar SortSam I=012m.sam O=012m.bam SO=coordinate
samtools index 012m.bam -@ 32
java -Xmx60G -jar /opt/biosoft/pilon/pilon-1.22.jar --genome 012m_FinisherSC.fasta --bam 012m.bam --output 012m_pilon --outdir pilon_out --fix all --threads 32
cd ../
################################################################## genome_feature_analysis#############################################################
# Assembly status statistics
mkdir genome_statistic
cd genome_statistic
genome_seq_clear.pl ../pilon/pilon_out/012m_pilon.fasta --seq_prefix A_012m_scafford_ > Armillaria_012m_genome.fasta
cal_seq_length.pl Armillaria_012m_genome.fasta > seq_length.txt
cd ../