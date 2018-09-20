################################################### gene prediction ####################################################################
mkdir RNA
#################Storing RNA-seq data##################
cd RNA
mkdir Trimmomatic
cd Trimmomatic
for i in `ls ~/data/012m/RNA/*.1.fq`
do
    i=${i/*\//}
    i=${i/.1.fq/}
java -jar /opt/biosoft/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 ~/data/012m/RNA/$i.1.fq ~/data/012m/RNA/$i.2.fq $i.1.fastq $i.1.unpaired.fastq $i.2.fastq $i.2.unpaired.fastq \
ILLUMINACLIP:/opt/biosoft/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 TOPHRED33
done
cd ../
# RNA-seq data quality control
mkdir hisat2_gene_p
ln -s ../Trimmomatic/??.?.fastq ./
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/012m_genome.softmask.fasta ./
hisat2-build -p 32 012m_genome.softmask.fasta genome
for i in `ls *.1.fastq`
do
    i=${i/.1.fastq/}
echo "hisat2 -x genome -p 32 --min-intronlen 20 --max-intronlen 5000 --rna-strandness RF -1 $i.1.fastq -2 $i.2.fastq -S $i.sam --new-summary --summary-file $i.hisat2.summary"
done > command.hisat2.list
sh command.hisat2.list
cd ../../
# Gene prediction using GeneMark-ES/ET
mkdir genemark_es_et
cd genemark_es_et
# obtain hints file
mkdir making_hints
cd making_hints
samtools merge -@ 32 rnaseq.bam ~/data/012m/RNA/hisat2_gene_p/*.sam
samtools sort -@ 32 -O bam -o rnaseq.sort.bam rnaseq.bam
bam2hints --intronsonly --in=rnaseq.sort.bam --out=hints.gff
cd ../
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/012m_genome.softmask.fasta genome.fasta
hints2genemarkETintron.pl genome.fasta making_hints/hints.gff > genemartET.intron.gff
gmes_petap.pl --sequence genome.fasta --ET genemartET.intron.gff --max_contig 10000000 --min_contig 1000 --fungus --et_score 4 --cores 32
/opt/biosoft/pasa-v2.2.0/misc_utilities/gtf_to_gff3_format.pl genemark.gtf genome.fasta > genemark.gff3
cd ../
rm -rf hisat2_gene_p
mkdir hisat2_gene_p1
cd hisat2_gene_p1
ln -s ../Trimmomatic/??.?.fastq ./
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/012m_genome.softmask.fasta ./
cp ~/data/012m/RNA/genemark_es_et/genemark.gtf ./
hisat2_extract_splice_sites.py genemark.gtf > splitSites.txt
hisat2_extract_exons.py genemark.gtf > exonSites.txt
hisat2-build -p 32 --ss splitSites.txt --exon exonSites.txt 012m_genome.softmask.fasta genome
for i in `ls *.1.fastq`
do
    i=${i/.1.fastq/}
echo "hisat2 -x genome -p 32 --min-intronlen 20 --max-intronlen 5000 --rna-strandness RF -1 $i.1.fastq -2 $i.2.fastq -S $i.sam --new-summary --summary-file $i.hisat2.summary"
done > command.hisat2.list
sh command.hisat2.list
cd ../
# Using Trinity to de novo assemble (RNA-seq)
mkdir Trinity
cd Trinity
ln -s ../Trimmomatic/??.?.fastq ./
Trinity --seqType fq --max_memory 30G --left C1.1.fastq,C2.1.fastq,C3.1.fastq,G1.1.fastq,G2.1.fastq,G3.1.fastq --right C1.2.fastq,C2.2.fastq,C3.2.fastq,G1.2.fastq,G2.2.fastq,G3.2.fastq \
--SS_lib_type RF --CPU 10 --jaccard_clip --normalize_reads --normalize_max_read_cov 100 --output trinity_denovo --bflyCalculateCPU &> trinity_denovo.log

# de novo assemble (summary)
/opt/biosoft/Trinity-v2.7.0/util/TrinityStats.pl trinity_denovo/Trinity.fasta > trinity_denovo/Trinity.fasta.stats


# Assembly with guided genome 
cd Trinity
samtools merge -@ 32 merged.bam ~/data/012m/RNA/hisat2_gene_p1/*.sam
samtools sort -@ 32 -O BAM -o merged.sort.bam merged.bam
Trinity --max_memory 30G --SS_lib_type RF --CPU 20 --jaccard_clip --normalize_reads --normalize_max_read_cov 100 --genome_guided_bam merged.sort.bam \
--genome_guided_max_intron 10000 --output trinity_genomeGuided --bflyCalculateCPU &> trinity_genomeGuided.log
cd ../../
## obtain best gene model (for HMM trainning )
mkdir pasa
cd pasa
cp ../genome_feature_analysis/repeat_analysis/012m_genome.softmask.fasta ./genome.fasta
cat ~/data/012m/RNA/Trinity/trinity_denovo/Trinity.fasta >> transcripts.fasta
cat ~/data/012m/RNA/Trinity/trinity_genomeGuided/Trinity-GG.fasta >> transcripts.fasta
perl -e 'while (<>) { print "$1\n" if />(\S+)/ }' ~/data/012m/RNA/Trinity/trinity_denovo/Trinity.fasta > tdn.accs
# end-trimming for transcripts sequences (vector, adaptor, primer, polyA/T tails)
seqclean transcripts.fasta -v /opt/biosoft/pasa-v2.2.0/seqclean/UniVec
# Generate alignment files
cp /opt/biosoft/pasa-v2.2.0/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
echo "perl -p -i -e 's/MYSQLDB=.*/MYSQLDB=pasa_201809/' alignAssembly.config" | sh
# Generate MySQL database and tables
/opt/biosoft/pasa-v2.2.0/scripts/create_mysql_cdnaassembly_db.dbi -r -c alignAssembly.config -S /opt/biosoft/pasa-v2.2.0/schema/cdna_alignment_mysqlschema
# obtain nr transcript sequence (transcript sequence => genome) Alternative splicing information
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -R -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta --ALIGNERS gmap,blat --CPU 32 \
--stringent_alignment_overlap 30.0 --transcribed_is_aligned_orient --TDN tdn.accs --MAX_INTRON_LENGTH 20000  --TRANSDECODER &> pasa.log
# Chain specific sequencing requires adding parameters --transcribed_is_aligned_orient
# Fungi and other small genomes, because of relatively dense genes, need to add parameters --stringent_alignment_overlap
# Constructing comprehensive transcriptome database
#perl -p -i -e 's#TR\\d\+\\\|#TRINITY_DN\\d\+_#' /opt/biosoft/pasa-v2.2.0/scripts/build_comprehensive_transcriptome.dbi
/opt/biosoft/pasa-v2.2.0/scripts/build_comprehensive_transcriptome.dbi -c alignAssembly.config -t transcripts.fasta.clean
# ORF prediction and gene prediction
/opt/biosoft/pasa-v2.2.0/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta pasa_201809.assemblies.fasta --pasa_transcripts_gff3 pasa_201809.pasa_assemblies.gff3 -S
#Extract complete genes (exon >= 3、cds lenth >=900bp、cds region accounted for exon region ratio>=0.6、sort intron by lenth (retain %95 )、sort cds by lenth (retain %95 )
pasa_extract_best_candidate_geneModels.pl pasa_201809.assemblies.fasta.transdecoder.genome.gff3 pasa_201809.assemblies.fasta.transdecoder.cds 3 900 0.60 0.95 0.95 > best_candidates.gff3
# all vs all identity < 70%
/opt/biosoft/EVidenceModeler-1.1.1/EvmUtils/gff3_file_to_proteins.pl best_candidates.gff3 genome.fasta prot > best_candidates.fasta
remove_redundant_high_identity_genes.pl best_candidates.gff3 best_candidates.fasta 32 0.70 > best_candidates.lowIdentity.gff3 2> remove_redundant_high_identity_genes.log
# best_candidates.lowIdentity.gff3 can be used for other gene prediction software (HMM trainning) 。

# AUGUSTUS Training Gene Model
mkdir -p augustus/training
cd ~/data/012m/augustus/training
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/012m_genome.softmask.fasta ./genome.fasta
gff2gbSmallDNA.pl ../../pasa/best_candidates.lowIdentity.gff3 genome.fasta 100 genes.raw.gb
# Remove the wrong genes
new_species.pl --species=for_bad_genes_removing
etraining --species=for_bad_genes_removing --stopCodonExcludedFromCDS=false genes.raw.gb 2> train.err
cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
filterGenes.pl badgenes.lst genes.raw.gb > genes.gb
# The genes used for training were divided into 2 parts。
randomSplit.pl genes.gb 200
new_species.pl --species=A_012m
# First training
etraining --species=A_012m genes.gb.train > train.out
perl -e 'open IN, "train.out"; while (<IN>) { $tag = $1 if m/tag:.*\((.*)\)/; $taa = $1 if m/taa:.*\((.*)\)/; $tga = $1 if m/tga:.*\((.*)\)/; } while (<>) { s#/Constant/amberprob.*#/Constant/amberprob                   $tag#; s#/Constant/ochreprob.*#/Constant/ochreprob                   $taa#; s#/Constant/opalprob.*#/Constant/opalprob                    $tga#; print }' /opt/biosoft/augustus/config/species/A_012m/A_012m_parameters.cfg > 11
mv 11 /opt/biosoft/augustus/config/species/A_012m/A_012m_parameters.cfg
augustus --species=A_012m genes.gb.test | tee firsttest.out
# Use optimize_augustus.pl to circulate training to find optimal parameters
# Genes.gb.train contains 2033 gene models. Segmenting it again, 1800 of them were used to detect the accuracy of HMM model parameters optimization, and 233 of them were only added to the training process。
ln -s genes.gb.train.test training.gb.onlytrain
optimize_augustus.pl --species=A_012m --rounds=5 --cpus=64 --kfold=64 --onlytrain=training.gb.onlytrain genes.gb.train.train > optimize.out 2> /dev/null
# second training
etraining --species=A_012m genes.gb.train
augustus --species=A_012m genes.gb.test | tee secondtest.out
cd /opt/biosoft/augustus/config/species/A_012m
cp A_012m_exon_probs.pbl A_012m_exon_probs.pbl.withoutCRF
cp A_012m_igenic_probs.pbl A_012m_igenic_probs.pbl.withoutCRF
cp A_012m_intron_probs.pbl A_012m_intron_probs.pbl.withoutCRF
cd -
etraining --species=A_012m --CRF=1 genes.gb.train
augustus --species=A_012m genes.gb.test | tee secondtest.out.withCRF
# Compare the accuracy of two cases of CRF and non CRF. In general, the accuracy of CRF training is higher. If the accuracy of CRF training is low, the backup parameter file can be restored back。
cd /opt/biosoft/augustus/config/species/A_012m
cp A_012m_exon_probs.pbl.withoutCRF A_012m_exon_probs.pbl 
cp A_012m_igenic_probs.pbl.withoutCRF A_012m_igenic_probs.pbl 
cp A_012m_intron_probs.pbl.withoutCRF A_012m_intron_probs.pbl
cd -
cd ../../

## get HMM model with snap
mkdir snap
cd snap
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/012m_genome.softmask.fasta ./genome.fasta
fathom genome.ann genome.dna -gene-stats &> gene-stats.log
fathom genome.ann genome.dna -validate &> validate.log
perl -ne 'print "$1\n" if /.*:\s+(\S+)\s+OK/' validate.log > zff2keep.txt
perl -e 'open IN, "zff2keep.txt"; while (<IN>) { chomp; $keep{$_} = 1; } while (<>) { if (m/>/) { print; } else { chomp; @_ = split /\t/; print "$_\n" if exists $keep{$_[-1]}; } }' genome.ann > out;
mv out genome.ann
fathom genome.ann genome.dna -categorize 200
rm alt.* err.* olp.* wrn.*
fathom genome.ann genome.dna -export 200 -plus
mkdir params; cd params
forge ../export.ann ../export.dna
cd ..
hmm-assembler.pl species params > species.hmm

mkdir maker
cd maker
######################################Download "Basidiomycota" uniprotKB proteins and (/opt/biosoft/BUSCO_V3.0.2b/database/basidiomycota_odb9/ singel copy proteins) as homolog proteins #############################################################
perl -pi -e 's/\|/_/g' homolog.fasta
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/genome.fasta ./genome.fasta
ln -s ~/data/012m/RNA/Trinity/trinity_denovo/Trinity.fasta ./
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/repeatModeler/RM_120853.TueSep111055302018/consensi.fa.classified consensi.fa.classified
ln -s ~/data/012m/genemark_es_et/output/gmhmm.mod ./
ln -s ~/data/012m/snap/species.hmm ./
perl -p -i -e 's/^genome=.*/genome=genome.fasta/; s/^est=.*/est=Trinity.fasta/; s/^protein=.*/protein=homolog.fasta/; s/^model_org=.*/model_org=Basidiomycota/; \
s/^rmlib=.*/rmlib=consensi.fa.classified/; s/^augustus_species=.*/augustus_species=A_012m/; s/^snaphmm=.*/snaphmm=species.hmm/; s/^gmhmm=.*/gmhmm=gmhmm.mod/; \
s/^est2genome=.*/est2genome=1/; s/^protein2genome=.*/protein2genome=1/; s/^trna=.*/trna=1/; s/^correct_est_fusion=.*/correct_est_fusion=1/; s/^keep_preds=.*/keep_preds=1/;' maker_opts.ctl
# 准备配置文件
maker -CTL
sudo hostname localhost
echo 'MPD_SECRETWORD=mr45-j9z' > ~/.mpd.conf
chmod 600 ~/.mpd.conf
mpd &
/usr/bin/mpiexec -n 64 maker
cd genome.maker.output
gff3_merge -d genome_master_datastore_index.log
grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3
############ find an error ######### remove #################
#perl -pi -e 's/Parent=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1,maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-2/Parent=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1/' genome.maker.gff3
#grep -v 'maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-2' genome.maker.gff3 > 11
#mv 11 genome.maker.gff3
#A_012m_scafford_07	maker	gene	1243672	1244989	.	-	.	ID=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2;Name=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2
#A_012m_scafford_07	maker	mRNA	1243672	1244989	1318	-	.	ID=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1;Parent=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2;Name=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1;_AED=0.01;_eAED=0.01;_QI=0|1|1|1|0.4|0.16|6|202|371
#A_012m_scafford_07	maker	mRNA	1243672	1244989	1318	-	.	ID=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-2;Parent=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2;Name=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-2;_AED=0.00;_eAED=0.00;_QI=0|1|1|1|0.5|0.4|5|202|371
#A_012m_scafford_07	maker	exon	1243672	1244989	.	-	.	ID=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1:exon:3421;Parent=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1,maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-2
#A_012m_scafford_07	maker	three_prime_UTR	1243672	1243873	.	-	.	ID=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1:three_prime_utr;Parent=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-1
#A_012m_scafford_07	maker	three_prime_UTR	1243672	1243873	.	-	.	ID=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-2:three_prime_utr;Parent=maker-A_012m_scafford_07-exonerate_est2genome-gene-12.2-mRNA-2
gff3_clear.pl --prefix A_012mGene genome.maker.gff3 > ../genome.maker.gff3
cd ../../
## Update the result of maker with Pasa
mkdir pasa_maker
cd pasa/
cp /opt/biosoft/pasa-v2.2.0/pasa_conf/pasa.annotationCompare.Template.txt annotationCompare.config
PasaMysqlDB=`perl -ne 'print $1 if m/MYSQLDB=(\S+)/' alignAssembly.config`
echo "perl -p -i -e 's/MYSQLDB=.*/MYSQLDB=$PasaMysqlDB/' annotationCompare.config" | sh
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c annotationCompare.config -A -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta -L --annots_gff3 ../maker/genome.maker.gff3
first_update_gff3=`ls *gene_structures_post_PASA_updates*.gff3 -t | head -n 1`
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c annotationCompare.config -A -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta -L --annots_gff3 $first_update_gff3
second_update_gff3=`ls *gene_structures_post_PASA_updates*.gff3 -t | head -n 1`
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c annotationCompare.config -A -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta -L --annots_gff3 $second_update_gff3
third_update_gff3=`ls *gene_structures_post_PASA_updates*.gff3 -t | head -n 1`
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c annotationCompare.config -A -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta -L --annots_gff3 $third_update_gff3
fourth_update_gff3=`ls *gene_structures_post_PASA_updates*.gff3 -t | head -n 1`
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c annotationCompare.config -A -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta -L --annots_gff3 $fourth_update_gff3
five_update_gff3=`ls *gene_structures_post_PASA_updates*.gff3 -t | head -n 1`
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c annotationCompare.config -A -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta -L --annots_gff3 $five_update_gff3
six_update_gff3=`ls *gene_structures_post_PASA_updates*.gff3 -t | head -n 1`
/opt/biosoft/pasa-v2.2.0/scripts/Launch_PASA_pipeline.pl -c annotationCompare.config -A -g genome.fasta -t transcripts.fasta.clean -T -u transcripts.fasta -L --annots_gff3 $six_update_gff3
seven_update_gff3=`ls *gene_structures_post_PASA_updates*.gff3 -t | head -n 1`
cp $seven_update_gff3 ../pasa_maker/pasa_maker.gff3
cd ../pasa_maker/
ln -s ~/data/012m/genome_feature_analysis/repeat_analysis/genome.fasta ./genome.fasta
gff3_clear.pl --prefix A_012mGene pasa_maker.gff3 > genome.maker.gff3
gff3ToGtf.pl genome.fasta genome.maker.gff3 > genome.gtf
eukaryotic_gene_model_statistics.pl genome.gtf genome.fasta out > out.gene_model_statistic.txt
bestGeneModels.pl genome.maker.gff3 > genome.bestGeneModels.gff3 2> geneModelsStatistic
gff3ToGtf.pl genome.fasta genome.bestGeneModels.gff3 > genome.bestGeneModels.gtf
eukaryotic_gene_model_statistics.pl genome.bestGeneModels.gtf genome.fasta bestGeneModels > bestGeneModels.gene_model_statistic.txt
cd ../
