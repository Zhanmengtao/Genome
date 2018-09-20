#################################################### genome_feature_analysis #######################################################
# Assembly status statistics
mkdir genome_statistic
cd genome_statistic
genome_seq_clear.pl ../pilon/pilon_out/012m_pilon.fasta --seq_prefix A_012m_scafford_ > Armillaria_012m_genome.fasta
cal_seq_length.pl Armillaria_012m_genome.fasta > seq_length.txt
cd ../
# genome_feature_analysis
mkdir genome_feature_analysis
cd genome_feature_analysis
#  RepeatMasker 
mkdir repeat_analysis
cd repeat_analysis
fasta_no_blank.pl ~/data/012m/genome_statistic/Armillaria_012m_genome.fasta > genome.fasta
mkdir repeatMasker
/opt/biosoft/RepeatMasker/RepeatMasker -pa 32 -e ncbi -species "Basidiomycota" -dir repeatMasker/ -gff genome.fasta
# repeatModeler 
mkdir repeatModeler
cd repeatModeler
/opt/biosoft/RepeatModeler-open-1.0.11/BuildDatabase -name Armillaria_012m -engine ncbi ../genome.fasta
/opt/biosoft/RepeatModeler-open-1.0.11/RepeatModeler -engine ncbi -pa 32 -database Armillaria_012m
/opt/biosoft/RepeatMasker/RepeatMasker -pa 32 -e ncbi -lib Armillaria_012m-families.fa -dir ./ -gff ../genome.fasta
cd ../
# merge
merge_repeatMasker_out.pl repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats
# masked genome
maskedByGff.pl genome.repeat.gff3 genome.fasta --mask_type softmask > 012m_genome.softmask.fasta
cd ../
# RNAmmer
mkdir RNAmmer
cd RNAmmer
ln -s ~/data/012m/genome_statistic/Armillaria_012m_genome.fasta ./
/opt/biosoft/rnammer-1.2/rnammer -S euk -multi -f rRNA.fasta -h rRNA.hmmreport -xml rRNA.ml -gff rRNA.gff2 Armillaria_012m_genome.fasta
rRNAmmer_gff2gff3.pl rRNA.gff2 > rRNA.gff3
cd ../
# tRNA
mkdir tRNAscan-SE
cd tRNAscan-SE
ln -s ~/data/012m/genome_statistic/Armillaria_012m_genome.fasta ./
tRNAscan-SE -o 012m_tRNA.out -f 012m_tRNA.ss -m 012m_tRNA.stats Armillaria_012m_genome.fasta
tRNAscanSE2GFF3.pl 012m_tRNA.out 012m_tRNA.ss > 012m_tRNA.gff3
cd ../../
