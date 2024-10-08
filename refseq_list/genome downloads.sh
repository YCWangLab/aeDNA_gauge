########################download the 350 bacteria, 100 Archaea, 50 fungi, 50 plants, 15 vertebrate mammalian, 15 vertebrate other and 20 invertebrate full genome assemblies from NCBI RefSeq

#download the assembly summary file from NCBI genome refseq manually
cd work_folder/genome/refseq_list

#make random number matix
shuf -i 1-249321 -n 350 > b_list.txt
shuf -i 1-1250 -n 100 > a_list.txt
shuf -i 1-444 -n 50 > f_list.txt
shuf -i 1-150 -n 50 > p_list.txt
shuf -i 1-190 -n 15 > vm_list.txt
shuf -i 1-280 -n 15 > vo_list.txt
shuf -i 1-296 -n 20 > i_list.txt

#subset info from the refseq summary list 
for file in *summary.txt
do
  bname=`echo $file | cut -d"_" -f1`ZZ 
  awk 'NR == FNR{a[$0]; next};FNR in a' ${bname}_list.txt $file >> genome_list.txt 
done
cut -f20 genome_list.txt > ftp_list.txt

#download genomes
cat refseq_list/ftp_list.txt | while read line 
do
  ftp=`echo $line | sed "s/https/ftp/"`
  acc=`echo $ftp | cut -d\\/ -f10`
  echo ${ftp}/${acc}_genomic.fna.gz
  wget ${ftp}/${acc}_genomic.fna.gz
done

gzip -d *.gz
