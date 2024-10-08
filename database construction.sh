
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

#subset mitochondrion genome
grep "mitochondrion" *.fna >> mito.txt
cat mito.txt | awk -F '[> ]' '{print $2}' >> mito_id.txt

for file in work_folder/genome/animal/vm/mit/*.fna
do
  bname=`basename $file | cut -d. -f1`
  cat $file | seqkit grep -f mito_id.txt -o $bname.fna
done

#subset chlo genome

for file in work_folder/genome/plant/mito/chlo/*.fna
do
  bname=`basename $file | cut -d. -f1`
  cat $file | seqkit grep -f work_folder/genome/plant/mito/chlo/chlo_id.txt -o $bname.fna
done

for file in work_folder/genome/plant/mito/chlo/*.fna
do
  bname=`basename $file | cut -d. -f1`
  cat $file | seqkit grep -f work_folder/genome/plant/mito/chlo/mito_id.txt -o $bname.fna
done

for file in work_folder/genome/plant/mit/chlo/*.fna
do
  bname=`basename $file | cut -d. -f1`
  cat $file | seqkit grep -v -f ../chlo_id.txt | seqkit grep -v -f ../mito_id.txt >> $bname.fna
done



######################### sampling
for file in *.fa
do
  samtools faidx $file
done

####modern sample
#add micro reads
a=1
while(($a<=30))
do
    c=1
    for file in work_folder/genome/originalgenome/*.fna
    do
        b=`sed -n ${c}p work_folder/result/sample/modern/$a.txt`
        gargammel/src/fragSim --seed $a -n $b -f work_folder/length/300_len_dis.txt -o $bname.fa.gz $file
        gzip -d $bname.fa.gz
        cat $bname.fa >> $a.fa
        rm $bname.fa
        let "c++"
    done
    let "a++"
done

#add animal reads
a=1
while(($a<=30))
do
  for file in work_folder/genome/animal/wholegenome/*.fna
  do
    bname=`basename $file | cut -d. -f1`
    b=`shuf -i 10-10000 -n 1`
	  echo $b >> $a.animal.txt
    gargammel/src/fragSim --seed $a -n $b -f work_folder/length/300_len_dis.txt -o $bname.fa.gz $file
	  gzip -d $bname.fa.gz
	  cat $bname.fa >> $a.animal.fa
	  rm $bname.fa
	  c=`expr $b / 100`
	  echo $c >> $a.animal.mit.txt
	  gargammel/src/fragSim --seed $a -n $c -f work_folder/length/300_len_dis.txt -o $bname.fa.gz work_folder/genome/animal/mit/$bname.fna
	  gzip -d $bname.fa.gz
	  cat $bname.fa >> $a.animal.mi.fa
	  rm $bname.fa
  done
  let "a++"
done

#add plant reads
a=1
while(($a<=30))
do
  for file in work_folder//genome/plant/wholegenome/*.fna
  do
    bname=`basename $file | cut -d. -f1`
    b=`shuf -i 10-10000 -n 1`
    echo $b >> $a.plant.txt
    gargammel/src/fragSim --seed $a -n $b -f work_folder/length/300_len_dis.txt -o $bname.fa.gz $file
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.plant.fa
    rm $bname.fa
    c=`expr $b / 50`
    echo $c >> $a.plant.mit.txt
    gargammel/src/fragSim --seed $a -n $c -f work_folder/length/300_len_dis.txt -o $bname.fa.gz work_folder/genome/plant/mit/$bname.fna
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.plant.mit.fa
    rm $bname.fa
    gargammel/src/fragSim --seed $a -n $b -f work_folder/length/300_len_dis.txt -o $bname.fa.gz work_folder/genome/plant/chlo/$bname.fna
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.plant.chlo.fa
    rm $bname.fa
  done
  let "a++"
done

a=1
while(($a<=30))
do
    cat work_folder/result/sample/modern/$a.fa work_folder/result/sample/modern/macro/$a.*.fa >> $a.modern.fa
    let "a++"
done

####young sample
#fragmentation
a=1
while(($a<=30))
do
    c=1
    for file in work_folder/genome/originalgenome/latest/*.fna	
    do
        b=`sed -n ${c}p work_folder/result/sample/modern/$a.txt`
        echo $b >> $a.txt
        gargammel/src/fragSim --seed $a -n $b -f work_folder/length/new/young_len_dis.txt -o $bname.fa.gz $file
	gzip -d $bname.fa.gz
	cat $bname.fa >> $a.fa
	rm $bname.fa
	let "c++"
    done
  let "a++"
done

a=1
while(($a<=30))
do
  c=1
  for file in work_folder//genome/animal/wholegenome/*.fna
  do
    bname=`basename $file | cut -d. -f1`
    b=`sed -n ${c}p work_folder/result/sample/modern/macro/$a.animal.txt`
    echo $b >> $a.animal.txt
    gargammel/src/fragSim --seed $a -n $b -f work_folder/length/new/young_len_dis.txt -o $bname.fa.gz $file
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.animal.fa
    rm $bname.fa
    e=`expr $b / 100`
    echo $e >> $a.animal.mit.txt
    gargammel/src/fragSim --seed $a -n $e -f work_folder/length/new/young_len_dis.txt -o $bname.fa.gz work_folder//genome/animal/mit/$bname.fna
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.animal.mit.fa
    rm $bname.fa
    let "c++"
  done
  let "a++"
done

a=1
while(($a<=30))
do
  c=1
  for file in work_folder/genome/plant/wholegenome/*.fna
  do
    bname=`basename $file | cut -d. -f1`
    b=`sed -n ${c}p work_folder/result/sample/modern/macro/$a.plant.txt`
    echo $b >> $a.plant.txt
    gargammel/src/fragSim --seed $a -n $b -f work_folder/length/young_len_dis.txt -o $bname.fa.gz $file
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.plant.fa
    rm $bname.fa
    e=`expr $b / 50`
    echo $e >> $a.plant.mit.txt
    gargammel/src/fragSim --seed $a -n $e -f work_folder/length/new/young_len_dis.txt -o $bname.fa.gz work_folder/genome/plant/mit/$bname.fna
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.plant.mit.fa
    rm $bname.fa
    gargammel/src/fragSim --seed $a -n $b -f work_folder/length/new/young_len_dis.txt -o $bname.fa.gz work_folder/genome/plant/chlo/$bname.fna
    gzip -d $bname.fa.gz
    cat $bname.fa >> $a.plant.chlo.fa
    rm $bname.fa	
    let "c++"
  done
  let "a++"
done

#simulate random mutations with Mutation-Simulator (0.0001 for Micro, 0.00001 for Macro and 0.0002 for Organelle reads)
cd work_folder/result/sample/young/mutation
for file in work_folder/result/sample/young/*.fa
do
  bname=`basename $file | cut -d. -f1`
  grep '>' $file >> $bname.header.txt
  awk 'NR%2==0' $file | tr '\n' '\t' | sed 's/\t/N/g' | sed '1i\>1' >> $bname.seq.fa
  mutation-simulator $bname.seq.fa args -sn 0.0001 -titv 1
  cat $bname.fa | sed '1d' | sed 's/N/\n/g' >> $bname.mutation.fa
  paste -d* $bname.header.txt $bname.mutation.fa | sed 's/*/\n/g' >> $bname.simulation.fa
  diff ../$bname.fa $bname.simulation.fa -y | grep -B 1 '|' | sed  '/--/d' | grep '>' | awk '{print $2}' | awk -F '[> ]' '{print $2}' >> acession.txt
  cat $bname.simulation.fa | seqkit grep -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed 's/$/:v/;n' | sed '$ s/$/ \n/' >> $bname.variation.fa
  cat $bname.simulation.fa | seqkit grep -v -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed '$ s/$/ \n/'>> $bname.invariation.fa
  cat $bname.variation.fa $bname.invariation.fa >> $bname.sample.fa
  rm acession.txt $bname.variation.fa $bname.invariation.fa $bname.header.txt $bname.seq.fa $bname.mutation.fa $bname.simulation.fa $bname.fa $bname.vcf $bname.seq.fa.fai
done

cd work_folder/result/sample/young/macro/mutation
for file in work_folder/result/sample/young/macro/*.animal.fa
do
  bname=`basename $file | cut -d. -f1-2`
  grep '>' $file >> $bname.header.txt
  awk 'NR%2==0' $file | tr '\n' '\t' | sed 's/\t/N/g' | sed '1i\>1' >> $bname.seq.fa
  mutation-simulator $bname.seq.fa args -sn 0.00001 -titv 1
  cat $bname.fa | sed '1d' | sed 's/N/\n/g' >> $bname.mutation.fa
  paste -d* $bname.header.txt $bname.mutation.fa | sed 's/*/\n/g' >> $bname.simulation.fa
  diff ../$bname.fa $bname.simulation.fa -y | grep -B 1 '|' | sed  '/--/d' | grep '>' | awk '{print $2}' | awk -F '[> ]' '{print $2}' >> acession.txt
  cat $bname.simulation.fa | seqkit grep -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed 's/$/:v/;n' | sed '$ s/$/ \n/' >> $bname.variation.fa
  cat $bname.simulation.fa | seqkit grep -v -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed '$ s/$/ \n/' >> $bname.invariation.fa
  cat $bname.variation.fa $bname.invariation.fa >> $bname.sample.fa
  rm acession.txt $bname.variation.fa $bname.invariation.fa $bname.header.txt $bname.seq.fa $bname.mutation.fa $bname.simulation.fa $bname.fa $bname.vcf $bname.seq.fa.fai
done

for file in work_folder/result/sample/young/macro/*.plant.fa
do
  bname=`basename $file | cut -d. -f1-2`
  grep '>' $file >> $bname.header.txt
  awk 'NR%2==0' $file | tr '\n' '\t' | sed 's/\t/N/g' | sed '1i\>1' >> $bname.seq.fa
  mutation-simulator $bname.seq.fa args -sn 0.00001 -titv 1
  cat $bname.fa | sed '1d' | sed 's/N/\n/g' >> $bname.mutation.fa
  paste -d* $bname.header.txt $bname.mutation.fa | sed 's/*/\n/g' >> $bname.simulation.fa
  diff ../$bname.fa $bname.simulation.fa -y | grep -B 1 '|' | sed  '/--/d' | grep '>' | awk '{print $2}' | awk -F '[> ]' '{print $2}' >> acession.txt
  cat $bname.simulation.fa | seqkit grep -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed 's/$/:v/;n' | sed '$ s/$/ \n/' >> $bname.variation.fa
  cat $bname.simulation.fa | seqkit grep -v -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed '$ s/$/ \n/' >> $bname.invariation.fa
  cat $bname.variation.fa $bname.invariation.fa >> $bname.sample.fa
  rm acession.txt $bname.variation.fa $bname.invariation.fa $bname.header.txt $bname.seq.fa $bname.mutation.fa $bname.simulation.fa $bname.fa $bname.vcf $bname.seq.fa.fai
done

cd work_folder/result/sample/young/macro/mutation
for file in work_folder/result/sample/young/macro/*.mit.fa
do
  bname=`basename $file | cut -d. -f1-3`
  grep '>' $file >> $bname.header.txt
  awk 'NR%2==0' $file | tr '\n' '\t' | sed 's/\t/N/g' | sed '1i\>1' >> $bname.seq.fa
  mutation-simulator $bname.seq.fa args -sn 0.0002 -titv 1
  cat $bname.fa | sed '1d' | sed 's/N/\n/g' >> $bname.mutation.fa
  paste -d* $bname.header.txt $bname.mutation.fa | sed 's/*/\n/g' >> $bname.simulation.fa
  diff ../$bname.fa $bname.simulation.fa -y | grep -B 1 '|' | sed  '/--/d' | grep '>' | awk '{print $2}' | awk -F '[> ]' '{print $2}' >> acession.txt
  cat $bname.simulation.fa | seqkit grep -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed 's/$/:v/;n' | sed '$ s/$/ \n/' >> $bname.variation.fa
  cat $bname.simulation.fa | seqkit grep -v -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed '$ s/$/ \n/' >> $bname.invariation.fa
  cat $bname.variation.fa $bname.invariation.fa >> $bname.sample.fa
  rm acession.txt $bname.variation.fa $bname.invariation.fa $bname.header.txt $bname.seq.fa $bname.mutation.fa $bname.simulation.fa $bname.fa $bname.vcf $bname.seq.fa.fai
done

cd work_folder/result/sample/young/macro/mutation
for file in work_folder/result/sample/young/macro/*.chlo.fa
do
  bname=`basename $file | cut -d. -f1-3`
  grep '>' $file >> $bname.header.txt
  awk 'NR%2==0' $file | tr '\n' '\t' | sed 's/\t/N/g' | sed '1i\>1' >> $bname.seq.fa
  mutation-simulator $bname.seq.fa args -sn 0.0002 -titv 1
  cat $bname.fa | sed '1d' | sed 's/N/\n/g' >> $bname.mutation.fa
  paste -d* $bname.header.txt $bname.mutation.fa | sed 's/*/\n/g' >> $bname.simulation.fa
  diff ../$bname.fa $bname.simulation.fa -y | grep -B 1 '|' | sed  '/--/d' | grep '>' | awk '{print $2}' | awk -F '[> ]' '{print $2}' >> acession.txt
  cat $bname.simulation.fa | seqkit grep -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed 's/$/:v/;n' | sed '$ s/$/ \n/' >> $bname.variation.fa
  cat $bname.simulation.fa | seqkit grep -v -f acession.txt | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed 's/\t/\n/' | sed 's/\t//g' | sed '$ s/$/ \n/' >> $bname.invariation.fa
  cat $bname.variation.fa $bname.invariation.fa >> $bname.sample.fa
  rm acession.txt $bname.variation.fa $bname.invariation.fa $bname.header.txt $bname.seq.fa $bname.mutation.fa $bname.simulation.fa $bname.fa $bname.vcf $bname.seq.fa.fai
done

#merge
a=1
while(($a<=30))
do
    cat work_folder/result/sample/young/mutation/$a.*.fa work_folder/result/sample/young/macro/mutation/$a.*.fa >> $a.fa
    let "a++"
done

#damage
for file in work_folder/result/sample/young/sample/*.fa
do
    bname=`basename $file`
    gargammel/src/deamSim -o $bname.gz -name -mapdamage work_folder/mapdamage/data/youngdata/results_cut/misincorporation.txt double $file
    gzip -d $bname.gz
done

# fa to fq
bbmap/reformat.sh in=simulation_d.fa out1=$bname.simulation.fq qfake=40


####old sample
#fragmentation (Apart from the old_len_dis.txt file, the rest are the same as the yong sample)
#mutation (0.001 for Micro, 0.0001 for Macro and 0.002 for Organelle reads)
#merge
#damage (Apart from the misincorporation.txt file, the rest are the same as the yong sample)
# fa to fq

####deep time sample
#fragmentation (Apart from the veryold_len_dis.txt file, the rest are the same as the yong sample)
#mutation (0.04 for Micro, 0.004 for Macro and 0.08 for Organelle reads)
#merge
#damage (Apart from the misincorporation.txt file, the rest are the same as the yong sample)
# fa to fq




