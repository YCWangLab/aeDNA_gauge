#####the analysis based on per reads

# getting reads information table from per reads
cd work_folder/analyse/fa/
for file in *.fa
do
  bname=$(basename "$file" | cut -d. -f1-2)
  grep ">" $file >> ${bname}.header.txt
  sort ${bname}.header.txt >> ${bname}.header.sort.txt
  awk '{ if ($1 ~ /:v/) { print $2, 1 } else { print $2, 0 } }' ${bname}.header.sort.txt > ${bname}.simulation
  awk '{ if ($1 ~ /,/) { print $2, 2 } else if ($1 ~ /_DEAM/) { print $2, 1 } else { print $2, 0 } }' ${bname}.header.sort.txt > ${bname}.damage.txt
  awk -F '[:_]' '{print $6}' ${bname}.header.sort.txt >> ${bname}.length
  paste -d '\t' ${bname}.header.sort.txt ${bname}.length ${bname}.simulation ${bname}.damage.txt > ${bname}.table.csv 
  rm ${bname}.length ${bname}.simulation ${bname}.damage.txt
done 

# sort table 
cd work_folder/analysis/

for file in *.txt
do
  sort -t $'\t' -k1 $file >> $file.table.csv
done

# getting lca information table from per reads
cd work_folder/lca/new/
for file in *.lca
do
  bname=$(basename "$file" | cut -d"." -f1-2)
  
  grep ":v_DEAM" $file | cut -d":" -f1-7 >> 1.txt 
  grep ":v_DEAM" $file | awk -F "\t" '{print $2}' >> 2.txt
  paste -d* 1.txt 2.txt | sed 's/*/\t/g' >> $file.table.txt
  rm 1.txt 2.txt

  grep -v ":v_DEAM" $file | awk -F "\t" '$1 ~ /:v/ {print $1"\t"$2}' | cut -d":" -f1-6 | sed 's/\t$//' >> 1.txt
  grep -v ":v_DEAM" $file | awk -F "\t" '$1 ~ /:v/ {print $1"\t"$2}' | awk -F "\t" '{print $2}' >> 2.txt
  paste -d* 1.txt 2.txt | sed 's/*/\t/g' >> $file.table.txt
  rm 1.txt 2.txt

  grep -v ":v_DEAM" $file | awk -F "\t" '!($1 ~ /:v/) {print $1"\t"$2}' | grep "_DEAM" | cut -d":" -f1-6 >> 1.txt 
  grep -v ":v_DEAM" $file | awk -F "\t" '!($1 ~ /:v/) {print $1"\t"$2}' | grep "_DEAM" | awk -F "\t" '{print $2}' >> 2.txt
  paste -d* 1.txt 2.txt | sed 's/*/\t/g' >> $file.table.txt
  rm 1.txt 2.txt


  sed '1d' $file | grep -v ":v_DEAM" | awk -F "\t" '!($1 ~ /:v/) {print $1"\t"$2}' | grep -v "_DEAM" | cut -d":" -f1-5 >> 1.txt
  sed '1d' $file | grep -v ":v_DEAM" | awk -F "\t" '!($1 ~ /:v/) {print $1"\t"$2}' | grep -v "_DEAM" | awk -F "\t" '{print $2}' >> 2.txt
  paste -d* 1.txt 2.txt | sed 's/*/\t/g' >> $file.table.txt
  rm 1.txt 2.txt
  
done

# sort table 
cd work_folder/analysis/

for file in *.txt
do
  sort -t $'\t' -k1 $file >> $file.table.csv
done

#merge table
install.packages('dplyr')
install.packages('tidyr')
library(dplyr)
library(tidyr)

# 读取文件1
df1 <- read.table("3.veryold.table.sort.csv", header=FALSE, sep="\t")
# 读取文件2
df2 <- read.table("3.veryold.min0max2.lca.table.txt.table.csv", header=FALSE, sep="\t")

# 使用外连接合并，id列作为连接键
merged <- merge(x = df1, y = df2, by.x=1, by.y=1, all = TRUE)

# 使用空格填充缺失值
merged[is.na(merged)] <- 0

# 将结果写入新文件
write.table(merged, "merged.txt", sep = "\t", quote = FALSE, row.names = FALSE)


library(dplyr)
library(tidyr)
# 定义文件名列表
file_names <- c("1.young.table.sort.csv", "1.young.min0max0.lca.table.txt.table.csv", "1.young.min0max1.lca.table.txt.table.csv", "1.young.min0max2.lca.table.txt.table.csv")

# 读取第一个文件，作为合并的基础
merged_data <- read.table(file_names[1], header = FALSE,  sep = "\t", stringsAsFactors = FALSE)

# 循环读取后续文件，并按照第一列合并
for (file_name in file_names[-1]) {
  # 读取文件
  data <- read.table(file_name, header = FALSE, quote="", sep = "\t", stringsAsFactors = FALSE)
  # 按照第一列合并
  merged_data <- merge(merged_data, data, by = 1, all = TRUE)
}
merged_data[is.na(merged_data)] <- 0
# 将结果写入文件
write.table(merged_data, file = "1.young.merged.txt", sep = "\t", quote = FALSE, row.names = FALSE)

"1.modern.table.sort.csv", "1.modern.min0max0.lca.table.txt.table.csv", "1.modern.min0max1.lca.table.txt.table.csv", "1.modern.min0max2.lca.table.txt.table.csv"

# Set 2 and Set 3
cd work_folder/analysis/bam

for file in 1.modern.merged.sorted.bam.sam
do 
  grep "@HD" $file > header_subset.1.txt
  grep "@SQ" $file > header_subset.2.txt
  grep "@PG" $file > header_subset.3.txt
  grep '^@' -v $file > alignment.txt
  awk 'NR==FNR{a[$0]=$0}FNR!=NR{if (!($3 in a)) print $0}' /home/dataset/chenjie/simulation/bowtie2/simulation.acession.txt alignment.txt >> alignment2.txt
  awk 'NR==FNR{a[$0]=$0}FNR!=NR{if (!($3 in a)) print $0}' /home/dataset/chenjie/simulation/bowtie2/simulation.acession2.txt alignment.txt >> alignment3.txt
  rm alignment.txt
  cat header_subset.1.txt header_subset.2.txt header_subset.3.txt alignment2.txt >> 1.modern.merged.sorted.bam2.sam
  samtools view -@ 40 -bS 1.modern.merged.sorted.bam2.sam -o 1.modern.2.bam
  cat header_subset.1.txt header_subset.2.txt header_subset.3.txt alignment3.txt >> 1.modern.merged.sorted.bam3.sam
  samtools view -@ 40 -bS 1.modern.merged.sorted.bam3.sam -o 1.modern.3.bam
  rm header_subset.1.txt header_subset.2.txt header_subset.3.txt alignment2.txt alignment3.txt
done

taxonkit lineage -t simulation.taxid  | taxonkit reformat -F -t | cut -f 1,4,5 >> simulation.tree.txt
awk -F'\t|;' '{print $1, $(NF-2), $(NF-1), $NF}' OFS='\t' simulation.tree.txt >> simulation.fgs.txt
cd /home/dataset/chenjie/simulation/analyse/fa/
cat 1.txaid.txt | while read line 
do
  awk -F'\t' -v line="$line" '$1 == line {print $0}' /home/dataset/chenjie/simulation/bowtie2/simulation.fgs.txt >> 1.fgs.txt
done

paste -d '\t' 1.young.merged.sort.csv.taxid.csv 1.fgs2.txt >> 1.young.merged.sort.csv.taxid.fgs.txt

awk '{ if ($5 == $8) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max0.csv
wc -l 1.young.max0.csv >> 1.young.result.txt
awk '{ if ($5 == $10) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max0.csv
wc -l 1.young.max0.csv >> 1.young.result.txt
awk '{ if ($5 == $9) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max0.csv
wc -l 1.young.max0.csv >> 1.young.result.txt
awk '{ if ($6 == $8) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max1.csv
wc -l 1.young.max1.csv >> 1.young.result.txt
awk '{ if ($6 == $10) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max1.csv
wc -l 1.young.max1.csv >> 1.young.result.txt
awk '{ if ($6 == $9) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max1.csv
wc -l 1.young.max1.csv >> 1.young.result.txt
awk '{ if ($7 == $8) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max2.csv
wc -l 1.young.max1.csv >> 1.young.result.txt
awk '{ if ($7 == $10) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max2.csv
wc -l 1.young.max1.csv >> 1.young.result.txt
awk '{ if ($7 == $9) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max2.csv
wc -l 1.young.max1.csv >> 1.young.result.txt

# species abundance
awk -F "\t" '{print $4}' simulation.fgs.txt | sort -u >> species.taxid.csv
awk -F "\t" 'NR==FNR{a[$1]=$1}FNR!=NR{if ($1 in a) print $0}' species.taxid.csv all_taxa_species.txt >> species.result.txt

awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $0}' species.result.txt species.taxid.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $0}' species.result.txt all_taxa_species.txt >> species.err.csv
111527

# genus abundance
awk -F "\t" '{print $3}' simulation.fgs.txt | sort -u >> genus.taxid.csv
awk -F "\t" 'NR==FNR{a[$1]=$1}FNR!=NR{if ($1 in a) print $0}' genus.taxid.csv all_taxa_genus.txt >> genus.result.txt
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $0}' genus.result.txt genus.taxid.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $0}' genus.result.txt all_taxa_genus.txt >> genus.err.csv
3025755

# family abundance
awk -F "\t" '{print $2}' simulation.fgs.txt | sort -u >> family.taxid.csv
awk -F "\t" 'NR==FNR{a[$1]=$1}FNR!=NR{if ($1 in a) print $0}' family.taxid.csv all_taxa_family.txt >> family.result.txt
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $0}' family.result.txt family.taxid.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $0}' family.result.txt all_taxa_family.txt >> family.err.csv

awk -F "\t" '{print $4}' fgs.txt | sort -n | uniq -c >> species.groundtruth.csv
awk -F "\t" '{print $3}' fgs.txt | sort -n | uniq -c >> genus.groundtruth.csv
awk -F "\t" '{print $2}' fgs.txt | sort -n | uniq -c >> family.groundtruth.csv

grep "@" classified_out >> 1.modern.header.txt
awk -F '[@ ]' '{print $2}' 1.modern.header.txt >> 1.txt
awk -F '[|]' '{print $2}' 1.modern.header.txt >> 2.txt 
paste -d '\t' 1.txt 2.txt | sort -t $'\t' -k1 >> 1.modern.kraken2.csv
sort 1.modern.kraken2.csv | uniq >> 1.modern.kraken2.0.8.csv
rm 1.modern.header.txt 1.txt 2.txt 1.modern.kraken2.csv



Rscript merge.R
sort -t $'\t' -k1 1.modern.kraken2.merged.txt >> 1.modern.kraken2.merged.sort.csv
paste -d '\t' 1.modern.kraken2.merged.sort.csv 1.fgs.txt >> 1.modern.kraken2.merged.sort.table.csv



awk '{ if ($5 == $13) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.2.csv
wc -l 1.modern.0.2.csv >> 1.modern.result.txt
awk '{ if ($5 == $12) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.2.csv
wc -l 1.modern.0.2.csv >> 1.modern.result.txt
awk '{ if ($5 == $11) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.2.csv
wc -l 1.modern.0.2.csv >> 1.modern.result.txt
awk '{ if ($6 == $13) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.4.csv
wc -l 1.modern.0.4.csv >> 1.modern.result.txt
awk '{ if ($6 == $12) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.4.csv
wc -l 1.modern.0.4.csv >> 1.modern.result.txt
awk '{ if ($6 == $11) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.4.csv
wc -l 1.modern.0.4.csv >> 1.modern.result.txt
awk '{ if ($7 == $13) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.6.csv
wc -l 1.modern.0.6.csv >> 1.modern.result.txt
awk '{ if ($7 == $12) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.6.csv
wc -l 1.modern.0.6.csv >> 1.modern.result.txt
awk '{ if ($7 == $11) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.6.csv
wc -l 1.modern.0.6.csv >> 1.modern.result.txt
awk '{ if ($8 == $13) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.8.csv
wc -l 1.modern.0.8.csv >> 1.modern.result.txt
awk '{ if ($8 == $12) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.8.csv
wc -l 1.modern.0.8.csv >> 1.modern.result.txt
awk '{ if ($8 == $11) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.0.8.csv
wc -l 1.modern.0.8.csv >> 1.modern.result.txt
awk '{ if ($9 == $13) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.1.0.csv
wc -l 1.modern.1.0.csv >> 1.modern.result.txt
awk '{ if ($9 == $12) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.1.0.csv
wc -l 1.modern.1.0.csv >> 1.modern.result.txt
awk '{ if ($9 == $11) print }' 1.modern.kraken2.merged.sort.table.csv >> 1.modern.1.0.csv
wc -l 1.modern.1.0.csv >> 1.modern.result.txt

awk '{ if ($5 == $8) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max0.csv
awk '{ sum += $2 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.fragmentation.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $2}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.fragmentation.unmatched.csv
awk '{ sum += $3 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.mutation.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $3}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.mutation.unmatched.csv
awk '{ sum += $4 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.damage.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $4}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.damage.unmatched.csv
awk '{ if ($5 == $10) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max0.csv
awk '{ sum += $2 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.fragmentation.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $2}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.fragmentation.unmatched.csv
awk '{ sum += $3 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.mutation.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $3}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.mutation.unmatched.csv
awk '{ sum += $4 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.damage.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $4}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.damage.unmatched.csv
awk '{ if ($5 == $9) print }' 1.young.merged.sort.csv.taxid.fgs.txt >> 1.young.max0.csv
awk '{ sum += $2 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.fragmentation.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $2}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.fragmentation.unmatched.csv
awk '{ sum += $3 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.mutation.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $3}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.mutation.unmatched.csv
awk '{ sum += $4 } END { avg = sum / NR; print avg }' 1.young.max0.csv >> 1.damage.matched.csv
awk 'NR==FNR{a[$1]=$1}FNR!=NR{if (!($1 in a)) print $4}' 1.young.max0.csv 1.young.merged.sort.csv.taxid.fgs.txt | awk '{ sum += $1 } END { avg = sum / NR; print avg }' >> 1.damage.unmatched.csv


data <- read.csv("1.csv")
colors <- c("blue", "red", "black", "blue", "red", "black", "blue", "red", "black")
pdf("1.pdf")
boxplot(data,
        names = c("Set 1", "Set 1", "Set 1", "Set 2", "Set 2", "Set 2", "Set 3", "Set 3", "Set 3"),
        col = "white",
		border = colors,
        range = 0)
dev.off()


data <- read.csv("2.csv")
colors <- c("blue", "red", "black", "purple", "blue", "red", "black", "purple", "blue", "red", "black", "purple")
pdf("2.pdf")
boxplot(data,
        names = c("No mismatch", "No mismatch", "No mismatch", "No mismatch", "Mismatch 0-1", "Mismatch 0-1", "Mismatch 0-1", "Mismatch 0-1", "Mismatch 0-2", "Mismatch 0-2", "Mismatch 0-2", "Mismatch 0-2"),
        col = "white",
		border = colors,
        range = 0)
dev.off()


data <- read.csv("3.csv")
colors <- c("blue", "red", "black", "purple", "blue", "red", "black", "purple", "blue", "red", "black", "purple")
pdf("3.pdf")
boxplot(data,
        names = c("Species", "Species", "Species", "Species", "Genus", "Genus", "Genus", "Genus", "Family", "Family", "Family", "Family"),
        col = "white",
		border = colors,
        range = 0)
dev.off()

data <- read.csv("4.csv")
colors <- c("blue", "red", "black", "purple", "blue", "red", "black", "purple", "blue", "red", "black", "purple")
pdf("4.pdf")
boxplot(data,
        names = c("Set 1", "Set 1", "Set 1", "Set 1", "Set 2", "Set 2", "Set 2", "Set 2", "Set 3", "Set 3", "Set 3", "Set 3"),
        col = "white",
		border = colors,
        range = 0)
dev.off()




