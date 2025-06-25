# ChIP-seq

## 目录 ####
- 0.配置环境
- 1.利用trimmomatic去除接头(Illumina)
- 2.比对到mm39
- 3.Sam转bam，然后处理bam文件
- 4.去重复

---
## 0.配置环境
跟atac的一致
```bash
mamba activate atac
```

## 1.利用trimmomatic去除接头(Illumina)  
```bash
vim pre_trim.sh

#!/bin/bash
## trimmomatic ##

cat filenames | while read i; 
do
nohup trimmomatic PE -phred33 -threads 4 \
./RawData/${i}/${i}*_R1_001.fastq.gz \
./RawData/${i}/${i}*_R2_001.fastq.gz \
./trim/${i}_forward_paired.fq.gz \
./trim/${i}_forward_unpaired.fq.gz \
./trim/${i}_reverse_paired.fq.gz \
./trim/${i}_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &

done
```

## 2.比对到mm39 
```bash
vim ch1_bw2.sh

#!/bin/bash
## Alignment to mm39 ##

mm39="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"

cat filenames | while read i; 
do
nohup bowtie2 -p 4 --very-sensitive -X 2000 \
-x ${mm39} \
-1 ./trim/${i}*_forward_paired.fq.gz \
-2 ./trim/${i}*_reverse_paired.fq.gz \
-S ./bam/${i}.sam 2> ./bam/${i}_map.txt & 
done
```

## 3. Sam转bam，然后处理bam文件
```bash
vim ch2_sam2bam.sh

#!/bin/bash

th=10  # 线程数设定
MAX_JOBS=4  # 最多允许同时运行的样本数
JOBS=0

# 读取文件名列表（每一行一个样本名）
while read i; do
  (
    # Step 1: 将SAM文件转换为BAM并排序（position-sorted）
    samtools view --threads ${th} -Su ./bam/${i}.sam | \
    samtools sort -o ./bam/${i}_sorted.bam

    # Step 2: 输出flagstat QC信息（name-sorted）
    samtools sort -n --threads ${th} ./bam/${i}_sorted.bam -O SAM | \
    SAMstats --sorted_sam_file - --outf ./bam/${i}_flagstat_qc

    # Step 3: 过滤unmapped、低质量、多重比对等reads，生成name-sorted中间BAM
    samtools view -F 1804 -f 2 -q 30 -u ./bam/${i}_sorted.bam | \
    samtools sort -n -o ./bam/${i}_sorted_name.bam

    # Step 4: 修复paired-end配对信息
    samtools fixmate -r ./bam/${i}_sorted_name.bam ./bam/${i}.fixmate.tmp

    # Step 5: 再次过滤并生成最终position-sorted BAM
    samtools view -F 1804 -f 2 -u ./bam/${i}.fixmate.tmp | \
    samtools sort -o ./bam/${i}.bam

    # Step 6: 删除中间临时文件

    rm ./bam/${i}.fixmate.tmp
    rm ./bam/${i}_sorted_name.bam
  ) > ./bam/log_${i}.log 2>&1 &

((JOBS++))
if [[ "$JOBS" -ge "$MAX_JOBS" ]]; then
  wait -n  # 等其中一个任务完成
  ((JOBS--))
fi

done < filenames

wait  # 等所有剩余任务完成

```

## 4.标记重复片段   
```bash
#!/bin/bash
## Mark duplicates ##

cat filenames | while read i; 
do
nohup picard MarkDuplicates \
  INPUT=./bam/${i}.bam \
  OUTPUT=./bam/${i}_dedup.bam \
  METRICS_FILE=./bam/${i}_dup.qc \
  VALIDATION_STRINGENCY=LENIENT \
  ASSUME_SORTED=true \
  REMOVE_DUPLICATES=false &
done
```

## 5.对标记重复的片段，进行处理
```bash
#!/bin/bash
MAX_JOBS=4  # 最多允许同时运行的样本数
JOBS=0

while read i; do
  (
  # Step 1: Filter BAM (valid, proper-paired, primary alignment)
  samtools view -F 1804 -f 2 -b ./bam/${i}_dedup.bam > ./bam/${i}_FINAL.bam

  # Step 2: Name-sort for BEDPE conversion
  samtools sort -n ./bam/${i}_FINAL.bam -o ./bam/${i}_FINAL.nmsrt.bam

  # Step 3: Index the filtered BAM
  samtools index ./bam/${i}_FINAL.bam

  # Step 4: Generate flagstat QC report
  samtools sort -n --threads 10 ./bam/${i}_FINAL.bam -O SAM | \
    SAMstats --sorted_sam_file - --outf ./bam/${i}_FINAL.flagstat.qc

  # Step 5: Calculate PBC metrics
  bedtools bamtobed -bedpe -i ./bam/${i}_FINAL.nmsrt.bam | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
    grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} 
         ($1==1){m1=m1+1} 
         ($1==2){m2=m2+1} 
         {m0=m0+1} 
         {mt=mt+$1} 
         END{
           printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt, m0, m1, m2, m0/mt, m1/m0, m1/m2
         }' > ./bam/${i}_FINAL.pbc.qc

  # Step 6: Generate BEDPE file
  bedtools bamtobed -bedpe -mate1 -i ./bam/${i}_FINAL.nmsrt.bam | \
    gzip -nc > ./bam/${i}_FINAL.nmsrt.bam.bedpe.gz

  # Step 7: Convert BEDPE to tagAlign format
  zcat ./bam/${i}_FINAL.nmsrt.bam.bedpe.gz | \
    awk 'BEGIN{OFS="\t"}{
      printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",
      $1,$2,$3,$9,$4,$5,$6,$10
    }' | gzip -nc > ./bam/${i}.FINAL_TA
  ) &

((JOBS++))
if [[ "$JOBS" -ge "$MAX_JOBS" ]]; then
  wait -n  # 等其中一个任务完成
  ((JOBS--))
fi

done < filenames

wait  # 等所有剩余任务完成

```


