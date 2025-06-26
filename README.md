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

# Bowtie2 index path
bwt2_idx="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"
nth_bwt2=4  # bowtie2线程数

# 最大并发任务数
MAX_JOBS=4
JOBS=0

# 主循环
while read i; do
(

  fastq1=trim/${i}_forward_paired.fq.gz
  fastq2=trim/${i}_reverse_paired.fq.gz
  bam_file=./bam/${i}.bam
  log_file=./logs/${i}.align.log
  flagstat_qc=./logs/${i}.flagstat.qc

  # Step 1: 比对 + 转 BAM + 排序
  bowtie2 -X2000 --mm --threads $nth_bwt2 -x $bwt2_idx \
    -1 $fastq1 -2 $fastq2 2> $log_file | \
    samtools view -Su - | samtools sort -o $bam_file -

  # Step 2: flagstat 样本质控
  samtools sort -n --threads $nth_bwt2 $bam_file -O SAM | \
    SAMstats --sorted_sam_file - --outf $flagstat_qc

) &

((JOBS++))
if [[ "$JOBS" -ge "$MAX_JOBS" ]]; then
  wait -n         # 等待任一任务结束
  ((JOBS--))      # 减去已结束的任务
fi

done < filenames

wait
```

## 3. Sam转bam，然后去重复
```bash
vim ch2_sam2bam_rmdup.sh

#!/bin/bash
nth=4  # 线程数
MAX_JOBS=4
JOBS=0

# 主循环
while read i; do
  (
    # Step 1: 过滤、fixmate、过滤、排序
    samtools view --threads $nth -F 1804 -f 2 -q 30 -u ./bam/${i}.bam | \
    samtools sort --threads $nth -n - | \
    samtools fixmate -r - - | \
    samtools view --threads $nth -F 1804 -f 2 -u - | \
    samtools sort --threads $nth -o ./bam/${i}_FILT.bam -

    # Step 2: 标记重复
    picard MarkDuplicates \
      INPUT=./bam/${i}_FILT.bam \
      OUTPUT=./bam/${i}_dupmark.bam \
      METRICS_FILE=./logs/${i}.dup.qc \
      VALIDATION_STRINGENCY=LENIENT \
      ASSUME_SORTED=true \
      REMOVE_DUPLICATES=false

    mv ./bam/${i}_dupmark.bam ./bam/${i}_FILT.bam

    # Step 3: 最终筛选、排序、索引
    samtools view --threads $nth -F 1804 -f 2 -b ./bam/${i}_FILT.bam -o ./bam/${i}_FINAL.bam
    samtools sort --threads $nth -n ./bam/${i}_FINAL.bam -o ./bam/${i}_FINAL_nmsrt.bam
    samtools index ./bam/${i}_FINAL.bam

    # Step 4: 统计信息
    samtools sort --threads $nth -n ./bam/${i}_FINAL.bam -O SAM | \
    SAMstats --sorted_sam_file - --outf ./logs/${i}_FINAL.flagstat.qc

    # Step 5: 计算库复杂度（PBC 指标）
    samtools sort --threads $nth -n ./bam/${i}_FILT.bam -o ./bam/${i}_FILT.srt.tmp.bam

    bedtools bamtobed -bedpe -i ./bam/${i}_FILT.srt.tmp.bam | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
    grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > ./logs/${i}_FILT.pbc.qc

    rm ./bam/${i}_FILT.srt.tmp.bam
    rm ./bam/${i}_FILT.bam

  ) &

  ((JOBS++))
  if [[ "$JOBS" -ge "$MAX_JOBS" ]]; then
    wait -n
    ((JOBS--))
  fi
done < filenames

wait

```

## 4. bam2bed
```bash
vim ch3_bam2bed.sh

#!/bin/bash

cat filenames | while read i; 
do
nohup bedtools bamtobed -bedpe -mate1 -i ./bam/${i}_FINAL_nmsrt.bam | gzip -nc > ./bed/${i}_FINAL_nmsrt.bedpe.gz && \
zcat ./bed/${i}_FINAL_nmsrt.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ./bed/${i}_FINAL_TA_FILE.bed &

done
```

##
```bash
#!/bin/bash
MAX_JOBS=4
JOBS=0

# 主循环
while read i; do
  (
# Step 1: Create fake BEDPE from tagAlign
zcat ./bed/${i}_FINAL_TA_FILE.bed | sed 'N;s/\n/\t/' > ./bed/${i}_temp.bedpe

# Step 2: Count read pairs
nlines=$(wc -l < ./bed/${i}_temp.bedpe)
nlines=$(( (nlines + 1) / 2 ))

# Step 3: Shuffle and split
shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat ./bed/${i}_FINAL_TA_FILE.bed | wc -c) -nosalt </dev/zero 2>/dev/null) ./bed/${i}_temp.bedpe \
  | split -d -l ${nlines} - ./bed/${i}.filt.nodup

# Step 4: Convert to tagAlign format
awk 'BEGIN{OFS="\t"}{
  printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",
  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
}' ./bed/${i}.filt.nodup00 | gzip -nc > ./bed/${i}.filt.PE2SE.pr1.tagAlign.gz

awk 'BEGIN{OFS="\t"}{
  printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",
  $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
}' ./bed/${i}.filt.nodup01 | gzip -nc > ./bed/${i}.filt.PE2SE.pr2.tagAlign.gz

rm -f ./bed/${i}_temp.bedpe ./bed/${i}.filt.nodup00 ./bed/${i}.filt.nodup01
) &

  ((JOBS++))
  if [[ "$JOBS" -ge "$MAX_JOBS" ]]; then
    wait -n
    ((JOBS--))
  fi
done < filenames

wait
```



## 计算单端
Trim R1 fastq to 50bp
```bash
python trimfastq.py $FASTQ_R1 50 | gzip -nc >  trim_50/${i}_forward_paired.fq.gz
```
```bash
#!/bin/bash

# Bowtie2 index path
bwt2_idx="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"
nth_bwt2=4  # bowtie2线程数

# 最大并发任务数
MAX_JOBS=4
JOBS=0

# 主循环
while read i; do
(

  fastq_file=trim_50/${i}_forward_paired.fq.gz
  bam_file=./bam/${i}_fastq_r1.bam
  filt_bam_file=./bam/${i}_fastq_r1_FILT.bam
  log_file=./logs/${i}_fastq_r1.log

  # 1. 比对 + BAM 排序
  bowtie2 --mm -x $mm39 --threads $nth_bwt2 -U $fastq_file 2> $log_file | \
    samtools view -Su - | samtools sort -o $bam_file -

  # 2. Filter BAM
  samtools view -F 1804 --threads $nth_bwt2 -q 30 -b $bam_file -o $filt_bam_file

) &

((JOBS++))
if [[ "$JOBS" -ge "$MAX_JOBS" ]]; then
  wait -n         # 等待任一任务结束
  ((JOBS--))      # 减去已结束的任务
fi

done < filenames

wait
```









## 安装phantompeakqualtools
[phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools)

```bash
git clone https://github.com/kundajelab/phantompeakqualtools
cd phantompeakqualtools
R
```
```r
install.packages("snow", repos="http://cran.us.r-project.org")
install.packages("snowfall", repos="http://cran.us.r-project.org")
install.packages("bitops", repos="http://cran.us.r-project.org")
install.packages("caTools", repos="http://cran.us.r-project.org")
install.packages("BiocManager")
BiocManager::install("Rsamtools")
devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)

library(snow)
library(snowfall)
library(bitops)
library(caTools)
library(Rsamtools)
library(spp)
q()
```





