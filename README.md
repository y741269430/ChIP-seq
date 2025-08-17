# ChIP-seq

参考：    
[Encode上的教程](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit?tab=t.0#heading=h.9ecc41kilcvq)     
[Encode Github](https://github.com/ENCODE-DCC/chip-seq-pipeline2/tree/master?tab=readme-ov-file)   
[MACS3](https://macs3-project.github.io/MACS/)     

实验大概流程：   
样本甲醛交联——细胞/组织破碎，提取细胞核——超声波打断DNA——用抗体富集目标蛋白（IP）——解交联，DNA提取——片段化，加接头——建库测序。  

生信分析路线：
- 1.当我们数据下机之后，得到的fastq文件。使用`FastQC`软件对raw data进行质量评估。后续clean data同样需要评估。
- 2.使用`fastp`软件对IgG（内参）和IP（目标蛋白）样本的原始数据，进行质控（这一步主要是去除3’端的接头污染、去除重复序列以及低质量序列（保留MPAQ >= 30））。得到clean data。（该软件能否去除重复序列存疑，我之前用的是`Trimmomatic`只有去接头和去除低质量序列的功能。）
- 3.将clean data使用`bowtie2`软件与基因组进行比对，得到的sam文件使用`samtools`转换成bam。
- 4.得到的bam文件，获取其唯一比对以及去重复reads的结果bam文件。
- 5.使用`Deeptools`绘制TSS, Peak center 或GeneBody富集热图（依组学而定），展示数据在这些区域及前后3kb上的富集情况。
- 6.使用`MACS2`或`MACS3`进行peak calling。
- 7.使用`IDR`软件进行样品间高可信度的峰筛选.
- 8.将bam文件转换成bigwig文件，使用`IGV`进行可视化。
- 9.使用r包`ChIPseeker`对peak进行注释。
- 10.使用`homer`或`MEME`进行motif预测。
- 11.使用`MAnorm`（无生物学重复）或`DiffBind`（有生物学重复）进行差异peak分析.

---
## 0.配置环境
跟atac的一致
```bash
mamba activate atac
```

## 1.利用trimmomatic去除接头(Illumina)  
```bash
vim pre_trim.sh
```
```bash
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
```bash
bash pre_trim.sh
```

## 2.Read alignment 双端比对到mm39
```bash
vim ch1_bw2.sh
```
```bash
#!/bin/bash

# Bowtie2 index path
bwt2_idx="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"
nth_bwt2=8  # bowtie2线程数

# 最大并发任务数
MAX_JOBS=4
JOBS=0

# 主循环
while read i; do
(
  #input
  fastq1=trim/${i}_forward_paired.fq.gz
  fastq2=trim/${i}_reverse_paired.fq.gz
  #output
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
```bash
nohup bash ch1_bw2.sh &
```

导出比对率   
```bash
cd logs
grep 'overall alignment rate' *.align.log | \
awk -F'[:%]' 'BEGIN{print "Sample,Alignment Rate (%)"} {gsub(/\.align\.log/, "", $1); printf "%s,%s\n", $1, $2}' > ../alignment_rates.csv
```

## 3. Post-alignment filtering
```bash
vim ch2_sam2bam_rmdup.sh
```
```bash
#!/bin/bash
nth=4  # 线程数
MAX_JOBS=4
JOBS=0

# 主循环
while read i; do
  (
    #input
    bam_file=./bam/${i}.bam
    #output
    bam_FILT=./bam/${i}_FILT.bam
    bam_dupmark=./bam/${i}_dupmark.bam
    bam_FINAL=./bam/${i}_FINAL.bam
    bam_FINAL_nmsrt=./bam/${i}_FINAL_nmsrt.bam
    bam_srt_tmp=./bam/${i}_FILT.srt.tmp.bam
    #logs
    log_file=./logs/${i}.dup.qc
    flagstat_file=./logs/${i}_FINAL.flagstat.qc
    pbc_file=./logs/${i}_FILT.pbc.qc

    # Step 1: 过滤、fixmate、过滤、排序
    samtools view --threads $nth -F 1804 -f 2 -q 30 -u $bam_file | \
    samtools sort --threads $nth -n - | \
    samtools fixmate -r - - | \
    samtools view --threads $nth -F 1804 -f 2 -u - | \
    samtools sort --threads $nth -o $bam_FILT -

    # Step 2: 标记重复
    picard MarkDuplicates \
      INPUT=$bam_FILT \
      OUTPUT=$bam_dupmark \
      METRICS_FILE=$log_file \
      VALIDATION_STRINGENCY=LENIENT \
      ASSUME_SORTED=true \
      REMOVE_DUPLICATES=false

    mv $bam_dupmark $bam_FILT

    # Step 3: 最终筛选、排序、索引
    samtools view --threads $nth -F 1804 -f 2 -b $bam_FILT -o $bam_FINAL
    samtools sort --threads $nth -n $bam_FINAL -o $bam_FINAL_nmsrt
    samtools index $bam_FINAL

    # Step 4: 统计信息
    samtools sort --threads $nth -n $bam_FINAL -O SAM | \
    SAMstats --sorted_sam_file - --outf $flagstat_file

    # Step 5: 计算库复杂度（PBC 指标）
    samtools sort --threads $nth -n $bam_FILT -o $bam_srt_tmp

    bedtools bamtobed -bedpe -i $bam_srt_tmp | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
    grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > $pbc_file
    
    rm $bam_srt_tmp
    rm $bam_FILT

  ) &

  ((JOBS++))
  if [[ "$JOBS" -ge "$MAX_JOBS" ]]; then
    wait -n
    ((JOBS--))
  fi
done < filenames

wait

```
```bash
nohup bash ch2_sam2bam_rmdup.sh
```

## 4. MACS3 call peak
分别写脚本分开运行
```bash
#!/bin/bash
## peak calling (macs3) ##

cat con_file | while read i; 
do
nohup macs3 callpeak -f BAMPE \
-t ./bam/${i}_FINAL.bam \
-c ./bam/E1_FINAL.bam ./bam/E2_FINAL.bam ./bam/E3_FINAL.bam \
-g mm -n ./macs3/${i} -B -q 0.1 --broad-cutoff 0.1 & 
done
```
```bash
#!/bin/bash
## peak calling (macs3) ##

cat conp_file | while read i; 
do
nohup macs3 callpeak -f BAMPE \
-t ./bam/${i}_FINAL.bam \
-c ./bam/E4_FINAL.bam ./bam/E6_FINAL.bam \
-g mm -n ./macs3/${i} -B -q 0.1 --broad-cutoff 0.1 & 
done
```

## 5. Convert PE BAM to tagAlign (BED 3+3 format)
```bash
vim ch3_bam2bed.sh
```
```bash
#!/bin/bash

cat filenames | while read i; 
do
    #input
    input_bam=./bam/${i}_FINAL_nmsrt.bam
    #output
    output_bedpe_gz=./bed/${i}_FINAL_nmsrt.bedpe.gz
    output_ta_file=./bed/${i}_FINAL_TA_FILE.bed

    nohup bedtools bamtobed -bedpe -mate1 -i $input_bam | gzip -nc > $output_bedpe_gz && \
    zcat $output_bedpe_gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > $output_ta_file &
    
done

```
```bash
bash ch3_bam2bed.sh
```



## 以下内容暂时先不做   
---
## 6.Generate self-pseudoreplicates for each replicate (PE datasets)

```bash
vim ch_pse.sh
```
```bash
#!/bin/bash
MAX_JOBS=4
JOBS=0

# 主循环
while read i; do
  (
    #input
    input_ta_file=./bed/${i}_FINAL_TA_FILE.bed
    #output
    temp_bedpe_file=./bed/${i}_temp.bedpe
    filt_nodup_prefix=./bed/${i}.filt.nodup
    pr1_tagalign_file=./bed/${i}.filt.PE2SE.pr1.tagAlign.gz
    pr2_tagalign_file=./bed/${i}.filt.PE2SE.pr2.tagAlign.gz

    # Step 1: Create fake BEDPE from tagAlign
    zcat $input_ta_file | sed 'N;s/\n/\t/' > $temp_bedpe_file

    # Step 2: Count read pairs
    nlines=$(wc -l < $temp_bedpe_file)
    nlines=$(( (nlines + 1) / 2 ))

    # Step 3: Shuffle and split
    shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat $input_ta_file | wc -c) -nosalt </dev/zero 2>/dev/null) $temp_bedpe_file \
    | split -d -l ${nlines} - $filt_nodup_prefix

    # Step 4: Convert to tagAlign format
    awk 'BEGIN{OFS="\t"}{
      printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",
      $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
    }' ${filt_nodup_prefix}00 | gzip -nc > $pr1_tagalign_file

    awk 'BEGIN{OFS="\t"}{
      printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",
      $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
    }' ${filt_nodup_prefix}01 | gzip -nc > $pr2_tagalign_file

    # 清理临时文件
    rm -f $temp_bedpe_file ${filt_nodup_prefix}00 ${filt_nodup_prefix}01
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





