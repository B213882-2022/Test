### 
# ArrayExpress数据 + 双端scRNA seq处理
###

home_dir="/home/s2321661/Test"
cd ${home_dir}

# 获取ArrayExpress获取测试的fastq数据（E-MTAB-2600），和GEO不同，下载下来就是fastq文件
mkdir -p ./raw_data/fastq
wget -P ./raw_data/fastq $(awk 'BEGIN{FS="\t"}{if(NR>1){print $28}}' ./E-MTAB-2600.sdrf.txt | sort)

# fastqc + multiqc
fastq_dir="./raw_data/fastq"
target_dir="./raw_data"
mkdir -p ${target_dir}/fastqc
mkdir -p ${target_dir}/multiqc
parallel --verbose -j 8 fastqc -t 16 -o ${target_dir}/fastqc ${fastq_dir}/{} ::: $(ls ${fastq_dir} | sort )
multiqc ${target_dir}/fastqc -o ${target_dir}/multiqc --interactive

# trimmomatic（双端测序）
input_dir="./raw_data/fastq"
# 由于需要用到bash中的变量字符串操作，而parallel无法直接使用，故借助一个函数来执行指令
run_trim(){
tool_dir="/home/s2321661/Dissertation/tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapter_dir="/home/s2321661/Dissertation/tools/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
input_dir="./raw_data/fastq"
output_dir="./trimmed_data"
name1=$1
name2=$2
mkdir -p ${output_dir}/fastq
java -jar ${tool_dir} PE -threads 16 -phred33 ${input_dir}/${name1} ${input_dir}/${name2} ${output_dir}/fastq/${name1/.fastq.gz/_paired.fq.gz} ${output_dir}/fastq/${name1/.fastq.gz/_unpaired.fq.gz} ${output_dir}/fastq/${name2/.fastq.gz/_paired.fq.gz} ${output_dir}/fastq/${name2/.fastq.gz/_unpaired.fq.gz} ILLUMINACLIP:${adapter_dir}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}
export -f run_trim  # 需要将函数export后才能被后续的指令识别（-f表示导出函数）
parallel --verbose -j 8 --xapply run_trim {1} {2} ::: $(ls ${input_dir} | sort | grep "1.fastq.gz") ::: $(ls ${input_dir} | sort | grep "2.fastq.gz")

# 双端测序fastp修剪
tool_dir="/home/s2321661/Dissertation/tools/fastp"
input_dir="./raw_data/fastq"
output_dir="./trimmed_data"
mkdir -p ${output_dir}/fastq
parallel --verbose -j 8 --xapply ${tool_dir} --thread 16 -i ${input_dir}/{1} -o ${output_dir}/fastq/{1} -I ${input_dir}/{2} -O ${output_dir}/fastq/{2} ::: $(ls ${input_dir} | sort | grep "1.fastq.gz") ::: $(ls ${input_dir} | sort | grep "2.fastq.gz")

# fastqc + multiqc
fastq_dir="./trimmed_data/fastq"
target_dir="./trimmed_data"
mkdir -p ${target_dir}/fastqc
mkdir -p ${target_dir}/multiqc
# fastp
parallel --verbose -j 8 fastqc -t 16 -o ${target_dir}/fastqc ${fastq_dir}/{} ::: $(ls ${fastq_dir} | sort )
# trimmomatic
parallel --verbose -j 8 fastqc -t 16 -o ${target_dir}/fastqc {} ::: $(find ${fastq_dir} -name "*_paired.fq.gz" | sort )
# multiqc
# fastp
multiqc ${target_dir}/fastqc -o ${target_dir}/multiqc --interactive  
# trimmomatic（multiqc会读取fastqc结果中的zip文件，这里要避免读取trimmomatic生成的unpaired文件）
multiqc ${target_dir}/fastqc -o ${target_dir}/multiqc --interactive --ignore "*unpaired*"

# 双端测序STAR（索引已在前面建立）
fastq_dir="./trimmed_data/fastq"
run_star(){
tool_dir="/home/s2321661/Dissertation/tools/STAR_2.7.10b/Linux_x86_64_static/STAR"
fastq_dir="./trimmed_data/fastq"
fastq_suffix="_1_paired.fq.gz"
ref_index_dir="/home/s2321661/Dissertation/ref_index_GRCm39"
output_dir="./align"
name1=$1
name2=$2
mkdir -p ${output_dir}/${name1/${fastq_suffix}/}
${tool_dir} --runThreadN 16 --genomeDir ${ref_index_dir} --outFileNamePrefix ${output_dir}/${name1/${fastq_suffix}/}/${name1/${fastq_suffix}/}_ --outSAMtype BAM SortedByCoordinate --readFilesIn ${fastq_dir}/${name1} ${fastq_dir}/${name2} --readFilesCommand zcat
}
export -f run_star
parallel --xapply --verbose -j 8 run_star {1} {2} ::: $(find ${fastq_dir} -name '*1_paired.fq.gz' | sort | sed 's#.*/##') ::: $(find ${fastq_dir} -name '*2_paired.fq.gz' | sort | sed 's#.*/##')

# featureCounts
align_dir="./align"
gtf_dir="/home/s2321661/Dissertation/GRCm39.gtf"
dirs=$(find ${output_dir} -name "*.out.bam" | tr '\n' ' ')
featureCounts -p --countReadPairs -a ${gtf_dir} -o counts.txt -T 60 -t exon -g gene_id ${dirs}

# 最后再用multiqc总结所有内容
multiqc . --ignore "./raw_data/*" --ignore "./trimmed_data/fastqc/*unpaired*" -o multiqc_summary --interactive

# 执行脚本
nohup bash XXX.sh > nohup.log 2>&1 & disown
