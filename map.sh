bowtie2 -p 12 -a --very-sensitive-local -L 5 --seed 740714 -x ../Ref/$2 -U $1.fq.gz 2>RAW/$1.log | samtools view -bS - > RAW/$1.bam
