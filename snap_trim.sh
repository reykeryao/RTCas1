for i in 1 6 11
do
	cutadapt -g ^NNNGCAATAATCTATACAATACAACACATACAAACAAATTCTTAAGGTAAAAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 2 7 12 17 18 19 20
do
	cutadapt -g ^NNNGCAATAATCTATACAATACAACACATACAAACAAATTCTTAAGGTCCCAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 3 8 13 16
do
	cutadapt -g ^NNNGCAATAATCTATACAATACAACACATACAAACAAATTCTTAAGGTGGGAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 4 9 14
do
	cutadapt -g ^NNNGCAATAATCTATACAATACAACACATACAAACAAATTCTTAAGGTTTTAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 5 10 15
do
	cutadapt -g ^NNNGCAATAATCTATACAATACAACACATACAAACAAATTCTTAAGGTCCGCAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 21 27
do
	cutadapt -g ^NNNAACACATACAAACAAATTCTTAAGGTCCCAAAAAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 22 28
do
	cutadapt -g ^NNNCAACACATACAAACAAATTCTTAAGGTCCCAAAAAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 23 29
do
	cutadapt -g ^NNNACAACACATACAAACAAATTCTTAAGGTCCCAAAAAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 24 30
do
	cutadapt -g ^NNNTACAACACATACAAACAAATTCTTAAGGTCCCAAAAAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 25 31
do
	cutadapt -g ^NNNACAACACATACAAACAAATTCTTAAGGTCCCAAAAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done

for i in 26 32
do
	cutadapt -g ^NNNTACAACACATACAAACAAATTCTTAAGGTCCCAAAA --discard-untrimmed -o $i.fq.gz ../joined/$i.fq.gz -j 12
done







