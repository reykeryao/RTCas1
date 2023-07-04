cd MAFFT_ana
for i in {17..20}
do
        samtools view -F4 RAW/$i.bam |grep -v "XS:i"|cut -f 10|sort |uniq -c|sort -k1,1nr|awk '{print ">ID_"$1"@"$2"\n"$2}' > test.$i.fa  
        mafft --auto --addfragments test.$i.fa --reorder --thread -1 --keeplength --preservecase ../../Ref/R50CCC-RC.fa > $i.output
        fasta_formatter -t -i $i.output -o $i.tmp
        mv $i.tmp $i.output
	rm test.$i.fa
done &

wait
echo "all done!"
cd ..
