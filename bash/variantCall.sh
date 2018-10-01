#!bin/bash

#sorting 

cd raw/test

for f in $(ls *.bam);
do ~/samtools/bin/samtools sort $f -o ${f/.bam/_s.bam};
done 

#splitting forward/reverse reads 

for f in $(ls *_s.bam);
do ~/samtools/bin/samtools view -F 0x10 -h $f | ~/samtools/bin/samtools view -bS - >  ${f/_s.bam/_sf.bam};
done

for f in $(ls *_s.bam);
do ~/samtools/bin/samtools view -f 0x10 -h $f | ~/samtools/bin/samtools view -bS - >  ${f/_s.bam/_sr.bam};
done

#generate consensus seq

for f in $(ls *_s.bam);
do ~/samtools/bin/samtools mpileup -d 40000 -uf ../seq/2001Ia $f | ~/bcftools/bin/bcftools call -cv -Oz -o ${f/_s.bam/.gz};
done

for f in $(ls *.gz);
do ~/htslib/bin/tabix $f;
done

for f in $(ls *.gz);
do cat ../seq/2001Ia | ~/bcftools/bin/bcftools consensus -e '(DP4[0]+DP4[1])>(DP4[2]+DP4[3])' $f > ${f/.gz/.fas};
done

#variant call by lofreq

for f in $(ls *_s.bam);
do ~/lofreq_star_v2.1.2/bin/lofreq call -f ${f/_s.bam/.fas} -o ${f/_s.bam/.vcf} $f;
done

for f in $(ls *_sf.bam);
do ~/lofreq_star_v2.1.2/bin/lofreq call -f ${f/_sf.bam/.fas} -o ${f/_sf.bam/_f.vcf} $f;
done

for f in $(ls *_sr.bam);
do ~/lofreq_star_v2.1.2/bin/lofreq call -f ${f/_sr.bam/.fas} -o ${f/_sr.bam/_r.vcf} $f;
done

