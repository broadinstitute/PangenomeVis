#!/bin/bash

set -euxo pipefail

FILES=$(find variants -name 'ref_allele.fasta')

for ref_file in $FILES
do
    alt_file=$(echo $ref_file | sed 's/ref_/alt_/')
    bn=$(basename $ref_file | sed 's/.fasta//')
    dn=$(dirname $ref_file | sed 's/refs/reads/')
    haps="${dn}/haps"
    reads="${dn}/reads"
    asm="${dn}/asm"
    rl=75
    d=50

    mkdir -p $haps
    ((echo '>ref_hap') && (paste -d'' flank_left.fasta ${ref_file} flank_right.fasta | grep -v '>')) > $haps/ref_hap.fasta
    ((echo '>alt_hap') && (paste -d'' flank_left.fasta ${alt_file} flank_right.fasta | grep -v '>')) > $haps/alt_hap.fasta

    samtools faidx $haps/ref_hap.fasta
    samtools faidx $haps/alt_hap.fasta

    mkdir -p $reads
    for sn in A B C
    do
        ~/repositories/readsim/readsim -r $haps/ref_hap.fasta -l $rl -d $d ${reads}/sample_${sn}.perfect
        ~/repositories/readsim/readsim -r $haps/ref_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.imperfect

        ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:perfect' $haps/ref_hap.fasta ${reads}/sample_${sn}.perfect.1.fa.gz ${reads}/sample_${sn}.perfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.perfect.bam
        ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:imperfect' $haps/ref_hap.fasta ${reads}/sample_${sn}.imperfect.1.fa.gz ${reads}/sample_${sn}.imperfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.imperfect.bam

        samtools index ${reads}/sample_${sn}.perfect.bam
        samtools index ${reads}/sample_${sn}.imperfect.bam

        mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.perfect.1.fa.gz:${reads}/sample_${sn}.perfect.2.fa.gz ${asm}/sample_${sn}.perfect.ctx

        mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz ${asm}/sample_${sn}.imperfect.ctx
        mccortex31 clean -f -m 1G -o ${asm}/sample_${sn}.imperfect.cleaned.ctx ${asm}/sample_${sn}.imperfect.ctx
        mccortex31 inferedges -f -m 1G -o ${asm}/sample_${sn}.imperfect.inferedges.ctx ${asm}/sample_${sn}.imperfect.cleaned.ctx
        mccortex31 thread -f -m 4G -W -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz -o ${asm}/sample_${sn}.imperfect.inferedges.ctp.gz ${asm}/sample_${sn}.imperfect.inferedges.ctx 
    done
    for sn in D E
    do
        ~/repositories/readsim/readsim -r $haps/alt_hap.fasta -l $rl -d $d ${reads}/sample_${sn}.perfect
        ~/repositories/readsim/readsim -r $haps/alt_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.imperfect

        ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:perfect' $haps/ref_hap.fasta ${reads}/sample_${sn}.perfect.1.fa.gz ${reads}/sample_${sn}.perfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.perfect.bam
        ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:imperfect' $haps/ref_hap.fasta ${reads}/sample_${sn}.imperfect.1.fa.gz ${reads}/sample_${sn}.imperfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.imperfect.bam

        samtools index ${reads}/sample_${sn}.perfect.bam
        samtools index ${reads}/sample_${sn}.imperfect.bam

        mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.perfect.1.fa.gz:${reads}/sample_${sn}.perfect.2.fa.gz ${asm}/sample_${sn}.perfect.ctx

        mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz ${asm}/sample_${sn}.imperfect.ctx
        mccortex31 clean -f -m 1G -o ${asm}/sample_${sn}.imperfect.cleaned.ctx ${asm}/sample_${sn}.imperfect.ctx
        mccortex31 inferedges -f -m 1G -o ${asm}/sample_${sn}.imperfect.inferedges.ctx ${asm}/sample_${sn}.imperfect.cleaned.ctx
        mccortex31 thread -f -m 4G -W -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz -o ${asm}/sample_${sn}.imperfect.inferedges.ctp.gz ${asm}/sample_${sn}.imperfect.inferedges.ctx 
    done

    mccortex31 join -f -m 4G -S -o ${asm}/final_joined.perfect.ctx ${asm}/sample_A.perfect.ctx ${asm}/sample_B.perfect.ctx ${asm}/sample_C.perfect.ctx ${asm}/sample_D.perfect.ctx ${asm}/sample_E.perfect.ctx
    mccortex31 join -f -m 4G -S -o ${asm}/final_joined.imperfect.ctx ${asm}/sample_A.imperfect.inferedges.ctx ${asm}/sample_B.imperfect.inferedges.ctx ${asm}/sample_C.imperfect.inferedges.ctx ${asm}/sample_D.imperfect.inferedges.ctx ${asm}/sample_E.imperfect.inferedges.ctx

    mccortex31 view -k ${asm}/final_joined.perfect.ctx > ${asm}/final_joined.perfect.txt
    mccortex31 view -k ${asm}/final_joined.imperfect.ctx > ${asm}/final_joined.imperfect.txt
done
