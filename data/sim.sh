#!/bin/bash

set -euxo pipefail

FILES=$(find variants -name 'ref_allele.fasta')

for ref_file in $FILES
do
    alt1_file=$(echo $ref_file | sed 's/ref_/alt1_/')
    alt2_file=$(echo $ref_file | sed 's/ref_/alt2_/')
    alt3_file=$(echo $ref_file | sed 's/ref_/alt3_/')

    bn=$(basename $ref_file | sed 's/.fasta//')
    dn=$(dirname $ref_file | sed 's/refs/reads/')
    haps="${dn}/haps"
    reads="${dn}/reads"
    asm="${dn}/asm"
    rl=75
    d=50

    mkdir -p $haps
    mkdir -p $reads

    ((echo '>ref_hap') && (paste -d'' flank_left.fasta ${ref_file} flank_right.fasta | grep -v '>')) > $haps/ref_hap.fasta && samtools faidx $haps/ref_hap.fasta

    for sn in A B
    do
        ~/repositories/readsim/readsim -r $haps/ref_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.ref_imperfect
        ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:ref_imperfect' $haps/ref_hap.fasta ${reads}/sample_${sn}.ref_imperfect.1.fa.gz ${reads}/sample_${sn}.ref_imperfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.ref_imperfect.bam && samtools index ${reads}/sample_${sn}.ref_imperfect.bam

        mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.ref_imperfect.1.fa.gz:${reads}/sample_${sn}.ref_imperfect.2.fa.gz ${asm}/sample_${sn}.ref_imperfect.ctx
        mccortex31 clean -f -m 1G -o ${asm}/sample_${sn}.ref_imperfect.cleaned.ctx ${asm}/sample_${sn}.ref_imperfect.ctx
        mccortex31 inferedges -f -m 1G -o ${asm}/sample_${sn}.ref_imperfect.inferedges.ctx ${asm}/sample_${sn}.ref_imperfect.cleaned.ctx
        mccortex31 thread -f -m 4G -W -2 ${reads}/sample_${sn}.ref_imperfect.1.fa.gz:${reads}/sample_${sn}.ref_imperfect.2.fa.gz -o ${asm}/sample_${sn}.ref_imperfect.inferedges.ctp.gz ${asm}/sample_${sn}.ref_imperfect.inferedges.ctx 
    done

    if [[ -f "$alt1_file" ]]; then
        ((echo '>alt1_hap') && (paste -d'' flank_left.fasta ${alt1_file} flank_right.fasta | grep -v '>')) > $haps/alt1_hap.fasta && samtools faidx $haps/alt1_hap.fasta

        for sn in C D
        do
            ~/repositories/readsim/readsim -r $haps/alt1_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.alt1_imperfect
            ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:alt1_imperfect' $haps/alt1_hap.fasta ${reads}/sample_${sn}.alt1_imperfect.1.fa.gz ${reads}/sample_${sn}.alt1_imperfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.alt1_imperfect.bam && samtools index ${reads}/sample_${sn}.alt1_imperfect.bam

            mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.alt1_imperfect.1.fa.gz:${reads}/sample_${sn}.alt1_imperfect.2.fa.gz ${asm}/sample_${sn}.alt1_imperfect.ctx
            mccortex31 clean -f -m 1G -o ${asm}/sample_${sn}.alt1_imperfect.cleaned.ctx ${asm}/sample_${sn}.alt1_imperfect.ctx
            mccortex31 inferedges -f -m 1G -o ${asm}/sample_${sn}.alt1_imperfect.inferedges.ctx ${asm}/sample_${sn}.alt1_imperfect.cleaned.ctx
            mccortex31 thread -f -m 4G -W -2 ${reads}/sample_${sn}.alt1_imperfect.1.fa.gz:${reads}/sample_${sn}.alt1_imperfect.2.fa.gz -o ${asm}/sample_${sn}.alt1_imperfect.inferedges.ctp.gz ${asm}/sample_${sn}.alt1_imperfect.inferedges.ctx 
        done
    fi

    if [[ -f "$alt2_file" ]]; then
        ((echo '>alt2_hap') && (paste -d'' flank_left.fasta ${alt2_file} flank_right.fasta | grep -v '>')) > $haps/alt2_hap.fasta && samtools faidx $haps/alt2_hap.fasta

        for sn in E F
        do
            ~/repositories/readsim/readsim -r $haps/alt2_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.alt2_imperfect
            ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:alt2_imperfect' $haps/alt2_hap.fasta ${reads}/sample_${sn}.alt2_imperfect.1.fa.gz ${reads}/sample_${sn}.alt2_imperfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.alt2_imperfect.bam && samtools index ${reads}/sample_${sn}.alt2_imperfect.bam

            mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.alt2_imperfect.1.fa.gz:${reads}/sample_${sn}.alt2_imperfect.2.fa.gz ${asm}/sample_${sn}.alt2_imperfect.ctx
            mccortex31 clean -f -m 1G -o ${asm}/sample_${sn}.alt2_imperfect.cleaned.ctx ${asm}/sample_${sn}.alt2_imperfect.ctx
            mccortex31 inferedges -f -m 1G -o ${asm}/sample_${sn}.alt2_imperfect.inferedges.ctx ${asm}/sample_${sn}.alt2_imperfect.cleaned.ctx
            mccortex31 thread -f -m 4G -W -2 ${reads}/sample_${sn}.alt2_imperfect.1.fa.gz:${reads}/sample_${sn}.alt2_imperfect.2.fa.gz -o ${asm}/sample_${sn}.alt2_imperfect.inferedges.ctp.gz ${asm}/sample_${sn}.alt2_imperfect.inferedges.ctx 
        done
    fi

    if [[ -f "$alt3_file" ]]; then
        ((echo '>alt3_hap') && (paste -d'' flank_left.fasta ${alt3_file} flank_right.fasta | grep -v '>')) > $haps/alt3_hap.fasta && samtools faidx $haps/alt3_hap.fasta

        for sn in G H
        do
            ~/repositories/readsim/readsim -r $haps/alt3_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.alt3_imperfect
            ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:alt3_imperfect' $haps/alt3_hap.fasta ${reads}/sample_${sn}.alt3_imperfect.1.fa.gz ${reads}/sample_${sn}.alt3_imperfect.2.fa.gz | samtools sort | samtools view -b > ${reads}/sample_${sn}.alt3_imperfect.bam && samtools index ${reads}/sample_${sn}.alt3_imperfect.bam

            mccortex31 build -f -m 1G -k 31 -S -s ${sn} -2 ${reads}/sample_${sn}.alt3_imperfect.1.fa.gz:${reads}/sample_${sn}.alt3_imperfect.2.fa.gz ${asm}/sample_${sn}.alt3_imperfect.ctx
            mccortex31 clean -f -m 1G -o ${asm}/sample_${sn}.alt3_imperfect.cleaned.ctx ${asm}/sample_${sn}.alt3_imperfect.ctx
            mccortex31 inferedges -f -m 1G -o ${asm}/sample_${sn}.alt3_imperfect.inferedges.ctx ${asm}/sample_${sn}.alt3_imperfect.cleaned.ctx
            mccortex31 thread -f -m 4G -W -2 ${reads}/sample_${sn}.alt3_imperfect.1.fa.gz:${reads}/sample_${sn}.alt3_imperfect.2.fa.gz -o ${asm}/sample_${sn}.alt3_imperfect.inferedges.ctp.gz ${asm}/sample_${sn}.alt3_imperfect.inferedges.ctx 
        done
    fi

    mccortex31 join -f -m 4G -S -o ${asm}/final_joined.imperfect.ctx ${asm}/sample_*imperfect.inferedges.ctx
    mccortex31 view -k ${asm}/final_joined.imperfect.ctx > ${asm}/final_joined.imperfect.txt
done
