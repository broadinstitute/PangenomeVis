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
    g=0

    mkdir -p $haps
    mkdir -p $reads

    ((echo '>ref_hap') && (cat flank_left.fasta ${ref_file} flank_right.fasta | grep -v '>' | sed '1{N;s/\n//;}' | sed '1{N;s/\n//;}')) >  $haps/ref_hap.fasta && samtools faidx $haps/ref_hap.fasta

    for sn in A B
    do
        ~/repositories/readsim/readsim -g $g -r $haps/ref_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.imperfect
        ((g=g+1))

        mccortex31 build -f -m 1G -k 31 -S -s $sn -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz ${asm}/${sn}.raw.ctx
        mccortex31 clean -f -m 1G -S -o ${asm}/${sn}.clean.ctx ${asm}/${sn}.raw.ctx
        mccortex31 inferedges -f -m 1G -o ${asm}/${sn}.inferedges.ctx ${asm}/${sn}.raw.ctx
    done

    if [[ -f "$alt1_file" ]]; then
        ((echo '>alt1_hap') && (cat flank_left.fasta ${alt1_file} flank_right.fasta | grep -v '>' | sed '1{N;s/\n//;}' | sed '1{N;s/\n//;}')) >  $haps/alt1_hap.fasta && samtools faidx $haps/alt1_hap.fasta

        for sn in C D
        do
            ~/repositories/readsim/readsim -g $g -r $haps/alt1_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.imperfect
            ((g=g+1))

            mccortex31 build -f -m 1G -k 31 -S -s $sn -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz ${asm}/${sn}.raw.ctx
            mccortex31 clean -f -m 1G -S -o ${asm}/${sn}.clean.ctx ${asm}/${sn}.raw.ctx
            mccortex31 inferedges -f -m 1G -o ${asm}/${sn}.inferedges.ctx ${asm}/${sn}.raw.ctx
        done
    fi

    if [[ -f "$alt2_file" ]]; then
        ((echo '>alt2_hap') && (cat flank_left.fasta ${alt2_file} flank_right.fasta | grep -v '>' | sed '1{N;s/\n//;}' | sed '1{N;s/\n//;}')) >  $haps/alt2_hap.fasta && samtools faidx $haps/alt2_hap.fasta

        for sn in E F
        do
            ~/repositories/readsim/readsim -g $g -r $haps/alt2_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.imperfect
            ((g=g+1))

            mccortex31 build -f -m 1G -k 31 -S -s $sn -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz ${asm}/${sn}.raw.ctx
            mccortex31 clean -f -m 1G -S -o ${asm}/${sn}.clean.ctx ${asm}/${sn}.raw.ctx
            mccortex31 inferedges -f -m 1G -o ${asm}/${sn}.inferedges.ctx ${asm}/${sn}.raw.ctx
        done
    fi

    if [[ -f "$alt3_file" ]]; then
        ((echo '>alt3_hap') && (cat flank_left.fasta ${alt3_file} flank_right.fasta | grep -v '>' | sed '1{N;s/\n//;}' | sed '1{N;s/\n//;}')) >  $haps/alt3_hap.fasta && samtools faidx $haps/alt3_hap.fasta

        for sn in G H
        do
            ~/repositories/readsim/readsim -g $g -r $haps/alt3_hap.fasta -l $rl -d $d -e 0.001 ${reads}/sample_${sn}.imperfect
            ((g=g+1))

            mccortex31 build -f -m 1G -k 31 -S -s $sn -2 ${reads}/sample_${sn}.imperfect.1.fa.gz:${reads}/sample_${sn}.imperfect.2.fa.gz ${asm}/${sn}.raw.ctx
            mccortex31 clean -f -m 1G -S -o ${asm}/${sn}.clean.ctx ${asm}/${sn}.raw.ctx
            mccortex31 inferedges -f -m 1G -o ${asm}/${sn}.inferedges.ctx ${asm}/${sn}.raw.ctx
        done
    fi

    for allele_type in ref alt1 alt2 alt3
    do
        fasta="${haps}/${allele_type}_hap.fasta"

        if [[ -f "$fasta" ]]; then
            for sn in A B C D E F G H
            do
                reads1="${reads}/sample_${sn}.imperfect.1.fa.gz"
                reads2="${reads}/sample_${sn}.imperfect.2.fa.gz"

                if [[ -f "$reads1" ]]; then
                    ~/repositories/minimap2/minimap2 -ayYL --MD --eqx -x sr -R '@RG\tSM:${sn}\tID:${allele_type}_imperfect' $fasta ${reads1} ${reads2} | samtools sort | samtools view -b > ${reads}/sample_${sn}.${allele_type}_imperfect.bam && samtools index ${reads}/sample_${sn}.${allele_type}_imperfect.bam
                fi
            done
        fi
    done

    mccortex31 join -f -m 4G -o ${asm}/all.joined.ctx ${asm}/*.inferedges.ctx
    mccortex31 view -k ${asm}/all.joined.ctx > ${asm}/all.joined.txt
done
