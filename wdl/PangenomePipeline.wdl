version 1.0

import "https://raw.githubusercontent.com/broadinstitute/PangenomeVis/master/wdl/MicrobialAlignmentPipeline.wdl" as AlignAndMarkDuplicates


workflow PangenomePipeline {

  meta {
    description: "Align all samples with all reads"
  }

  input {
    Array[File] references
    Array[File] fq1_samples
    Array[File] fq2_samples
    String? destination

    #Optional runtime arguments
    Int? preemptible_tries
  }

  scatter (ref in references) {
    String ref_basename = basename(basename(basename(basename(ref, ".gz"), "sta"), ".fa"), ".fna")
    call IndexReference { input: ref_fasta = ref }
    
    scatter (index in range(length(fq1_samples))) {
        String fq_basename=basename(basename(fq1_samples[index], "1.fa.gz"), "_R1.fastq.gz")
        call AlignAndMarkDuplicates.MicrobialAlignmentPipeline as AlignToRef {
            input:
                fastq1 = fq1_samples[index],
                fastq2 = fq2_samples[index],
                basename = fq_basename + "_" + ref_basename,
                ref_dict = IndexReference.ref_dict,
                ref_fasta = ref,
                ref_fasta_index = ref + ".fai",
                ref_amb = IndexReference.ref_amb,
                ref_ann = IndexReference.ref_ann,
                ref_bwt = IndexReference.ref_bwt,
                ref_pac = IndexReference.ref_pac,
                ref_sa = IndexReference.ref_sa,
                preemptible_tries = preemptible_tries
        }   

        if (defined(destination)) {
           call moveOutputs {
               input: 
                  destination = destination,
                  bam = AlignToRef.aligned_bam,
                  bai = AlignToRef.aligned_bai
           }
        }
        
    }
  }
}

task moveOutputs {
    input {
      File bam
      File bai
      String? destination
    }
    command <<<
        gsutil cp ~{bam} ~{destination}/
        gsutil cp ~{bai} ~{destination}/
    >>>
  runtime {
    memory: "1 GB"
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
  }
}


task IndexReference {
  input {
    File ref_fasta    
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(ref_fasta, "GB") * 2.5) + 20
  String basename = basename(ref_fasta)
  
  command <<<
      set -e
      cp ~{ref_fasta} .
      /usr/gitc/bwa index ~{basename}
      java -jar /usr/gitc/picard.jar CreateSequenceDictionary REFERENCE=~{ref_fasta} OUTPUT=~{basename}.dict
      
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }

  output {
    File ref_amb = "~{basename}.amb"
    File ref_ann = "~{basename}.ann"
    File ref_bwt = "~{basename}.bwt"
    File ref_pac = "~{basename}.pac"
    File ref_sa = "~{basename}.sa"
    File ref_dict = "~{basename}.dict"
  }
}

