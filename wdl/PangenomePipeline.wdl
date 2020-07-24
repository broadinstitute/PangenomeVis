version 1.0

import "https://raw.githubusercontent.com/broadinstitute/PangenomeVis/ah-wdl/wdl/MicrobialAlignmentPipeline.wdl" as AlignAndMarkDuplicates


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

  parameter_meta {
    ref_aligned_bam: "Output is aligned duplicate marked coordinate sorted bam."
  }

  scatter (ref in references) {
    call IndexReference { input: ref_fasta = ref }
    
    scatter (fq1 in fq1_samples) {
        String fq_basename=basename(fq1_samples, ".1.fa.gz")
        String fq2 = fq_basename + ".2.ga.gz" 
        call AlignAndMarkDuplicates.MicrobialAlignmentPipeline as AlignToRef {
            input:
                fastq1 = fq1,
                fastq2 = fq2,
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
               input: destination = destination
           }
        }
        
    }
  }
}

task moveOutputs {
    input {
        String? destination
    }
    command <<<
        gsutil cp AlignToRef.output_bam destination
        gsutil cp AlignToRef.output_bam_index destination
    >>>
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
      ls -al
      find . -name *.pac -print
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
  }
}

