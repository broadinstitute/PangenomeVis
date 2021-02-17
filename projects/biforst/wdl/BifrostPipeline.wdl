version 1.0

workflow BifrostPipeline {

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

  
  scatter (index in range(length(fq1_samples))) {
        # String fq_basename=basename(basename(fq1_samples[index], "1.fa.gz"), "_R1.fastq.gz")
        call Concat {
            input:
                fq1 = fq1_samples[index],
                fq2 = fq2_samples[index]
        }
  }

 
  call RunBifrost {
    input:
      references = references,
      samples = Concat.out,
      destination = destination
  }
}

task Concat {
    input {
      File fq1
      File fq2
    }
    String fq_basename=basename(basename(fq1, "1.fa.gz"), "_R1.fastq.gz")
        
    command <<<
      touch ~{fq_basename}.fa
      cat ~{fq1} >> ~{fq_basename}.fa
      cat ~{fq2} >> ~{fq_basename}.fa
      #echo ~{fq_basename}.fa >> fqList.txt
    >>>
  runtime {
    memory: "1 GB"
    docker: "us.gcr.io/broad-gatk/gatk:4.1.7.0"
  }

  output {
      File out = "~{fq_basename}.fa"
    }
}

task RunBifrost {
    input {
        Array[File] references
        Array[File] samples
        String? destination
    }
    command <<<
        touch references.txt
        for ref in ~{sep=' ' references}  ; do
            echo $ref >> references.txt
        done

        cat references.txt

        touch samples.txt
        for sample in ~{sep=' ' samples} ; do
            echo $sample >> samples.txt
        done

        cat samples.txt
        #String fq_basename=basename(basename(fq1_samples[index], "1.fa.gz"), "_R1.fastq.gz")

        #Bifrost build -c -i -d -r references.txt -s samples.txt -o graph
    >>>
  runtime {
    memory: "4 GB"
    docker: "us.gcr.io/broad-dsde-methods/bifrost-snapshot:bifrost"
  }
}

  