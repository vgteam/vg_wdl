version 1.0

task indexReference {
    input {
        File in_reference_file
        Int in_index_mem = 4
        Int in_index_disk = 2 * round(size(in_reference_file, "G")) + 10
    }

    command <<<
        set -eux -o pipefail
        
        ln -s ~{in_reference_file} ref.fa
                
        # Index the subset reference
        samtools faidx ref.fa 
        
        # Save a reference copy by making the dict now
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
          R=ref.fa \
          O=ref.dict
    >>>
    output {
        File reference_index_file = "ref.fa.fai"
        File reference_dict_file = "ref.dict"
    }
    runtime {
        preemptible: 2
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task indexVcf {
    input {
        File in_vcf
        Int in_index_mem = 4
        Int in_index_disk = 2 * round(size(in_vcf, "G")) + 10
    }

    command <<<
        set -eux -o pipefail
        bcftools index -t -o variants.vcf.gz.tbi ~{in_vcf}
    >>>
    output {
        File vcf_index_file = "variants.vcf.gz.tbi"
    }
    runtime {
        preemptible: 2
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task fixVCFContigNaming {
    input {
        File in_vcf
        String in_prefix_to_strip = "GRCh38."
        Int in_index_mem = 4
        Int in_index_disk = 2 * round(size(in_vcf, "G")) + 10
    }

    command <<<
    set -eux -o pipefail

    bcftools view -h ~{in_vcf} | grep contig | awk -v PREF="~{in_prefix_to_strip}" '{split($0, a, ","); gsub("##contig=<ID=", "", a[1]); chr=a[1]; gsub(PREF, "", a[1]); print chr"\t"a[1]}' > chrsrename.txt
    
    bcftools annotate --rename-chrs chrsrename.txt -o variants.vcf.gz -Oz ~{in_vcf}
    >>>
    output {
        File vcf = "variants.vcf.gz"
    }
    runtime {
        preemptible: 2
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

