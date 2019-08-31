shell.prefix("source ~/.bash_profile; set -euo pipefail;")

#wildcard rule for samples  directory
SAMPLES, = glob_wildcards(config['datadirs']['bam'] + "/" + "{sample}_Aligned.sortedByCoord.out.md.bam")

rule all:
     input:
         expand(config['datadirs']['bai'] + "/" + "{sample}_Aligned.sortedByCoord.out.md.bam.bai", sample=SAMPLES),
         expand(config['datadirs']['verifybam'] + "/" + "{sample}.selfSM", sample=SAMPLES),
         expand(config['datadirs']['parse'] + "/" + "{sample}-results.txt", sample=SAMPLES),
         config['datadirs']['verifybamFinal'] + "/" +  "verify.txt"


rule index:
     input:
        bam = config['datadirs']['bam'] + "/" + "{sample}_Aligned.sortedByCoord.out.md.bam"
     output:
        bai = config['datadirs']['bai'] + "/" + "{sample}_Aligned.sortedByCoord.out.md.bam.bai"
     shell: """
          samtools index {input.bam} {ouput.bai}
          """

rule verifybam_id:
   input:
       vcf = config['reference']['vcf'][['hg38'],
       bam = config['datadirs']['bam'] + "/" + "{sample}_Aligned.sortedByCoord.out.md.bam".
       index = config['datadirs']['bai'] + "/" + "{sample}_Aligned.sortedByCoord.out.md.bam.bai"
   output:
       selfSM = config['datadirs']['verifybam'] + "/" + "{sample}.selfSM",
       bestSM = config['datadirs']['verifybam'] + "/" + "{sample}.bestSM",
       depthSM = config['datadirs']['verifybam'] + "/" + "{sample}.depthSM"
   params:
       prefix = config['datadirs']['verifybam'] + "/" + "{sample}",
       individual = "{sample}"
   resources:
       mem_mb = 30000
   shell: """
        verifyBamID --vcf {input.vcf} --bam {input.bam} --best --ignoreRG --smID {params.individual} --out {params.prefix} --verbose
         """     

 rule parse_verify:
    input: 
       selfSM = config['datadirs']['verifybam'] + "/" + "{sample}.selfSM",
       bestSM = config['datadirs']['verifybam'] + "/" + "{sample}.bestSM",
       depthSM = config['datadirs']['verifybam'] + "/" + "{sample}.depthSM"
    output: config['datadirs']['parse'] + "/" + "{sample}-results.txt"
    run:
        selfSM = open(input.selfSM, "rt")
        bestSM = open(input.bestSM, "rt")
        depthSM = open(input.depthSM, "rt")
        results = open(output[0], "w")

        for line in selfSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[0] == "#SEQ_ID" and cols[2] == "CHIP_ID" and \
                       cols[3] == "#SNPS" and cols[4] == "#READS" and \
                       cols[5] == "AVG_DP" and cols[6] == "FREEMIX" and \
                       cols[11] == "CHIPMIX", "selfSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                seq_id = cols[0]
                chip_id = cols[2]
                chip_id == seq_id
                snps = cols[3]
                reads = cols[4]
                avg_dp = cols[5]
                selfSM_freemix = cols[6]
                selfSM_chipmix = cols[11]

        for line in bestSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[0] == "#SEQ_ID" and cols[2] == "CHIP_ID" and \
                       cols[3] == "#SNPS" and cols[4] == "#READS" and \
                       cols[5] == "AVG_DP" and cols[6] == "FREEMIX" and \
                       cols[11] == "CHIPMIX", "bestSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                seq_id = cols[0]
                bestSM_chip_id = cols[2]
                match = str(bestSM_chip_id == seq_id).upper()
                bestSM_freemix = cols[6]
                bestSM_chipmix = cols[11]

         # Report the number of SNPs that had more than a minimum read depth
        depth_min = 10
        snps_w_min = 0
        for line in depthSM:
            # Confirm the header columns
            if line[0] == "#":
                cols = line.strip().split("\t")
                assert cols[1] == "DEPTH" and cols[2] == "#SNPs", \
                       "depthSM columns are as expected"
            else:
                cols = line.strip().split("\t")
                depth = int(cols[1])
                n_snps = int(cols[2])
                if depth >= depth_min:
                    snps_w_min = snps_w_min + n_snps

        out_header = ["id", "best", "match",
                      "self_freemix", "self_chipmix",
                      "best_freemix", "best_chipmix",
                      "snps", "reads", "avg_dp",
                      "min_dp", "snps_w_min"]
        out_cols = [seq_id, bestSM_chip_id, match,
                    selfSM_freemix, selfSM_chipmix,
                    bestSM_freemix, bestSM_chipmix,
                    snps, reads, avg_dp,
                    str(depth_min), str(snps_w_min)]
        results.write("\t".join(out_header) + "\n")
        results.write("\t".join(out_cols) + "\n")

        selfSM.close()
        bestSM.close()
        depthSM.close()
        results.close()

rule combine_verify:
    input: 
       files = expand(config['datadirs']['parse'] + "/" + "{sample}-results.txt", sample = SAMPLES)
    output: config['datadirs']['verifybamFinal'] + "/" +  "verify.txt"
    shell:
        "head -n 1 {input[0]} > {output};"
        "cat {input.files} | grep -v \"id\" | sort -k1n >> {output}"                 
