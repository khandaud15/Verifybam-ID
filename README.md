# Verifybam-ID
Snakemake based python script of verifybamID  that verifies whether the reads in particular file match previously known genotypes for an individual (or group of individuals), and checks whether the reads are contaminated as a mixture of two samples and it can detect sample contamination.

#### Required tools
* conda install -c bioconda verifybamid
* conda install -c samtools

#### To Run the script on cluster using this command 'modify cluster.json  parameters according to your cluster configuration 
```
snakemake -j 999 --configfile config.yaml --use-conda --nolock --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  -c {cluster.ncpus} -n {cluster.tasks}  -t {cluster.time} --mem {cluster.mem}"
```
