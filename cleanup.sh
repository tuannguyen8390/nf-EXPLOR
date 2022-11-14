# Bash script to clean-up Nextflow temp file after analysis
rm -rf work/
rm -rf .nextflow/
rm -rf .nextflow.log*
rm -rf slurm-*.out
rm -rf test_results/
