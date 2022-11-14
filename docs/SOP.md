Steps to run AgVic-UQ CLRC pipeline

1. Obtain & install Nextflow
- Installation guide for Nextflow can be found [here](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. Obtain & install Docker/Shifter/Singularity 
- Installation guide for Docker can be found [here](https://docs.docker.com/get-docker/)
- Installation guide for Shifter can be found [here](https://www.nersc.gov/users/software/nersc-software/shifter/)
- Installation guide for Singularity can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html)
3. Pull assets (genome) and perform some initial setup
- Run the following command to pull assets (genome) and perform some initial setup (choose 1 among Shifter/Docker/Singularity only)
```
nextflow run setup.nf -profile shifter/docker/singularity
```
4. Test run the pipeline ((choose 1 among Shifter/Docker/Singularity only)
```
nextflow run setup.nf -profile shifter/docker/singularity,test
```
5. Run the pipeline
```
nextflow run main.nf -profile shifter/docker/singularity
```
5*. If you run AWS, you can use the following command to run the pipeline
```
nextflow run main.nf -profile shifter/docker/singularity,awsbatch
```

---

The pipeline works using 2 metadata spreadsheet in the `meta` folder: in which:
metadata_SR.csv : metadata for short-read data
metadata_LR.csv : metadata for long-read data

