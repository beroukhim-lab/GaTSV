# GaTSV
The GaTSV (Germline and Tumor SV) classifier is an SVM that is able to distinguish germline and somatic structural variants (SVs) in samples with no matched normal. In order to run this, you must have run your WGS samples through the SvABA SV caller. Outputs from other callers may work given that their outputs match the format of SvABA outputs, but the GaTSV classifier is trained on SvABA outputs, so there may be a decrease in performance. Please refer to our paper for more information: https://doi.org/10.1101/2023.10.09.561462

## Figures
In order to recreate the figures from our paper, you can run the code in the `/figures` directory. The code is organized by figure for easy reference. Some data used to generate figures require access to TCGA patient data, so these are not included. Please reach out to us for more information.

## Running the classifier
### In R
All preprocessing and classification scripts are given in the `/scripts` directory. The GaTSV rda object is contained in the `/svm` directory. The `process_classify.R` script will process a given metadata file and SvABA vcf into a bedpe file, and it will classify each variant as germline or somatic. Currently, an example vcf is given within the script, as seen in lines 51 and 52 or the `process_classify.R ` script, so running this code will classify the variants in this vcf. Replacing this file with a similar file containing variants of interest will classify SVs and output a bedpe. 
### Using Docker
The GaTSV docker can be accessed using the pull command `docker pull wchukwu/gatsv_docker:latest`. The syntax for a suitable docker run using the default file mounts to the docker container is `docker run -it -v '\local_path\to\metadata:/data/metadata.txt' -v '\local_path\to\svaba_vcf:/data/input_vcf.vcf' -v '\local_path\to\output_folder:/out/' wchukwu/gatsv_docker:latest /scripts/gaTSV_run.sh sample_name genome cores`. 

An example run is given as `docker run -it -v '$(pwd)\GaTSV\data\example_metadata.txt:/data/metadata.txt' -v '$(pwd)\GaTSV\data\example.sv.vcf:/data/input_vcf.vcf' -v '$(pwd)\GaTSV\out\:/out/' wchukwu/gatsv_docker:latest /scripts/gaTSV_run.sh example.sv hg19 1` 


## Script Dependencies
At the time of the following package versions were used to develop our script. Note, this may not include our figure codes. 

- `BiocGenerics` 0.44.0
- `caTools` 1.18.2 
- `data.table` 1.15.2
- `e1071` 1.7-14
- `GenomeInfoDb` 1.34.9
- `GenomicRanges` 1.50.2
- `gUtils` 0.2.0
- `IRanges` 2.32.0
- `parallel` 4.2.3
- `rlang` 1.1.3
- `ROCR` 1.0-11
- `rstudioapi` 0.15.0
- `S4Vectors` 0.36.2
- `stats4` 4.2.3
- `stringr` 1.5.1
- `here` 1.0.1
- `optparse` 1.7.4
- `rtracklayer` 1.58.0
