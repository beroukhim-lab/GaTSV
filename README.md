# GaTSV
The GaTSV (Germline and Tumor SV) classifier is an SVM that is able to distinguish germline and somatic structural variants (SVs) in samples with no matched normal. In order to run this, you must have run your WGS samples through the SvABA SV caller. Outputs from other callers may work given that their outputs match the format of SvABA outputs, but the GaTSV classifier is trained on SvABA outputs, so there may be a decrease in performance. Please refer to our paper for more information: https://pubmed.ncbi.nlm.nih.gov/40081367/

## Figures
In order to recreate the figures from our paper, you can run the code in the `/figures` directory. The code is organized by figure for easy reference. Some data used to generate figures require access to TCGA patient data, so these are not included. Please reach out to us for more information.

## Running the classifier
### In R
All preprocessing and classification scripts are given in the `/scripts` directory. The GaTSV rda object is contained in the `/svm` directory. The `process_classify.R` script will process a given metadata file and SvABA vcf into a bedpe file, and it will classify each variant as germline or somatic. Currently, an example vcf is given within the script, as seen in lines 51 and 52 or the `process_classify.R ` script, so running this code will classify the variants in this vcf. Replacing this file with a similar file containing variants of interest will classify SVs and output a bedpe. 
### Using Docker
The GaTSV docker can be accessed using the pull command `docker pull wchukwu/gatsv_docker:latest`. The syntax for a suitable docker run using the default file mounts to the docker container is `docker run -it -v '\local_path\to\metadata:/data/metadata.txt' -v '\local_path\to\svaba_vcf:/data/input_vcf.vcf' -v '\local_path\to\output_folder:/out/' wchukwu/gatsv_docker:latest /scripts/gaTSV_run.sh sample_name genome cores`. 

An example run is given as `docker run -it -v '$(pwd)\GaTSV\data\example_metadata.txt:/data/metadata.txt' -v '$(pwd)\GaTSV\data\example.sv.vcf:/data/input_vcf.vcf' -v '$(pwd)\GaTSV\out\:/out/' wchukwu/gatsv_docker:latest /scripts/gaTSV_run.sh example.sv hg19 1` 

## Expected Outcomes
1. The GaTSV classifier generates two output files per sample. These are:
- `[sample name]_processed.bedpe`: This file contains SVs that passed our internal QC metrics (as described in our manuscript) and the feature annotations.
- `[sample name]_classified.bedpe`: This file contains the class label for each classified SV under the `predicted_class` column. An SV is labeled as either `GERMLINE` or `SOMATIC`.

NOTE: The number of SVs in `[sample name]_processed.bedpe` may not match the number of SVs in `[sample name]_classified.bedpe` because we require that SVs that are not translocations must be at least 1000bp in length before they can be confidently classified by the GaTSV classifier. SVs that do not meet these criteria will not be in `[sample name]_classified.bedpe`.

2. These outputs will be stored in `C:\Users\your_profile\outputs` or in `/path/to/outputs`.

3. These results can be compared with the TCGA and external pHGG dataset results from our manuscript. Although the specific features can vary at an individual level, we can typically expect over a 10:1 ratio of germline to somatic events, the SPAN of somatic SVs to be much larger than germline events on average, and most of the germline events to be deletions, while somatic events are more evenly distributed across all SV types. As mentioned previously, these analyses were conducted on a population level, so individual SVs may not always follow these trends. The following figures were generated from the GaTSV calls on the pHGG dataset:



<div align="center">
  <img src="https://u.cubeupload.com/snlee12/Screenshot20251215at.png" height="400" alt="SV Type Pie Chart">
  <img src="https://u.cubeupload.com/snlee12/9a7Screenshot20251215at.png" height="400" alt="SV Count Distribution">
</div>


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
