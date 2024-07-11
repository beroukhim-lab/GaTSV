FROM rocker/tidyverse:4.2

LABEL AUTHORS='Wolu Chukwu:woluchukwu@broadinstitute.org, Siyun Lee:slee@broadinstitute.org. Alexander Crane: ajc406@case.edu, Shu Zhang: shu@broadinstitute.org' 

RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(version="3.16");' \
&& Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("BiocGenerics")' \
&& Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("GenomeInfoDb")' \
&& Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("GenomicRanges")' \
&& Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("IRanges")' \
&& Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("S4Vectors")' \
&& R -e "install.packages('caTools',dependencies=TRUE, version='1.18.2',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('data.table', dependencies = TRUE, version='1.15.2',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('e1071', dependencies = TRUE, version='1.7-14',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('rlang', dependencies = TRUE, version='1.1.3',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('ROCR', dependencies = TRUE, version='1.0-11',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('rstudioapi', dependencies = TRUE, version='0.15.0',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('stringr', dependencies = TRUE, version='1.5.1',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('optparse',dependencies=TRUE, version='1.7.4',repos='https://cran.rstudio.com')" \
&& R -e "install.packages('here',dependencies=TRUE, version='1.0.1',repos='https://cran.rstudio.com')" \
&& installGithub.r mskilab/gUtils \
#mimic github architecture
&& mkdir -p /scripts \
&& mkdir -p /data \
&& mkdir -p /svm \
&& mkdir -p /out

COPY data/ /data
COPY scripts/ /scripts
COPY svm/ /svm

RUN chmod +x /scripts/gaTSV_run.sh