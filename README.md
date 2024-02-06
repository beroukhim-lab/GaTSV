# GaT-SV
The GaT-SV (Germline and Tumor SV) classifier is an SVM that is able to distinguish germline and somatic structural variants (SVs) in samples with no matched normal. In order to run this, you must have run your WGS samples through the SvABA SV caller. Outputs from other callers may work given that their outputs match the format of SvABA outputs, but the Gat-SV classifier is trained on SvABA outputs, so there may be a decrease in performance. Please refer to our paper for more information: XXX

## Figures
In order to recreate the figures from our paper, you can run the code in the `/figures` directory. All data used to create our figures are provided. The code is organized by figure for easy access. 

## Running the classifier
You can run the code in the `/svm` directory to run our classifier. 