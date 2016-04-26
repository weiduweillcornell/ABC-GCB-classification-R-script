# ABC-GCB-classification-R-script

The "ABC-GCB-classifier.zip" contains the R-script and training data(microarray data and sample label) needed to perform the classification. 

To use the code, please download the zip file and modify the working directory to where the unzipped folder is located. Then please supply your to-be-classified data in the ## load RNA-seq data section and approporiately prepare the data for classification. 

The program would generate 3 outputs under the working directory: 

- a heatmap of predictor gene expression level across samples sorted by the probability of belonging to ABC. 

- a plot of probability belonging to ABC/GCB across samples.

- a text file containing the classification label of each sample("ABC", "GCB" or "unclassified").
