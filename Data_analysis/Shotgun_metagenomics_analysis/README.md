# Overview
After the Step1 of classification, the Step2 is to parse taxonomic information from the output files. After this step we could plot figure 1 using Plot_figure1 script.

Step3 is to apply progressive cut-offs on the blast result against fungal database. After this step we could plot figure 2 using the Plot_figure2 script.

Step4 is to assess and select cut-offs for the blast results from publically available datasets. This step is done in the Prepare_Table2_data script.

Step5 is to generate the community composition for each datasets, and compare the similarity of them with the gold standard analysis. This step is mostly done in the Plot_figure4 script.

Step6 is to apply progressive cut-offs on the query coverage and re-calculate the change of similarities between the gold standard classification and the 'best practice' of each dataset. Analysis of this step including plotting figure 4C, 4D and 4E are all summarized in Plot_figure4 script.
