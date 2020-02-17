AbundanceQt.py: Python script for classifying lipid measurements according to lipid class. For users new to Python, a setup guide is provided separately as 'Setup guide.pdf'. 


The script takes in two input files, 'ClassSpecs.txt' and 'data.txt' (hard-coded file names). ClassSpecs.txt contains upper and lower bounds for retention times and mass as well as the response factor for each lipid class in a tabular format where rows correspond to lipid classes and the columns to each of the aforementioned quantities. 'Data.txt' contains the lipid measurements, with a header line describing the contents of each column, and one row per compound. For each compound, columns contain the compound ID, mass, retention time (1 column each), and a variable number of columns corresponding to the sample measurements. Templates for each of 'ClassSpecs.txt' and 'data.txt' are provided in this repository.



Running the script provides five output files:

-compounds-abundances.txt: Similar to 'data.txt', but with each compound labeled by lipid class in the 4th column

-compounds-abundances-RF-corrected: Similar to compounds-abundances.txt, but measured abundances are corrected according to the response factor

-LC-abundances.txt: Sum of abundances for all compounds in each lipid class for each sample

-LC-abundances-RF-corrected.txt: Sum of abundances for all compounds in each lipid class for each sample, corrected according to the response factor

-compoundCounts.txt: Number of compounds from each lipid class present in each sample
