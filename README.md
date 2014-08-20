README

The package is used to identify identify the joint effect of multiple Copy Number Variants (CNVs) on categorical phenotypes. This version of code is used to demonstrate how the result in the paper is computed. Please contact the authors of the paper to report bugs and any other problems.


---- Requirement ----

To run the programs in this package, you need to correctly install those packages.

-Perl 5.0 or newer
-Matlab 2010 or newer
-SLEP 4.0 (http://www.public.asu.edu/~jye02/Software/SLEP/download.htm) 

---- Installation ----

1. Untar all the files in the same folder.
2. Add the home path of SLEP and all its sub-folders in the searching path of Matlab.
3. Run test.pl to test the settings.

---- Running the program and input/output files ----

[1. Generating the candidate CNV groups by random walking]

The script, genonet.GetSubnetRandom.iter.pl, generates a set of candidate groups by random walking on a given network/graph. Each group contains two or more CNVs. The maximum size of group can be set by the parameter MaxStep. 

>genonet.GetSubnetRandom.iter.pl [NodeListFile] [EdgeListFile] [MapSeqGene] [SeqList] [MaxStep] [Prefix]

NodeListFile: This file defines all the node names in the given network. Each line of this file contains the name of one node. Node name can be a combination of any letters, numbers and "_". Node name is case sensitive. 

EdgeListFile: This file defines the given network by listing all the edges and their weights in the network. Each line in this file defines an edge with the following format:

[Node name 1] [Node name 2] [weight]  

The fields of [Node name 1] and [Node name 2] should appear in the NodeListFile. The field of [weight] is a real positive number, which is the weight of the edge between the two nodes. These three fields are separated by one or more spaces.

MapSeqGene: This file defines the mapping from CNV to genes by listing one CNV ID (or CNV location) and one gene name in each line. In our paper, we map CNVs to genes by one base pair overlap. 

SeqList: This file lists all the CNVs mapped (genic CNVs or gCNVs), which are mapped to some genes in the given network. The line number of each copy number variant is used as index of the corresponding CNV in the output files.

MaxStep: This parameter is a natural number, which indicates the maximal number of steps of a random walk on the network.

Prefix: This parameter is a string indicating the prefix of output files. Non-alphabet characters are not allowed in this string.


The program output two files, the [Prefix].idx.txt and [Prefix].map.txt. They contain all the information about the candidate groups. These two files are used in the following group LASSO analysis. 

[2. Group LASSO analysis]

In this step, we perform group LASSO analysis on a given set of data with the candidate groups generated in the first step. The main matlab script can be called in the following format.

>> run_cv_pvalue_group_iter_maxoccu(datafile,groupIdxFile,groupMapFile)

Input:
	datafile: This parameter is the file name of the input data file. The input data file is a table, where the first column is phenotypic label and the remaining columns are CNV features. Each row represents a single sample of the input data. The value of each cell in the table represents the genotypes (integer copy number) of the corresponding CNV in the sample.
	groupIdxFile: The file name of the .idx.txt file generated in the first step.
	groupMapFile: The file name of the .map.txt file generated in the first step.

Output: 
	The script, run_cv_pvalue_group_iter_maxoccu.m, outputs three files for each label, when the number of labels are larger than 2. Each set of 3 files are associated with the resultant groups distinguishing one label from other labels. The three files are named in the following way.

The [Prefix].[label].nzgroup.csv contains the indices of groups with non-zero weights, which are the CNV groups selected by our model.

The [Prefix].[label].pvalue.csv contains the p-value of each CNV.

The [Prefix].[label].weight.csv contains the weight of each CNV.

To obtain the decent sparsity of the model, you can change the range of the parameter of muRange in the script of run_cv_pvalue_group_iter_maxoccu.m . To compute p-value more precisely, you can increase the variable of nrepeat in the script of run_gllassopvalue.m. The recommended value is 1000, NOT the default value.


[3. Post analysis]

In this step, we use the script, checkgroup.pl, to translate the result from CNV index numbers to CNV IDs(Locations) and gene names. 

The usage of checkgroup.pl is as following.

> checkgroup.pl [Prefix].[label].nzgroup.csv

The script outputs groups and their members' CNV IDs/Locations and gene names to standard output. It also writes the network edges between CNVs in each group to a file named as [Prefix].[label].nzgroup.csv.group.[group index].edges.csv . In these files, each line contains an edge in the network that lists the two gene names and the weight between them. 

Besides CNV IDs/Locations in selected CNV groups, this script also outputs the Ensemble IDs and refSeq IDs of the corresponding genes.



