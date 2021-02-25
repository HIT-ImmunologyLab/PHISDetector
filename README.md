# PHISDetector: a great tool to detect and systematically study diverse in silico phage-host interaction signals
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/iwillnever/PHISDector)

## Table of Contents

- [Background](#background)
- [Requirements](requirements)
- [Install](#install)
- [Usage](#usage)
- [Examples](#example)
- [Contributors](#contributors)
- [License](#license)

## Background

Phage-host interactions are appealing systems to study bacterial adaptive evolution and are increasingly recognized as playing an important role in human health and diseases, which may contribute to novel therapeutic agents, such as phage therapy to combat multi-drug resistant infections.PHISDetector is the first comprehensive tool to detect and systematically study diverse in silico phage-host interaction signals, including analyses for oligonucleotide profile/sequence composition, CRISPR-targeting, prophages,alignment-based similarity, protein-protein interactions, special gene check and co-occurrence/co-abundance patterns.Further more,PHISDector provides fancy visualizations for users in [http://www.microbiome-bigdata.com/PHISDetector/index/](http://www.microbiome-bigdata.com/PHISDetector/index/)


## Requirements ##
The source code is written by python3,c++ and wrapped by PyInstaller. Thus it requires python3 and c++ compiler. It works under Linux environment. The following python packages 
are required.
Requirements: 

1. numpy
 
2. Biopython
 
3. joblib

## Development & Funding ##
The development of PHISDector was supported as a collaboration of School of Life Science and Technology, Harbin Institute of Technology.

## Install ##

- step1:Download the whole packages and partial profiles from [https://github.com/HIT-ImmunologyLab/PHISDetector](https://github.com/HIT-ImmunologyLab/PHISDetector)

- step2:Download the rest of large profiles and databases from [http://www.microbiome-bigdata.com/PHISDetector/index/download](http://www.microbiome-bigdata.com/PHISDetector/index/download)


- step3:Unpack the corresponding files downloaded from step2 and put these large profiles and databases in the directory named database in step1.


- step4:When you organize the whole files well,the corresponding directory structure are displayed as shown below. 
![](https://github.com/gancao/PHISDector/tree/master/document/images/standalone_directory_structure.png)


## Usage

The standalone version of PHISDector can predict the hosts of query phages in fasta ,multi-fasta or GenBank format,combining 5 different signals including sequence composition, CRISPR-targeting, prophages, regions of genetic homology and protein-protein or domain-domain interactions.The program of prediction of hosts was divided into two steps.First,detect all the candidate hosts in terms of CRISPR-targeting, prophages and regions of genetic homology,satisfying the **loose** criteria referred in the paper [http://www.microbiome-bigdata.com/PHISDetector/index/](http://www.microbiome-bigdata.com/PHISDetector/index/).Second,the candidate hosts satisfying the strict criteria will be directly output and the rest of candidate hosts will be sent to the machine learning models including RandomForest(RF), Decision Trees (DT), Logistic Regression (LR), and Support Vector Machines (SVM) with RBF kernel and linear kernel, Gaussian Naive Bayes, and Bernoulli Naive Bayes.Finally,the consensus results will be written to a text file and intermediate results of each signal will be saved in corresponding directory.

- step1:Add the PHIS directory named bin into your environment PATH or you can use the absolute path to run the program.
- step2:To run a prediction, you should proceed like the following instructions.The required input is input file path in fasta or multi-fasta or Genbank format and output directory.The optional model input is signal name(all,crispr,prophage,protein_protein_interaction,blast),the default option is all.The visualization html page is generated if you run all signals.If you want to see more parameters,please excute `PHIS` in command line.<br>
In addition,the optional threshold parameters are listed as described below,the program will return the results that meet these thresholds:<br>
for *crispr*: the CRISPR module<br>
(1) min\_mis\_crispr <x> : Minimal mismatch of a Blastn hit on hit spacers(default: 2)<br>
(2) min\_cov\_crispr <x> : Minimal % coverage of a Blastn hit on hit spacers(default: 70)<br>
for *prophage*: the Prophage module<br>
(3) min\_per\_prophage <x> : Minimal % percentage of hit proteins on hit prophage region(default:30)<br>
(4) min\_id\_prophage <x> : Minimal % identity of hit region on hit prophage region by making blastn(default:70)<br>
(5) min\_cov\_prophage <x> : Minimal % coverage of hit region on hit prophage region by making blastn(default:30)<br>
for *protein\_protein\_interaction*: the Protein_Protein_Interaction module<br>
(6) min\_PPI <x> : Minimal PPI number of a pair of phage-host pair(default:1)<br>
(7) min\_DDI <x> : Minimal DDI number of a pair of phage-host pair(default:5)<br>
for *blast*: the GeneHomolog module<br>
(8) min\_per\_blast <x> : Minimal % percentage of hit proteins on query phage and host by making blastp(default:10)<br>
(9) min\_id\_blast <x> : Minimal % identity of hit region on query phage and host by making blastn(default:70)<br>
(10) min\_cov_blast <x> : Minimal % coverage of hit region on query phage and host by making blastn(default:10)
- step3:Predict the hosts of your query phages using all the signals with a simple command:<br>
`PHIS --input <file path> --output <folder name>`
- step4:Predict the hosts of your query phages using single signal with a simple command:
`PHIS --input <file path> --output <folder name> --model all/crispr/prophage/protein_protein_interaction/blast`
- step5:Predict the hosts of your query phages using the criterias you want with a simple example command:<br>
`PHIS --input <file path> --output <folder name> --min_mis_crispr <x> --min_cov_crispr <x>`
`PHIS --input <file path> --output <folder name> --model crispr --min_mis_crispr <x> --min_cov_crispr <x>`

##Examples
You can find a directory named "test" in the PHISDector package. Five test examples, all,crispr,prophage,protein_protein_interaction and blast have been prepared for users.All the example input is "Staphylococcus_phage_47.gb"  in Genbank format.<br>
(1) The folder named "all" contains the prediction results when the input file is Staphylococcus_phage_47.gb. The root directory contains "PHIS_result.html",the visualization page which visualizes and analyzes the prediction hosts of Staphylococcus phage 47 through browser.The file overall_result.txt in results folder lists all the possibly interactive hosts with Staphylococcus phage 47,tab delimited.other files are Intermediate files.There are 13 columns:

1. phage_id: the ID of input phage
2. phage_def: the definition of input phage
3. host_id: the NCBI accession number of the host
4. host_def: the definition of host
5. predict_signal: the signals that support the interaction between the input phage and the host,including crispr,prophage,blast,protein_protein_interaction or sequence_composition 
6. overall_score: PHIS score to rate the interaction between the input phage and the host
7. bayes_bnb: the possibility of the interaction between the input phage and the host based on bayes bnb machine learning model output
8. bayes_gnb: the possibility of the interaction between the input phage and the host based on bayes gnb machine learning model output
9. randomforest_rfc: the possibility of the interaction between the input phage and the host based on randomforest machine learning model output
10. randomforest_dtc: the possibility of the interaction between the input phage and the host based on decision tree machine learning model output
11. svm_rbf: the possibility of the interaction between the input phage and the host based on svm rbf machine learning model output
12. svm_linear: the possibility of the interaction between the input phage and the host based on svm linear machine learning model output
13. logistic: the possibility of the interaction between the input phage and the host based on logistic machine learning model output

![](https://github.com/gancao/PHISDector/tree/master/document/images/PHIS_result_visualization_add_description1.png)

![](https://github.com/gancao/PHISDector/tree/master/document/images/PHIS_result_visualization_add_description2.png)

The folder named "module_results" contains the corresponding results in each kind of signal.

----------

The crispr\_result.txt,tab delimited in the folder named crispr lists all the possibly interactive hosts with Staphylococcus phage 47 based on crispr.There are 8 columns:

1. phage_id: the ID of input phage
2. phage_def: the definition of input phage
3. host_id: the NCBI accession number of the host
4. host_def: the definition of host
5. hit_spacer_num: the number of shared CRISPR spacers between the phage and the host.
6. hit_spacer_identity: the best identity over all the hits between the host spacers and the phage
7. hit_spacer_coverage: the best coverage and best coverage over all the hits between the host spacers and the phage
8. mode: strong/weak,strongly interactive or weakly interactive

----------
The prophage\_result.txt,tab delimited in the folder named prophage lists all the possibly interactive hosts with Staphylococcus phage 47 based on prophage.There are 8 columns:

1. phage_id: the ID of input phage
2. phage_def: the definition of input phage
3. host_id: the NCBI accession number of the host
4. host_def: the definition of host
5. prophage_homolog_percent: the percentage of shared homologous proteins of prophage between the phage and the prophage region based on blastp homologous protein sequences search.
6. prophage_alignment_identity: the best identity over all blastn hits between the host prophage regions and the phage nucleotide sequences
7. prophage_alignment_coverage: the best coverage over all blastn hits between the host prophage regions and the phage nucleotide sequences
8. mode: strong/weak,strongly interactive or weakly interactive


----------
The blast\_result.txt,tab delimited in the folder named blast lists all the possibly interactive hosts with Staphylococcus phage 47 based on genetic homology.There are 8 columns:

1. phage_id: the ID of input phage
2. phage_def: the definition of input phage
3. host_id: the NCBI accession number of the host
4. host_def: the definition of host
5. blastp_score: the similarity between the phage and host based on homologous proteins comparison.
6. blastn_identity: the average identity by merging all blastn hit regions between the phage and host nucleotide sequences.
7. blastn_coverage: the overall coverage by merging all blastn hit regions between the phage and host nucleotide sequences.
8. mode: strong/weak,strongly interactive or weakly interactive

----------
The pro\_int\_result.txt,tab delimited in the folder named pro\_pro\_int lists all the possibly interactive hosts with Staphylococcus phage 47 based on crispr,prophage and genetic homology and corresponding feature values in protein-protein interaction.There are 10 columns:

1. phage_id: the ID of input phage
2. phage_def: the definition of input phage
3. host_id: the NCBI accession number of the host
4. host_def: the definition of host
5. pro_pro_int_num: the number of protein-protein interactions (PPIs) between bacteria proteins and phage proteins 
6. phage_int_pro_prcent: the percentage of phage proteins involved in PPIs
7. bac_int_pro_percent: the percentage of bacterial proteins involved in PPIs
8. pro_pro_int_num_domain: the number of domain-domain interactions (DDIs) between bacteria proteins and phage proteins
9. domain_phage_int_pro_percent: the percentage of phage proteins involved in DDIs
10. domain_bac_int_pro_percent: the percentage of bacterial proteins involved in DDIs

----------
The genehomo\_result.txt,tab delimited in the folder named genehomo lists all the possibly interactive hosts with Staphylococcus phage 47 based on crispr,prophage and genetic homology and corresponding feature values in Oligonucleotide profile/sequence composition.There are 7 columns:

1. phage_id: the ID of input phage
2. phage_def: the definition of input phage
3. host_id: the NCBI accession number of the host
4. host_def: the definition of host
5. virhostmatcher_score: S2* score,the oligonucleotide frequency (ONF) patterns similarity between the input phage andthe host
6. wish_scoreco: WIsH score,the likelihood of the input phage under the Markov models of order k trained on the phage genome. 
7. condon_score: the Euclidean distance of the phage-host pair based on the counts of 64 condon usage profiles in both phage and bacterial coding regions.

(2) The folder named "crispr" contains the prediction results when the input file is Staphylococcus_phage_47.gb,using the crispr signal.the result file "crispr_result.txt" description is described in (1).<br>

(3) The folder named "prophage" contains the prediction results when the input file is Staphylococcus_phage_47.gb,using the prophage signal.the result file "prophage_result.txt" description is described in (1).

(4) The folder named "genetic_homolgy" contains the prediction results when the input file is Staphylococcus_phage_47.gb,using the blast signal.the result file "blast_result.txt" description is described in (1).

(5) The folder named "protein-protein_interaction" contains the prediction results when the input file is Staphylococcus_phage_47.gb,using the protein_protein_interaction signal.the result file "pro_int_result.txt" description is described in (1).

If you want to known more detailed description,please refer to the paper.

## Contributors
This project exists thanks to all the people who contribute.

## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
