Author: Pascal Schläpfer  
Date: 2021/06/21  
Version: 1.0.0

# Python script to filter vcf files: filter_vcf
This script was produced in accordance of the work involving the manuscript that can be found at https://www.biorxiv.org/content/10.1101/2022.04.13.487913v1.full.pdf. The script searches for an above noise signal (SNPs) in an in-group of samples and filters out signals that are above signal in out-groups. It is a simplification of other methods that use machine learning. For our purpose, Machine learning was not needed, because 1) our cause was known to be monogenetic, and 2) our phenotypes were known with almost 100% certainty. If either of the two conditions are not met, modify the script and use a machine learning algorithm to filter SNPs.

## Function description:
This function parses a vcf file produced by the software freebayes (If you use any other software, you might have to modify the script) containing SNPs and filters SNPs according to parameters given. The function was designed to find allelic SNPs present in an in group sub population of samples, but absent in a out-group sub population. Groups where only a subset of samples have to be hit can also be defined.

## Usage:
python3 filter_vcf.py -in vcf_file_name [-argument argument_value]  
### Arguments:  
#### Mandatory:  
-in [vcf_file]: Input vcf file that should be searched for a SNP variant of interest.  
#### Optional:  
-p [path]: Path where input file is stored. If empty, script uses current working directory to look up input file.  
-out [file]: Output file name, where results should be stored (in tsv format). (default results_filtering_vcf.tsv)  
-rf [path]: Path where results should be stored. If not given, current working directory is used.  
-mtrc [int]: Sets the parameter minimum absolute read coverage of a SNP in a sample of interest. (default is 20)  
-mnrc [int]: Sets the parameter maximum absolute reads of a SNP in an out group sub population, so that if this number of reads is found for the out group sample, it is assumed to be noise and the SNP is still flagged as interesting. (default is 2).  
-mrrc [float between 0 and 1]: Sets the parameter minimum relative read coverage of a variant of a SNP of interest in an in group sub  population. (default is 0.33)  
-mq [float > 0]: Sets the parameter minimum quality a SNP has to reach to be flagged as interesting. (default is 100)  
-i [string]: Defines the samples that should be counted as in group sub population. Sample ids have to be separated by ";".  
-ni [int]: Defines the number of samples given in -i (can be lower than actual number).  
-o [string]: Defines the samples that should be counted as out group sub population. Sample ids have to be separated by ";".  
-no [int]: Defines the number of samples given in -o (can be lower than actual number).  
-fi [string]: Defines the samples that should be counted as optional in group sub population. Sample ids have to be separated by ";".
-nfi [int]: Defines the number of samples given in -fi (can be lower than actual number).  
-fo [string]: Defines the samples that should be counted as optional out group sub population. Sample ids have to be separated by ";".
-nfo [int]: Defines the number of samples given in -fo (can be lower than actual number).
    
#### Auxiliary:
-h: Display this help message  
-v [True, False]: Set verbose to True or False and give screen output.  
-s TME204 or -s TME14 or -s 91-02324 Special mode to reproduce SNP findings of the associated manuscript.  

See publication methods for more information.  
Please use script at own risk. Interprete the results and make changes to parameters or script as you see fit. The author does not take any responsibility for correctness of results of script.
