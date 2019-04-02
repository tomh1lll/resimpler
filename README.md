A simple script for finding differences in statistics for a particular set of genes and the background surrounding genes.

Input should be a tab delimited bed file with a header in the following format:
    
    chr	start	stop	class	group	STAT
    X 100 500 mel AMP NA

Columns must be named as above, gene name is not required.
If a stat is going to be resampled, its column name should be called with -v
If a group column is present and multiple groups (e.g. species, populations) are found within the column, the script will iterate over each group seperately.
class column should contain the classes you are planning on finding differences for. Any classes not defined under -c will be considered part of the background class.

here is an example script for the provided test input, assuming that all 4 statistics (-v snipre,snipreb,alpha,dos) want to be resampled for AMPs and Antiviral RNAi (-c AMP,Antiviral RNAi) for 1000 bootstraps each (-b 1000). In this example, we'll only include background genes within 100kb of the focal gene (-w 100000), but won't care about size differences (-s not used)

    python resimpler.py -b 1000 -i test_input.txt -o test_output.txt -v snipre,snipreb,alpha,dos -c AMP,Antiviral RNAi

More statistics can be assessed by increasing the list, or including a text file with columns that should be resampled over. Classes of genes can also be included as a text file.
