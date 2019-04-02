A simple script for finding differences in statistics for a particular set of genes and the background surrounding genes.

Input should be a tab delimited bed file with a header in the following format:
    
    chr	start	stop	class	group	STAT
    X 100 500 mel AMP NA

Columns must be named as above, gene name is not required.
If a stat is going to be resampled, its column name should be called with -v
If a group column is present and multiple groups (e.g. species, populations) are found within the column, the script will iterate over each group seperately.
class column should contain the classes you are planning on finding differences for. Any classes not defined under -c will be considered part of the background class.

for example, the test input provided:
    group	chr	start	end	class	snipre	snipreb	alpha	dos
    mel	3R	26042206	26042660	AMP	-0.386202589	-1.156999555	-2	-0.25
    mel	2R	10634867	10635720	AMP	-0.243477792	-1.024424195	0	0
    mel	2L	3539247	3548423	AMP	0.044868257	-1.545892965	0.75	0
    mel	2R	14272072	14272416	AMP	-0.225929233	-1.049411862	-0.5	-0.1
    mel	2R	14272329	14273094	AMP	-0.741004016	-0.932324151	-5.111111111	-0.285714286
    mel	2R	14754896	14755400	AMP	-0.18498487	-1.200699716	-0.5	-0.166666667


here is an example script for the provided test input, assuming that all 4 statistics (-v snipre,snipreb,alpha,dos) want to be resampled for AMPs and Antiviral RNAi (-c AMP,Antiviral RNAi) for 1000 bootstraps each (-b 1000). In this example, we'll only include background genes within 100kb of the focal gene (-w 100000), but won't care about size differences (-s not used)

    python resimpler.py -b 1000 -i test_input.txt -o test_output.txt -v snipre,snipreb,alpha,dos -c AMP,Antiviral RNAi

Gives the following output:
    group	class	rep	snipre	snipreb	alpha	dos
    mel	AMP	0	-0.2352687808888889	0.7204883728888889	-0.938230933	0.10225491233333332
    mel	AMP	1	-0.3494176825555555	0.979108613111111	-0.19754234700000012	0.08438851266666667
    mel	AMP	2	-0.13725086766666664	0.3217980618888889	-0.023144964555555654	0.059921065111111095
    mel	AMP	3	-0.2591489235555555	0.41454101455555553	-1.0316762664444445	0.02733241666666665
    mel	AMP	4	0.02467841322222221	0.4750443162222222	-0.2696062705555555	0.09886944844444442

More statistics can be assessed by increasing the list, or including a text file with columns that should be resampled over. Classes of genes can also be included as a text file.
