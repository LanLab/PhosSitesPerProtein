# Count residues in peptides

Here, we present a script to count the number of S,T,Y phosphorylated residues above a given threshold from MaxQuant data.  



### Example

We provide an example input file Phospho(STY)Sites.txt to use our script with. To run the script on this data:

```sh
python3 phospho_maxquant_process.py ---maxQuantInFile Phospho(STY)Sites.txt

```

The output is written to by default to the file Output.txt (this can be changed using the option --outfile below).


The output file contains the S, T, Y counts in the last 6 columns. The first column shows the peptide identifier, and the second column shows the quantified phosphopeptide.

### Usage
```sh

$ python3 phospho_maxquant_process.py -h
usage: phospho_maxquant_process.py [-h] --maxQuantInFile MAXQUANTINFILE
                                   [--protColumn PROTCOLUMN]
                                   [--locProbColumn LOCPROBCOLUMN]
                                   [--locProbTh LOCPROBTH]
                                   [--styProbColumn STYPROBCOLUMN]
                                   [--styProbTh STYPROBTH] [--outfile OUTFILE]

Extract S, T and Y site counts (which are above thresholds) from MaxQuant
data.


optional arguments:
  -h, --help            show this help message and exit
  --maxQuantInFile MAXQUANTINFILE
                        An output file from MaxQuant, as input to this script.
  --protColumn PROTCOLUMN
                        Column containing the protein identifier. Default is
                        0.
  --locProbColumn LOCPROBCOLUMN
                        Column containing the localization probability.
                        Default is 5.
  --locProbTh LOCPROBTH
                        Threshold to filter localization probability (between
                        0 and 1). Default is 0.75.
  --styProbColumn STYPROBCOLUMN
                        Column containing the phospho-(STY)-probabilities.
                        Default is 29.
  --styProbTh STYPROBTH
                        Threshold to filter phospho-(STY)-probabilities
                        (between 0 and 1). Default is 0.75.
  --outfile OUTFILE     Output file to save the results to. Default is
                        Output.txt.
```
