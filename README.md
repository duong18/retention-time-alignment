# retention-time-alignment
This repository code is from our research project of large scale retention time (RT) alignment

File iRT_alignment_4_from_reference_args.py:

- Installed Python packages: argparse, csv, sys, random, scipy, numpy, sklearn, statistics
- Running command format: python iRT_alignment_4_from_reference_args.py --evidence <evidence file with Maxquant output format> --reference <references set file> --output <output file exporting the alignment results>
Example: python iRT_alignment_4_from_reference_args.py --evidence ./data/2f0a28e8c8fe43ec8c50256e117f52a2_evidence.tsv --reference 29tissues_ref_set.txt --output ./data/2f0a28e8c8fe43ec8c50256e117f52a2_mean_absolute_errors.txt
Where the evidence file should be in the Maxquant software's output format. An example of the file can be found in the data folder of this repository. The references set file will contain reference precursors selected and their assigned RTs, each line will contain a pair of precursor and RT. The output file will contain the alignment results - MAE (Mean Absolute Error) values for the run pairs in the evidence file. For extra evaluation, besides the main MAEs in column "A and B after aligning", we also output the alignment results from the other approaches in the other columns. Each line in the output file will contain MAE values for each pair of runs.