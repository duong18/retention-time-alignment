# retention-time-alignment
This repository code is from our research project of large-scale retention time (RT) alignment. This contains two main components: (1) build a high-quality set of peptide references together with their consensus RTs (by graph_alignment_with_consensus_weight_ver2.py) and (2) perform alignments and calculate RT errors for pair-wise runs given a set of peptide references constructed (by iRT_alignment_4_from_reference_args.py)

File iRT_alignment_4_from_reference_args.py:

- Installed Python packages: argparse, csv, sys, random, scipy, numpy, sklearn, statistics
- Running command format:
```
python iRT_alignment_4_from_reference_args.py --evidence <evidence file with Maxquant output format> --reference <file containing peptide references> --output <output file exporting the alignment results>
Example: python iRT_alignment_4_from_reference_args.py --evidence ./data/2f0a28e8c8fe43ec8c50256e117f52a2_evidence.tsv --reference ./data/29tissues_ref_set.txt --output ./data/2f0a28e8c8fe43ec8c50256e117f52a2_mean_absolute_errors.txt
or: python iRT_alignment_4_from_reference_args.py --evidence ./data/fdc3b305405e48508fe1fbf176186ba3_sample_evidence.tsv --reference ./data/29tissues_ref_set.txt --output ./data/fdc3b305405e48508fe1fbf176186ba3_sample_mean_absolute_errors.txt
```

Where the evidence file contains the information of identified peptides and should be in the Maxquant software's [output format](http://www.coxdocs.org/doku.php?id=maxquant:table:evidencetable). An example of the files can be found in the data folder of this Github repository or you can download one from [this workflow](https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task=fdc3b305405e48508fe1fbf176186ba3) (by navigating to the button "Download results with mapped identifications" or "Download results directly from MaxQuant output").

The references set file contains reference precursors selected and their assigned RTs that would be served as anchor precursors for alignment (you can find the detail of the alignment approach in our paper). Each line in the file should contain a pair of precursor and RT.

The output file will report the alignment result for all the run pairs in the evidence file. Each line in the output file will contain MAE value before (in column `"A and B originally"`) and after (in column `"A and B after aligning"`) for each pair.