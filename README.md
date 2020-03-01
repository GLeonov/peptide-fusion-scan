# Peptide Fusion Scan App

Determines which peptides from MS/MS experimental data can be explained by peptide fusion.

Input:
* Peaks Software DB-matched peptides file
* Peaks Software de novo peptides file
* FASTA file of a reference proteome

Output:
* Pre-scan:
    * DB-matched peptides excluding razor peptides
    * DB-matched peptides excluding razor peptides, removing duplicate protein IDs
    * de novo peptide candidates for peptide fusion scanning
* Post-scan:
    * all scan results, and then further filtered for:
        * matched de novo peptides (forward-fused)
        * matched de novo peptides (reverse-fused)
        * matched de novo peptides (not fused, DB-matched)
        * unmatched de novo peptides (not fused, not DB-matched)

Example usage (single run):
```
python main.py --experiment_name N1_Normal --all_denovo_file data/N1_Normal/N1_Normal_all_de_novo.csv --db_matched_file data/N1_Normal/N1_Normal_DB-matched_peptides.csv --reference_proteome_file data/SPHu.fasta --min_peptide_length 7 --max_peptide_length 25 --alc_score_cutoff 50 --max_fusion_distance 40 --include_PTMs True
```

Example usage (sequential run, multiple experiments):
```
nohup ./start_peptide_fusion_scans.sh &
```

Author: Dr German Leonov, University of York, leonov89@gmail.com
