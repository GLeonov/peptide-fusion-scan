# User configured experiment name and paths to DB-matched and de novo files
EXPERIMENTS=('N1_Normal' \
             'N2_Normal')
DENOVO_FILES=('data/N1_Normal/N1_Normal_all_de_novo.csv' \
              'data/N2_Normal/N2_Normal_all_de_novo.csv')
DB_MATCHED_FILES=('data/N1_Normal/N1_Normal_DB-matched_peptides.csv' \
                  'data/N2_Normal/N2_Normal_DB-matched_peptides.csv')

# User configured values that will be used for each run
MIN_PEPTIDE_LENGTH=7
MAX_PEPTIDE_LENGTH=25
ALC_SCORE_CUTOFF=50
MAX_FUSION_DISTANCE=40
REFERENCE_PROTEOME='data/SPHu.fasta'
INCLUDE_PTMS='True'

# execute runs sequentially
for i in "${!EXPERIMENTS[@]}"; do
python3 -u main.py --experiment_name "${EXPERIMENTS[i]}"\
                   --all_denovo_file "${DENOVO_FILES[i]}"\
                   --db_matched_file "${DB_MATCHED_FILES[i]}"\
                   --reference_proteome_file $REFERENCE_PROTEOME\
                   --min_peptide_length $MIN_PEPTIDE_LENGTH\
                   --max_peptide_length $MAX_PEPTIDE_LENGTH\
                   --alc_score_cutoff $ALC_SCORE_CUTOFF\
                   --max_fusion_distance $MAX_FUSION_DISTANCE\
                   --include_PTMs $INCLUDE_PTMS
done
