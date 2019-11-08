# LDM-RF
A pipeline for generating pangenome sequence tags
Authors: Shang Gao, Jiri Stiller @ A&F, CSIRO, St Lucia, QLD, Australia

The pileline consists of four perl scirpts:

1. create_unique_records.pl : creates an integrated unique record (UR) hash and split UR hashes for parallel computing.

2. calculate_sum.pl :  This scirpt computes and descending sorts the statistic ‘SUM’ for each UR in split hash against whole SNP matrix.

3. creat_index.pl : generates the one-to-many links between unique record and GBS tags from a hash.

4. creat_ML_features.pl : generates features used for training machine learning models.

The python scirpt recored the codes perform data cleaning, model training and result visulizaiton.

The final GBS tags, including Non-PAV and PAV tags, for domesticated and wild panels were archived in the file: GBS_tags_collection_Gao_etal_2019.7z.
