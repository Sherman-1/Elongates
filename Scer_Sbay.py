import polars as pl

cov = 0.5
seuil = 4

elongates = pl.read_csv(f"output/{cov}/{cov}_elongates.csv", infer_schema_length =10000)

########
# Get the clusters that are elongated in Scer but not in Sbay
# Knowing well that such clusters include clusters without Sbay, or only with Scer or any other scenario 
########

Nter_scer_conditions = ((pl.col("species") == "Scer_NCBI") & (abs(pl.col("max_Nter") - pl.col("Nter_nb_aa")) < seuil) & (pl.col("max_Nter") >= 10))
Nter_sbay_conditions = ((pl.col("species") == "Sbay") & (pl.col("Nter_nb_aa") > seuil))

Nter_scer_clusters = set(elongates.filter(Nter_scer_conditions)["cluster_name"].to_list())
Nter_sbay_clusters = set(elongates.filter(Nter_sbay_conditions)["cluster_name"].to_list())

Nter_clusters = elongates.filter(pl.col("cluster_name").is_in(Nter_scer_clusters- Nter_sbay_clusters))



Cter_scer_conditions = ((pl.col("species") == "Scer_NCBI") & (abs(pl.col("max_Cter") - pl.col("Cter_nb_aa")) < seuil) & (pl.col("max_Cter") >= 10))
Cter_sbay_conditions = ((pl.col("species") == "Sbay") & (pl.col("Cter_nb_aa") > seuil))

Cter_scer_clusters = set(elongates.filter(Cter_scer_conditions)["cluster_name"].to_list())
Cter_sbay_clusters = set(elongates.filter(Cter_sbay_conditions)["cluster_name"].to_list())

Cter_clusters = elongates.filter(pl.col("cluster_name").is_in(Cter_scer_clusters- Cter_sbay_clusters))


########
# From those clusters, we want to keep only the ones that have Scer AND Sbay
########


# Nter
scer_df_nter = set(Nter_clusters.filter(pl.col("species") == "Scer_NCBI")["cluster_name"].to_list())
sbay_df_nter = set(Nter_clusters.filter(pl.col("species") == "Sbay")["cluster_name"].to_list())
common = scer_df_nter.intersection(sbay_df_nter)

full_filtered_Nter = Nter_clusters.filter(pl.col("cluster_name").is_in(common)).sort("cluster_name")

# Cter
scer_df_cter = set(Cter_clusters.filter(pl.col("species") == "Scer_NCBI")["cluster_name"].to_list())
sbay_df_cter = set(Cter_clusters.filter(pl.col("species") == "Sbay")["cluster_name"].to_list())
common = scer_df_cter.intersection(sbay_df_cter)

full_filtered_Cter = Cter_clusters.filter(pl.col("cluster_name").is_in(common)).sort("cluster_name")


# Question : are there some clusters that correspond to both Nter and Cter elongation criteria ?

# Answer : 

# full_filtered_Nter.filter(pl.col("cluster_name").is_in(full_filtered_Cter["cluster_name"].unique().to_list()))