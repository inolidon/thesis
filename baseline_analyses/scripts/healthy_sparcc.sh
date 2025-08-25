# Step 1: SparCC Estimation of the species
# data input is the main species file in txt format (rows = species, columns = samples)

python Compute_SparCC.py \
-n Experiment_Healthy \
-di "/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/data/Healthy_sp_sparcc.txt" \
-xi 50 \
-ni 100 \
--save_cor="/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/h_cor_sparcc.csv"


# Step 2: Make Bootstraps with species

python MakeBootstraps.py "/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/data/Healthy_sp_sparcc.txt" \
-n 100 \
-t "h_permutation_#.csv" \
-p "/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/pvals/"

# Step 3: SparCC Estimation for all files from Step 2 Bootstrapping

for i in `seq 0 99`; do
    python Compute_SparCC.py \
    --name Experiment_PVals \
    -di "/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/pvals/h_permutation_$i.csv" \
    --save_cor "/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/pvals/h_perm_cor_$i.csv" \
    --verbose False >> "/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/Pval_healthy.log"
done

# Step 4: Estimation of P-tests on all SparCC files from sub samples

python PseudoPvals.py \
"/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/h_cor_sparcc.csv" \
"/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/pvals/h_perm_cor_#.csv" \
100 -o "/Users/inolishennon/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/python/SparCC/project/output/healthy/pvals/h_pvals_one_sided.csv" -t one_sided

