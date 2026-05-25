```bash
# Run with both mono- and polyclonals
Rscript scripts/malecot_wrapper/malecot_wrapper.R \
    --allele_table data/example2_allele_table.tsv \
    --threads 5 \
    --model_results_output MALECOT_res.rds

# Monoclonals only
Rscript scripts/malecot_wrapper/malecot_wrapper.R \
    --allele_table data/example2_allele_table_monos_only.tsv \
    --use_provided_mean_COI \
    --COI_mean 1 \
    --COI_max 1 \
    --threads 5 \
    --model_results_output MALECOT_res.rds
```
