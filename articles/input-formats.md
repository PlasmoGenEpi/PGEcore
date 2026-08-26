# Standard input and output formats

PGEcore tools read and write **TSV** tables with agreed column names so
outputs from one step can be inputs to another. Function and CLI
arguments use the same schema names as the sections below
(`allele_table`, `aa_calls`, `snp_calls`, `coi_table`, `loci_groups`,
…). Tool help pages describe *which* tables they need; this vignette is
the place for column layouts and examples.

Bundled examples live in `inst/extdata/` after you install the package
(or in that folder in the source tree).

``` r

library(PGEcore)
extdata <- system.file("extdata", package = "PGEcore")
list.files(extdata, pattern = "^example")[1:8]
#> [1] "example_aa_calls.tsv"            "example_allele_table.tsv"       
#> [3] "example_coi_table.tsv"           "example_coiaf_plmaf.tsv"        
#> [5] "example_collapsed_snp_calls.tsv" "example_ecoi.tsv"               
#> [7] "example_heterozygosity.tsv"      "example_loci_groups.tsv"
```

The snippets below use those examples only to show **column structure**.
In normal use, pass paths to **your** files on the CLI or in R.

## COI table (`coi_table`)

**Columns:** `specimen_name`, `coi`

**Example:** `example_coi_table.tsv`  
**Used by:** `count_samples_by_coi`, and as optional COI input to
several wrappers (`dcifer_*`, …).  
**Exception:** `FreqEstimationModel_wrapper` uses `--coi`, which accepts
either a path to this table or a single numeric average COI.

``` r

readr::read_tsv(
  file.path(extdata, "example_coi_table.tsv"),
  show_col_types = FALSE,
  n_max = 3
)
#> # A tibble: 3 × 2
#>   specimen_name                                                              coi
#>   <chr>                                                                    <dbl>
#> 1 PARAV3-ENV-MH04-7S1-7C1-1000-parasitedensity-sampleDB-4064911117_S43_L0…     3
#> 2 PARAV3-ENV-MH04-DS2-DC11-1000-parasitedensity-sampleDB-4064911661_S35_L…     2
#> 3 PARAV3-ENV-MH04-DS2-DC11-10000-parasitedensity-sampleDB-4064911565_S30_…     3
```

``` bash
Rscript exec/count_samples_by_coi \
  --coi_table my_coi_table.tsv \
  --output coi_distribution.tsv
```

## SNP calls (`snp_calls`)

**Minimum columns:** `specimen_name`, `snp_name`, `reads`, `seq_base`

Collapsed SNP tables may include extra columns (`target_name`, `chrom`,
`pos`, …). Tools that need only the minimum set ignore the rest.

**Example:** `example_collapsed_snp_calls.tsv`  
**Used by:** `coiaf_wrapper`, `THEREALMcCOIL_wrapper`,
`snp_calls_to_vcf`, SNP filtering and pileup tools, …

``` r

readr::read_tsv(
  file.path(extdata, "example_collapsed_snp_calls.tsv"),
  show_col_types = FALSE,
  n_max = 3
) |>
  dplyr::select(specimen_name, snp_name, reads, seq_base)
#> # A tibble: 3 × 4
#>   specimen_name snp_name           reads seq_base
#>   <chr>         <chr>              <dbl> <chr>   
#> 1 Laos2017-01   Pf3D7_01_v3-181574  4648 G       
#> 2 Laos2017-01   Pf3D7_01_v3-181658  4648 G       
#> 3 Laos2017-01   Pf3D7_01_v3-528906  5032 G
```

``` bash
Rscript exec/THEREALMcCOIL_wrapper \
  --snp_calls my_snp_calls.tsv \
  --slaf_output slaf.tsv \
  --coi_output coi.tsv
```

## Allele table (`allele_table`)

**Columns:** `specimen_name`, `target_name`, `seq`, `reads`

**Example:** `example_allele_table.tsv`  
**Used by:** `estimate_coi_naive`, `dcifer_slaf_wrapper`,
`moire_wrapper`, `malariaem_wrapper`, `allele_per_locus_summary`,
`per_locus_popgen_summary`, …

``` r

readr::read_tsv(
  file.path(extdata, "example_allele_table.tsv"),
  show_col_types = FALSE,
  n_max = 3
)
#> # A tibble: 3 × 4
#>   specimen_name                                          target_name seq   reads
#>   <chr>                                                  <chr>       <chr> <dbl>
#> 1 PARAV3-ENV-MH04-DS2-DC3-1000-parasitedensity-sampleDB… Pf3D7_01_v… GATA…    47
#> 2 PARAV3-ENV-MH04-DS2-DC3-1000-parasitedensity-sampleDB… Pf3D7_01_v… GATA…   822
#> 3 PARAV3-ENV-MH04-DS2-DC11-1000-parasitedensity-sampleD… Pf3D7_01_v… GATA…   510
```

``` bash
Rscript exec/dcifer_slaf_wrapper \
  --allele_table my_allele_table.tsv \
  --slaf_output slaf.tsv
```

## AA calls (`aa_calls`)

**Core columns:** `specimen_name`, `gene_id`, `aa_position`, `aa`,
`reads`  
Often also: `target_name`, `aa_locus`, `gene`, `ref_aa`

**Example:** `example_aa_calls.tsv`  
**Used by:** `estimate_allele_frequency_naive`,
`estimate_allele_prevalence_naive`, `filter_biallelic_calls`,
`multilocus_prevfreq_*`, `snpslice_wrapper`, …

``` r

readr::read_tsv(
  file.path(extdata, "example_aa_calls.tsv"),
  show_col_types = FALSE,
  n_max = 3
)
#> # A tibble: 3 × 9
#>   specimen_name target_name  reads gene    aa_locus   gene_id aa_position ref_aa
#>   <chr>         <chr>        <dbl> <chr>   <chr>      <chr>         <dbl> <chr> 
#> 1 specimen1     pfdhfr_1_150     5 dhfr-ts PF3D7_041… PF3D7_…          51 N     
#> 2 specimen1     pfdhfr_1_150     5 dhfr-ts PF3D7_041… PF3D7_…          59 C     
#> 3 specimen1     pfdhfr_1_150     5 dhfr-ts PF3D7_041… PF3D7_…         108 S     
#> # ℹ 1 more variable: aa <chr>
```

## Loci groups (`loci_groups`)

**Columns:** `group_id`, `gene_id`, `aa_position` (often `aa_locus` too)

Defines which loci are analysed together for multi-locus tools.

**Example:** `example_loci_groups.tsv`  
**Used by:** `multilocus_prevfreq_*`, `snpslice_wrapper`,
`FreqEstimationModel_wrapper`, `MultiLociBiallelicModel_wrapper`

``` r

readr::read_tsv(
  file.path(extdata, "example_loci_groups.tsv"),
  show_col_types = FALSE,
  n_max = 4
)
#> # A tibble: 4 × 4
#>   group_id      aa_locus            gene_id         aa_position
#>   <chr>         <chr>               <chr>                 <dbl>
#> 1 pfdhfr_pfdhps PF3D7_0417200.1:51  PF3D7_0417200.1          51
#> 2 pfdhfr_pfdhps PF3D7_0417200.1:59  PF3D7_0417200.1          59
#> 3 pfdhfr_pfdhps PF3D7_0417200.1:108 PF3D7_0417200.1         108
#> 4 pfdhfr_pfdhps PF3D7_0810800.1:437 PF3D7_0810800.1         437
```

## Multilocus allele frequency (MLAF)

**Columns:** `group_id`, `variant`, `freq`

`variant` is often a STAVE-style string (`gene_id:aa_position:aa`).

**Example:** `example_mlaf.tsv`  
**Used by:** `slaf_from_stave_mlaf` (and produced by several multi-locus
tools)

``` r

readr::read_tsv(
  file.path(extdata, "example_mlaf.tsv"),
  show_col_types = FALSE,
  n_max = 3
)
#> # A tibble: 3 × 3
#>   group_id      variant                                                     freq
#>   <chr>         <chr>                                                      <dbl>
#> 1 pfdhfr_pfdhps PF3D7_0417200.1:51_59_108:N_R_N;PF3D7_0810800.1:437_540:… 0.333 
#> 2 pfdhfr_pfdhps PF3D7_0417200.1:51_59_108:N_C_S;PF3D7_0810800.1:437_540:… 0.267 
#> 3 pfdhfr_pfdhps PF3D7_0417200.1:51_59_108:N_R_S;PF3D7_0810800.1:437_540:… 0.0667
```

## Population-level minor allele frequency (PLMAF)

**Columns:** `snp_name`, `seq_base`, `plmaf`

**Example:** `example_coiaf_plmaf.tsv`  
**Used by:** `coiaf_wrapper` (optional; otherwise PLMAF is estimated
from the SNP table)

``` r

readr::read_tsv(
  file.path(extdata, "example_coiaf_plmaf.tsv"),
  show_col_types = FALSE,
  n_max = 3
)
#> # A tibble: 3 × 3
#>   snp_name           seq_base  plmaf
#>   <chr>              <chr>     <dbl>
#> 1 Pf3D7_01_v3-181574 G        0.450 
#> 2 Pf3D7_01_v3-181658 G        0.369 
#> 3 Pf3D7_01_v3-528906 A        0.0779
```

## Microhaplotype SLAF

**Columns:** `target_name`, `seq`, `freq`, `sample_total`

**Example:** `example_mhaps_slaf.tsv`  
**Used by:** `slaf_from_mhaps_freqs`

``` r

readr::read_tsv(
  file.path(extdata, "example_mhaps_slaf.tsv"),
  show_col_types = FALSE,
  n_max = 3
)
#> # A tibble: 3 × 4
#>   target_name                 seq                              freq sample_total
#>   <chr>                       <chr>                           <dbl>        <dbl>
#> 1 Pf3D7_01_v3-0181544-0181729 TTTCATTATTGTTTTCATTCTTTTTTTAAC… 0.254           25
#> 2 Pf3D7_01_v3-0181544-0181729 TTTCATTATTGTTTTCATTCTTTTTTTAAC… 0.199           25
#> 3 Pf3D7_01_v3-0181544-0181729 TTTCATTATTGTTTTCATTCTTTTTTTAAC… 0.347           25
```

## Common outputs

| Output | Typical columns | Example producers |
|----|----|----|
| COI per specimen | `specimen_name`, `coi` | `estimate_coi_naive`, McCOIL, snpslice, … |
| COI distribution | `coi`, `n`, `proportion` | `count_samples_by_coi` |
| Allele frequency | `variant`, `freq` (sometimes counts) | naive AF, McCOIL SLAF, … |
| Prevalence | `variant`, `prev`, … | `estimate_allele_prevalence_naive` |

When you add a new tool, reuse these layouts whenever possible so
downstream steps do not need format converters.
