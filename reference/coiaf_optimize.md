# Call `coiaf::optimize_coi()` and resolve its monoclonal sentinel

When the frequency method finds no variant loci, coiaf returns `NaN`
carrying an `estimated_coi` attribute (1) rather than the estimate
itself. Reading that attribute keeps monoclonal specimens as COI 1
instead of dropping them. This wrapper therefore returns the estimate
itself otherwise this script would return NA for all monoclonal samples

## Usage

``` r
coiaf_optimize(...)
```
