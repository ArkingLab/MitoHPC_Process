# MitoHPC_Process

Run MitoHPC_process script for summary files

## Filters applied
### Variant level filters applied  
 -- FILTER column: no "strict_strand|strand_bias|base_qual|map_qual|weak_evidence|slippage|position|Homopolymer|clustered|fragment|haplotype"  
 -- INFO column: no "INDEL|HP(homopolymer region)", no Read_depth < 300  
 -- no biallelic het variants (for each sample only one het at each POS)  

### Sample level filters generated (all samples are kept in the per_sample file, with different filter flags)  
 -- suspicious contamination from MitoHPC  
 -- Copy number <= 40  
 -- het_count > 5  

## 05/08/2025
    - Add mito_genome_reference file for all annotation
    - Update process script for the new reference file

## 06/26/2024
    - Update script for mMSS calculation; add new_mito_score_m file

## 02/21/2024
    - Add script for MitoHPC output processing; add new_mito_score file
