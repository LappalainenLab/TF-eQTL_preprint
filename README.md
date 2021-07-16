# TF-eQTL
Data tables and scripts for <i>Transcription factor regulators of eQTL effects across individuals and tissues</i>.

## Overview of repository
<ul>
  <li><a href="https://github.com/LappalainenLab/TF-eQTL/tree/master/data_tables">data_tables/</a></li>: Data tables organized by analysis
  <li><a href="https://github.com/LappalainenLab/TF-eQTL/tree/master/plots">plots/</a></li>: Plots from manuscript and scripts to generate plots
  <li><a href="https://github.com/LappalainenLab/TF-eQTL/tree/master/scripts">scripts/</a></li>: Scripts organized by analysis
</ul>

## Analyses included
<ul>
  <li><b>Variant sets</b></li>:
  <ul>
    <li>Overlap of >=5% MAF GTEx variants with transcription factor ChIPseq and motif information</li>
    <li>Overlap with Caviar fine-mapped credible sets of GTEx eQTLs</li>
    <li>Effect size (aFC) calculation for variant-gene pairs</li>
  </ul>
  <li><b>Within-tissue TF-eQTLs</b></li>:
  <ul>
    <li>TF-eQTL interactions discovered with tensorqtl</li>
  </ul>
  <li><b>Cross-tissue TF-eQTLs</b></li>:
  <ul>
    <li>Cross-tissue expression-based TF-eQTL correlations discovered with Spearman correlations</li>
    <li>Cross-tissue protein-based TF-eQTL correlations discovered with Spearman correlations</li>
  </ul>
  <li><b>Overlap</b></li>:
  <ul>
    <li>Permutations and calculations for ChIPseq and motif overlap of discovered TF-eQTLs</li>
    <li>Regulon gene overlap of discovered TF-eQTL genes</li>
    <li>Permutations and calculations for ADASTRA allele-specific TF binding overlap of discovered TF-eQTLs</li>
  </ul>
  <li><b>IRF1 knockdown</b></li>:
  <ul>
    <li>Allele specific expression in an IRF1 knockdown experiment</li>
  </ul>
  <li><b>Gene-by-environment interactions and GWAS</b></li>:
  <ul>
    <li>TF-eQTL overlap with gene-by-environment interactions</li>
    <li>TF-eQTL overlap with GWAS-eQTL colocalizations</li>
  </ul>
</ul>
