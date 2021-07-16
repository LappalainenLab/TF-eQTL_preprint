# TF-eQTL
Data tables and scripts for <i>Transcription factor regulators of eQTL effects across individuals and tissues</i>.

## Overview of repository
<ul>
  <li><a href="https://github.com/LappalainenLab/TF-eQTL/tree/master/data_tables">data_tables/</a>: Data tables organized by analysis</li> 
  <li><a href="https://github.com/LappalainenLab/TF-eQTL/tree/master/plots">plots/</a>: Plots from manuscript and scripts to generate plots</li>
  <li><a href="https://github.com/LappalainenLab/TF-eQTL/tree/master/scripts">scripts/</a>: Scripts organized by analysis</li>
</ul>

## Analyses included
<ul>
  <li><b>Variant sets</b>:</li>
  <ul>
    <li>Overlap of >=5% MAF GTEx variants with transcription factor ChIPseq and motif information</li>
    <li>Overlap with Caviar fine-mapped credible sets of GTEx eQTLs</li>
    <li>Effect size (aFC) calculation for variant-gene pairs</li>
  </ul>
  <li><b>Within-tissue TF-eQTLs</b>:</li>
  <ul>
    <li>TF-eQTL interactions discovered with tensorqtl</li>
  </ul>
  <li><b>Cross-tissue TF-eQTLs</b>:</li>
  <ul>
    <li>Cross-tissue expression-based TF-eQTL correlations discovered with Spearman correlations</li>
    <li>Cross-tissue protein-based TF-eQTL correlations discovered with Spearman correlations</li>
  </ul>
  <li><b>Overlap</b>:</li>
  <ul>
    <li>Permutations and calculations for ChIPseq and motif overlap of discovered TF-eQTLs</li>
    <li>Regulon gene overlap of discovered TF-eQTL genes</li>
    <li>Permutations and calculations for (<a href="https://adastra.autosome.ru/">ADASTRA</a>) allele-specific TF binding overlap of discovered TF-eQTLs</li>
  </ul>
  <li><b>IRF1 knockdown</b>:</li>
  <ul>
    <li>Allele specific expression in an IRF1 knockdown experiment (<a href="https://www.biorxiv.org/content/10.1101/2020.02.21.959734v1">Brandt et al., 2020</a>)</li>
  </ul>
  <li><b>Gene-by-environment interactions and GWAS</b>:</li>
  <ul>
    <li>TF-eQTL overlap with gene-by-environment interactions (<a href="10.7554/eLife.67077">Findley et al., 2021</a>)</li>
    <li>TF-eQTL overlap with GWAS-eQTL colocalizations (<a href="10.1186/s13059-020-02252-4">Barbeira et al, 2021</a>)</li>
  </ul>
</ul>
