1 Biological backdrop: melanin drives perceived eye colour
Iris stromal melanocytes vary chiefly in total melanin load and in the eumelanin : pheomelanin ratio.

Dark irides contain ≈ 3–5 × more total melanin than light irides and a higher proportion of eumelanin (Wakamatsu et al., 1998; Imesch et al., 1997).

Histology confirms the number of melanocytes is constant; colour differences arise from pigment-granule density, not cell counts.

Therefore perceived brightness tracks a single latent quantity (melanin), justifying an ordered “light → dark” scale.

2 Large-scale GWAS evidence

| Study                                       | N (discovery + replication)       | Main findings
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
| **Simcoe et al., 2021** *Science Advances*  | 192 986 Europeans                 | 124 independent signals at 61 loci; > 50 new. Melanogenesis genes (OCA2, TYR, TYRP1, SLC45A2…) dominate; 53 % of colour variance explained. |
| **Liu et al., 2010** *PLOS Genetics*        | 5 951 Dutch (+ 3 543 replication) | First to use photo-derived hue/saturation PCs (continuous). Three new loci (LYST, DSCR9, TTC3) detected only with quantitative phenotyping. |
| **Sturm et al., 2008** *Am. J. Hum. Genet.* | ≈ 3 000 Europeans                 | Single intronic SNP *HERC2* rs12913832 silences OCA2, explaining most blue-vs-brown contrast.                                               |


3 Phenotype scales already in use
Martin–Schultz 20-step chart (light 1a → 16 black-brown).

Digital image PCs / CIELAB scores – nearly linear with chemical melanin concentration.

4-level survey bins (blue → green → hazel → brown) – used by UK Biobank & 23andMe; order is monotonic in melanin load.

4 GWAS modelling strategy
4.1 Treat 4-level code (0 = blue … 3 = brown) as quantitative
PhenotypeScore = β₀ + β₁·Genotype + covariates + ε (using plink --linear)

Pros: simplest, maximum power.

Assumption: equal spacing between categories. Check by plotting marginal means; if clearly non-linear, switch to ordinal.

4.2 Robustness: ordinal logistic
Command: plink2 --glm ordinal (proportional-odds model).

No equal-spacing assumption; retains ordering with one parameter per SNP.

Test proportional-odds via Brant test; relax for violating SNPs if needed.

4.3 Fallback: multinomial logistic
Only if ordering hypothesis fails badly. (nnet::multinom in R).

3 × parameters per SNP → lower power; otherwise results usually concur with ordinal/linear.

5 Recommended workflow
Primary linear additive model (PLINK 1.9) on 0–3 score.

Secondary ordinal model (PLINK 2.0) to verify top hits; compare Manhattan & QQ plots.

Sensitivity multinomial model if proportional-odds assumption is violated.

QC plots: Manhattan & QQ for each model; correlation of −log₁₀P values across models.

Effect-size sanity check: plot genotype vs. mean phenotype category to demonstrate monotonicity.

References
Simcoe M, Valdes A, Liu F, et al. (2021). Genome-wide association study in almost 195 000 individuals identifies 50 previously unidentified genetic loci for eye colour. Science Advances 7(11): eabd1239. https://doi.org/10.1126/sciadv.abd1239

Liu F, Wollstein A, Hysi PG, et al. (2010). Digital quantification of human eye colour highlights genetic association of three new loci. PLOS Genetics 6(5): e1000934. https://doi.org/10.1371/journal.pgen.1000934

Sturm RA, Duffy DL, Zhao ZZ, et al. (2008). A single SNP in an evolutionary conserved region within intron 86 of HERC2 determines human blue–brown eye colour. American Journal of Human Genetics 82(2): 424-431. https://doi.org/10.1016/j.ajhg.2007.11.005

Wakamatsu K, Hu DN, Ito S, McCormick SA. (1998). Characterization of melanins in human irides and cultured uveal melanocytes from eyes of different colours. Experimental Eye Research 67(3): 293-299. https://doi.org/10.1006/exer.1998.0538

Imesch PD, Wallow IHL, Albert DM. (1997). The colour of the human eye: a review of morphologic correlates and of some conditions that affect iridial pigmentation. Survey of Ophthalmology 41(Suppl 2): S117-S123. https://doi.org/10.1016/S0039-6257(97)80018-5

Martin–Schultz scale. (2025, March 20). In Wikipedia. Retrieved May 16 2025, from https://en.wikipedia.org/wiki/Martin%E2%80%93Schultz_scale