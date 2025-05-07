# Single-Cell Transcriptomic and Chromatin Accessibility Profiling Identify Heterogenous Alternative Polyadenylation Gene Expression in Cell Types of Breast Cancer Microenvironment


<br>

## Table of Contents

I. Introduction

II. Methods

III. Results

IV. Discussion

V. Technologies

VI. Abbreviation

VII. Acknowledgements

VIII. References


<br>

## I. Introduction

Single-cell sequencing is a useful technique to study gene regulation of individual cells.  With publicly available single-cell RNA-seq (scRNA-seq) and single-cell ATAC-seq (scATAC-seq) datasets, we further explore alternative polyadenylation, one of the main mechanisms of gene regulation, in a breast cancer model.


<br>

## II. Methods

Datasets. <i>scRNA-seq</i>: data were obtain from NCBI GEO acession GSE176078 (Wu <i>et al.</i>, 2021).  Ten samples (ER+, primary breast tumor: CID3941, CID4040, CID4463, CID4535; HER2+, primary breast tumor: CID3838, CID3921, CID45171; triple-negative breast cancer (TNBC), primary breast tumor: CID3946, CID4465, CID44041) were selected in this study.  <i>scATAC-seq</i>: data were from ().  Single sample () was used.

<br>

Data processing.


<br>

## III. Results

PASs in the terminal regions, i.e. the 3'UTR regions, were examined in this study.  There are 247,852 unique TR PASs, which are associated with 29,434 unique genes.  Correlations between PAS features were compared (Figure 1).  The correlations bewtween the three models were as follow: PolyAID vs PolyAStrength: 0.61; PolyAID vs PolyA_SVM: 0.55; and PolyAStrength vs PolyA_SVM: 0.52.  

There are many diseases associated with cleavage and polyadenlyation activites, including cancer, aging, and skin-related diseases.  For example, FIP1L1 gene, which plays a role in leukemia (Ali <i>et al.</i>, 2023), enhances usage of proximal PASs (global 3'UTR shortening), while knockdown of FIP1L1 expression leads to usage of distal PASs (global 3'UTR lengthening) (Davis <i>et al.</i>, 2022, Li <i>et al.</i>, 2015).  By profiling FIP1L1 based on expression level and modeling, PAS genomic locations chr4:+:53459611 (PolyAID: 0.9986, PolyAStrength: 0.5618, PolyA_SVM: 0.9911) and chr4:+:53459667 (PolyAID: 0.9768, PolyAStrength: -7.8462, PolyA_SVM: 0.9750) ranked at the top two (Figure 2).  The result from APARENT2 shows that FIP1L1 has a preference for proximal PASs (delta log-odds-narrow < 0) (Figure 5), which confirms with the previous findings.  Similarly, CSTF2 gene, which involved in cancer, and proliferating cells (Xia <i>et al.</i>, 2014), favors 3'UTR shortening and is in accordance with the APARENT2 result (preference for proximal PASs; genomic locations: chrX:+:100840916 (polyAID: 0.9984, polyAStrenght: -4.5279, PolyA_SVM: 0.6831); chrX:+:100841518 (polyAID: 0.9908, polyAStrenght: -6.7179, PolyA_SVM: 0.7773); delta log-odds-narrow < 0) (Figure 2, 5).  In addition to cancer, RRAS2 gene is also associated with aging.  It is reported that RRAS2 promotes 3'UTR lengthening in aging process, as in senescent cells (chen <i>et al.</i>, 2018).  Profiling of RRAS2 and comparison of proximal (chr11:-:14278745, polyAID: 0.7174, polyAStrenght: -6,5240, PolyA_SVM: 0.8469) and distal PASs (chr11:-:14277918, polyAID: 0.7682, polyAStrenght: -5.1554, PolyA_SVM: 0.8949) shows that RRAS2 has the same preference for distal PASs (delta log-odds-narrow > 0) (Figure 3, 5).  In skin-related diseases, such as systemic sclerosis, NUDT21 gene acts as an important regulator of alternative polyadenylation, and reduced exprssion of NUDT21 leads to predominant global 3'UTR shortening (Weng <i>et al.</i>, 2019).  Profiling of NUDT21 shows a preference for distal PASs (proximal PAS (chr16:-:56431674, polyAID: 0.9889, polyAStrenght: -5.3467, PolyA_SVM: 0.9730); distal PAS (chr16:-:56429141, polyAID: 0.9997, polyAStrenght: -1.3675, PolyA_SVM: 0.9925); delta log-odds-narrow > 0) (Figure 4, 5).


<br>

Figure.  .

<img src="figure/tr_corrmatrix.png" style="width: 50%; height: 50%;">


<br>

Figure.  Profiling of cancer-associated genes, FIP1L1 and CSTF2.

<img src="figure/genes_FIP1L1_CSTF2.png" style="width: 50%; height: 50%;">






## IV. Conclusion





<br>

## V. Technologies

Bioinformatics


<br>

## VI. Abbreviation

PAS: polyA site <br>




<br>

## VII. Acknowledgements




<br>

## VIII. References

Bogard N, Linder J, Rosenberg AB, and Seelig G. 2019.  A Deep Neural Network for Predicting and Engineering Alternative Polyadenylation.  Cell, 178(1):91-106.e23.  doi: 10.1016/j.cell.2019.04.046.

Linder J, Koplik SE, Kundaje A, and Seelig G. 2022.  Deciphering the impact of genetic variation on human polyadenylation using APARENT2.  Genome Biol, 23(1):232.  doi: 10.1186/s13059-022-02799-4.

Stroup EK, and Ji Z. 2023. Deep learning of human polyadenylation sites at nucleotide resolution reveals molecular determinants of site usage and relevance in disease.  Nature Commun, 14(1):7378:1-17.  doi: 10.1038/s41467-023-43266-3.


