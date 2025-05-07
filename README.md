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

Single-cell sequencing is a useful technique to study gene regulation of individual cells.  With single-cell RNA-seq (scRNA-seq) and single-cell ATAC-seq (scATAC-seq) datasets publicly available, we further explore alternative polyadenylation, one of the main mechanisms of gene regulation, in a breast cancer model.


<br>

## II. Methods

Datasets. <i>scRNA-seq</i>: data were obtained from NCBI GEO accession GSE176078 (Wu <i>et al.</i>, 2021).  Ten samples (ER+, primary breast tumor: CID3941, CID4040, CID4463, CID4535; HER2+, primary breast tumor: CID3838, CID3921, CID45171; triple-negative breast cancer (TNBC), primary breast tumor: CID3946, CID4465, CID44041) were selected in this study.  <i>scATAC-seq</i>: data (MCF-7 DMSO) were obtained from GEO accession GSE190162 (Bommi-Reddy <i>et al.</i>, 2022).

<br>

Data processing.  <i>scRNA-seq</i>: data in each sample were first filtered based on minimum number of cells required for each gene (removed if < 10), highly variable genes (kept only top 2,000 genes using Seurat v3 method), doublet detection using SOLO model (Bernstein <i>et al.</i>, 2020), minimum number of genes required for each cell (removed if < 200), and outliers (removed if gene count per cell exceeds 98 percentile).  Cells were also removed if mitochondrial and ribosomal gene percentages were greater than or equal to 20%.  Samples were then combined afterward.  Data were normalized and filtered additionally for highly variable genes prior to cell type clustering.


<br>

## III. Results

PASs in the terminal regions, i.e. the 3'UTR regions, were examined in this study.  There are 247,852 unique TR PASs, which are associated with 29,434 unique genes.  Correlations between PAS features were compared (Figure 1).  The correlations bewtween the three models were as follow: PolyAID vs PolyAStrength: 0.61; PolyAID vs PolyA_SVM: 0.55; and PolyAStrength vs PolyA_SVM: 0.52.  

There are many diseases associated with cleavage and polyadenlyation activites, including cancer, aging, and skin-related diseases.  For example, FIP1L1 gene, which plays a role in leukemia (Ali <i>et al.</i>, 2023), enhances usage of proximal PASs (global 3'UTR shortening), while knockdown of FIP1L1 expression leads to usage of distal PASs (global 3'UTR lengthening) (Davis <i>et al.</i>, 2022, Li <i>et al.</i>, 2015).  By profiling FIP1L1 based on expression level and modeling, PAS genomic locations chr4:+:53459611 (PolyAID: 0.9986, PolyAStrength: 0.5618, PolyA_SVM: 0.9911) and chr4:+:53459667 (PolyAID: 0.9768, PolyAStrength: -7.8462, PolyA_SVM: 0.9750) ranked at the top two (Figure 2).


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

Bernstein NJ, Fong NL, Lam I, Roy MA, Hendrickson DG, and Kelley DR.  2020.  Solo: Doublet Identification in Single-Cell RNA-Seq via Semi-Supervised Deep Learning.  Cell Syst, 11(1):95-101.e5. doi: 10.1016/j.cels.2020.05.010.

Bogard N, Linder J, Rosenberg AB, and Seelig G. 2019.  A Deep Neural Network for Predicting and Engineering Alternative Polyadenylation.  Cell, 178(1):91-106.e23.  doi: 10.1016/j.cell.2019.04.046.

Bommi-Reddy A, Park-Chouinard S, Mayhew DN, Terzo E, Hingway A, Steinbaugh MJ, Wilson JE, Sims RJ 3rd, and Conery AR.  2022.  CREBBP/EP300 acetyltransferase inhibition disrupts FOXA1-bound enhancers to inhibit the proliferation of ER+ breast cancer cells.  PLoS One, 30;17(3):e0262378.  doi: 10.1371/journal.pone.0262378.

Linder J, Koplik SE, Kundaje A, and Seelig G. 2022.  Deciphering the impact of genetic variation on human polyadenylation using APARENT2.  Genome Biol, 23(1):232.  doi: 10.1186/s13059-022-02799-4.

Stroup EK, and Ji Z. 2023. Deep learning of human polyadenylation sites at nucleotide resolution reveals molecular determinants of site usage and relevance in disease.  Nature Commun, 14(1):7378:1-17.  doi: 10.1038/s41467-023-43266-3.

Wu SZ, Al-Eryani G, Roden DL, Junankar S, Harvey K, Andersson A, Thennavan A, Wang C, Torpy JR, Bartonicek N, Wang T, Larsson L, Kaczorowski D, Weisenfeld NI, Uytingco CR, Chew JG, Bent ZW, Chan CL, Gnanasambandapillai V, Dutertre CA, Gluch L, Hui MN, Beith J, Parker A, Robbins E, Segara D, Cooper C, Mak C, Chan B, Warrier S, Ginhoux F, Millar E, Powell JE, Williams SR, Liu XS, O'Toole S, Lim E, Lundeberg J, Perou CM, and Swarbrick A.  2021.  A single-cell and spatially resolved atlas of human breast cancers.  Nat Genet, 53(9):1334-1347.  doi: 10.1038/s41588-021-00911-1.
