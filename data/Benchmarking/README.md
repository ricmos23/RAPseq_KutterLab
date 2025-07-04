
To compare the motifs identified by our method with those obtained from other approaches (eCLIP-seq, RNA-Compete, and RBNS), we profiled RBP–RNA interactomes using HepG2 RNA. 
- eCLIP-seq data were retrieved from the ENCODE Project (https://www.encodeproject.org/) with the following accession numbers: HNRNPA1 (ENCFF797GSK), HNRNPC (ENCFF440ROZ), PTBP1 (ENCFF726SQU), RBFOX2 (ENCFF871NYM), and YBX3 (ENCFF185OEI). 
- RBNS data were also obtained from the ENCODE Project, with accession numbers: HNRNPC (ENCSR569UIU) and RBFOX2 (ENCSR441HLP) (Sequence, Structure, and Context Preferences of Human RNA Binding Proteins). 
- RNA-Compete data were retrieved from the publication, A Compendium of RNA-binding Motifs for Decoding Gene Regulation.

> We recomputed the position weight matrices (PWMs) for both RAP-seq and eCLIP-seq using the DREME tool (DREME: Motif Discovery in Transcription Factor ChIP-seq Data). As PWMs for RBNS were not available, motif files were directly downloaded from the ENCODE Project. RNA-Compete PWMs were obtained from the supplementary materials of the corresponding publication.

> https://github.com/ricmos23/RAPseq_KutterLab/blob/main/code/02_benchmarking_analysis_RAPSeq.Rmd