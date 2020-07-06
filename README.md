# ntcirseq
In-house python script for detecting non-tandem repeats in cir-seq data

## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk/)

## Status
Active Development

## Introduction
Circular sequencing was initially developed by [Andino Group at UCSF](https://andino.ucsf.edu/) in 2014 for detecting low-frequency variants in the study of RNA virus by reducing sequencing errors [1]. The authors also released their data analysis software, namely CirSeq, which is open-access at https://andino.ucsf.edu/CirSeq.

However, when sequencing transcriptomes of advanced organisms, circular sequencing generates non-tandem repeats due to excessively long RNA fragments. Therefore, we present our solution, namely **ntcirseq**, to address this problem by rewriting the Consensus Generation Module of [CirSeq](https://andino.ucsf.edu/CirSeq) using a seed-extension strategy in consensus sequence generation. This code was in part adapted in our EmPC-seq study [2].

## Usage
1. Download the original CirSeq software from https://andino.ucsf.edu/CirSeq.
2. Replace `ConsensusGeneration.py` and `ConsensusModule.pyx` with our ntcirseq scripts in this repository.
3. Compile the codes by typing `python setup.py build_ext --inplace`.
4. Call the function using `python ConsensusGeneration.py ${WORKDIR} ${FASTQS}`.

## Reference
[1] Acevedo, A., Brodsky, L., & Andino, R. (2014). Mutational and fitness landscapes of an RNA virus revealed through population sequencing. Nature, 505, pp.686-690.

[2] Cheung, P.\*, Jiang, B.\*, Booth, G.T., Chong, T.H., Unartar, I.C., Wang, Y., Suarez, G.D., Wang, J.#, Lis, J.T.#, Huang, X.#. Identifying Transcription Error-Enriched Genomic Loci Using Nuclear Run-On Circular-Sequencing Coupled with Background Error Modeling. J. Mol. Biol., 432(13):3933-3949, 12 June 2020.

## Contact
For technical questions, please contact Biaobin via email: biaobinjiang@gmail.com
