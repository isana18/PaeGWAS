# PaeGWAS
Source code for publication "The Pseudomonas aeruginosa accessory genome elements influence virulence towards Caenorhabditis elegans" by Alejandro Vasquez-Rifo, Isana Veksler-Lublinsky, Zhenyu Cheng, Frederick M. Ausubel and Victor Ambros.


The code provides a statistical test for association of genetic elements with virulence.

To run the code we generated a presence/absence matrix for 11,731 homologues gene clusters across 52 strains. We also supplied worm median survival (corresponds to virulence) for the 52 strains.

The Mann-Whitney (MW) ranking test and linear-regression (LR) analysis were applied to every gene to test the association of the presence/absence pattern with virulence.

We then performed a permutation test to control for multiple hypothesis testing. 10,000 permutations of the virulence values and their assignment to strains were generated and the MW and LR association tests were repeated for each permutation. Then, for each gene the number of times that it received a better p-value using the shuffled virulence data compared to the original one was recorded, separately for MW and LR. The above count was divided by 10,000 to obtain the permutation corrected p-value for the MW and LR tests.
