# Supplementary-material
Supplementary R-code to "Structure without semantics? A reanalysis of Berthet et al. (2025)":
1. R-code for the full permutation test of the MCA appearing in Berthet et al. (2025) (Full MCA permutation.R)
2. the manipulated R-code of Berthet et al. (2015) feeding it with random data. (Random.R)
3. the manipulated R-code of Berthet et al. (2015) feeding it with random data with retained NA:s. (Random with NAs retained.R)


Full MCA permutation test

This script performs a permutation test of call-type structure in the MCA space. The observed ratio of between-call-type inertia to total inertia is first calculated from a weighted MCA. Call-type labels are then randomly permuted across observations, a new MCA is recomputed for each permutation, and the inertia ratio is recalculated.

Repeating this procedure generates a null distribution corresponding to the hypothesis that call-type labels are unrelated to feature profiles. The Monte Carlo p-value is the proportion of permutations producing an inertia ratio at least as large as the observed value.

Unlike tests conducted within a fixed MCA embedding, this approach allows the entire MCA geometry to vary across permutations and therefore evaluates whether the observed clustering of call types exceeds chance expectations while respecting the dependence structure of the ordination.


Random-data validation

This script tests the inferential procedure of Berthet et al. (2025) using simulated data containing no association between call types and contextual features. For each iteration, a random binary feature matrix is generated and passed through the same MCA → barycenter distance → linear model → emmeans workflow used in the original study.

Because the data contain no true signal, any significant result is by definition a false positive. Repeating the procedure thousands of times provides an empirical estimate of the false-positive rate and evaluates whether the analytical pipeline can produce apparently meaningful semantic structure from random data alone.

An additional option preserves the original NA structure while randomizing all observed values, allowing assessment of how patterns of missingness influence the resulting MCA geometry.

Link to Berthet et al.(2005): https://www.science.org/doi/10.1126/science.adv1170

Corresponding author: wartelandreas@hotmail.com
