# Community detection for interaction networks

## Data Folder

This folder contains both of the data examples used in the paper
- Senate dataset : Records the voting alignment of all 100 senators for each amendment considered by the 107th U.S. Senate
- Zachary Karate network : Records the number of social interactions of thirty-four members in a karate club that experienced a split among its two leaders.

## Functions Folder

This folder contains the functions necessary for computation of the 
joint and conditional distributions as well as the search algorithm
for finding the maximum a posteriori clustering given the observed
interaction network.
- Probability distribution functions: p(G,B) and p(G,B,X)
- Search algorithm : \arg \max_B p(B | G)

## Analysis Folder

This folder contains the data analysis for the two network examples
- Senate dataset : Estimation of the maximum a posterior clustering with AND without covariates
- Karate network : Estimation of the maximum a posterior clustering
