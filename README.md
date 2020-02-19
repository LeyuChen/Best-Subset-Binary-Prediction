# Best-Subset-Binary-Prediction
Matlab implementation of the best subset maximum score binary prediction method proposed by Chen and Lee (2018).
Description of this prediction method and its computation details can be found in the paper:

Chen, Le-Yu and Lee, Sokbae (2018), "Best Subset Binary Prediction". 

The paper has been published at Journal of Econometrics. See https://www.sciencedirect.com/science/article/pii/S0304407618300770. The latest working paper version of this work can be found in this repository. 

The matlab function max_score_constr_fn can be used to compute the the best subset maximum score prediction rule via the mixed integer optimization (MIO) approach. A warm-start strategy based on refining the input parameter space can be employed to improve the MIO computational performance. The matlab function warm_start_max_score implements this warm-start strategy. These codes require the matlab version of the Gurobi solver for solving the MIO problems. The Gurobi solver is freely available for academic purposes.
