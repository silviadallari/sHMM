The sHMM repository provides the R implementation of a supervised Hidden Markov Model for the classification of tractor operational states based on compressed segment-level data with mixed-type covariates.

The implementation follows the methodology described in *A Scalable and Interpretable Supervised Hidden Markov Model for Tractor Operational State Prediction*

and includes the core estimation and prediction routines for the supervised model.

The model is formulated within a generative Hidden Markov framework in which the latent states follow a first-order homogeneous Markov chain and the observed covariates are modeled through state-dependent emission distributions. 

Parameter estimation is carried out in a fully supervised setting via direct maximization of the complete-data likelihood, yielding closed-form estimators for emission parameters and empirical estimators for the initial and transition probabilities. State prediction on unlabeled sequences is performed through forward–backward recursions implemented on the log scale for numerical stability.

The file sHMM.R contains the full implementation of the model, including estimation of emission parameters, transition probabilities, and posterior state decoding.
