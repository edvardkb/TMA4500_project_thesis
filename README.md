This repository contains the code used in my project thesis:

**Parameter estimation of time-irreversible stochastic volatility models with TMB**

## Notes / current status

- `checkConsistency` is primarily set up for `SV_SHASH/SV_Model_SHASH_inference` with  
  `normal_transformed_latent_variables = TRUE`.

- `SV_SHASH_BURNIN` includes extra debugging for the case  
  `normal_transformed_latent_variables = TRUE`, where it logs Hessian eigenvalues right before the
  run fails (this typically happens when `burnin_time >= 2`).

- These Hessian-debug outputs usually *do not* appear when `normal_transformed_latent_variables = FALSE`,
  because inference can often run through in that setting as long as the SHASH shape parameters are fixed.

- `fixed_SV_params` is a list (up to 6 parameters) specifying which parameters to hold fixed during
  inference. Use this to:
  - restrict to the Gaussian special case, or
  - run controlled experiments where only selected SHASH/SV parameters are estimated.
  - The parameter names are: mu = E(y_i), alpha, phi, sigma > 0, delta > 0, epsilon.
