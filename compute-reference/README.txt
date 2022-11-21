Updated methods to calculate reference vectors. Below are the changes compared to the past method.
- 1) randomly sample half of the normal samples from current cohort, and use PC1 vectors of them only when calculating reference vector
  - this is to prevent information leak, so that held-out normal samples can be evaluated by the reference vector from which the information regarding held-out normal samples are excluded.
- 2) does not use minmax_scaling to the range of [0,1]
  - since distance is always >= 0, all samples which are represented as a point (score2, score4) will be located in the first quadrant. 
  - therefore, there is no need to scale the scores. 
