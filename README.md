# Dynamic k-t MRI Image Reconstruction

Example script of some low-rank k-t MR image reconstruction with
sensitivity encoding and non-Cartesian sampling
Mark Chiew
mark.chiew@ndcn.ox.ac.uk

Compiled 24/07/19

These examples illustrate a few different methods for reconstructing
under-sampled k-t data, including:
  - un-regularised reconstruction (linear)
  - temporal finite difference L2 regularisation (linear)
    see Chiew & Miller, NeuroImage 2019 (https://doi.org/10.1016/j.neuroimage.2019.116165)
  - spatial finite difference L2 regularisation (linear)
  - spatio-temporal finite difference L2 regularisation (linear)
  - low-rank recon using soft thresholding (nonlinear)
  - low-rank recon using hard thresholding (nonlinear)
    see Chiew et al., MRM 2015 (https://doi.org/10.1002/mrm.25395) 
    and Chiew et al., MRM 2016 (https://doi.org/10.1002/mrm.26079)

NB: This example relies on MR encoding operators (for the sensitivity and non-Cartesian 
k-space encoding) that can be found at https://github.com/mchiew/recon-tools-matlab
