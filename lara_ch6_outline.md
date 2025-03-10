# Document Title
1. Reconstruction requirements 
    - Contractual spatial resolution requirements and reporting interval (these are not precisely defined in some ways) 
    - Precise mathematical interpretation of requirements 
2. Implementation approach justification (separating major topics with subsubsections or at least bold font headers as I have done below will make this easier to reference later.  My order below is intentional -- flows from L1C image processing choices to reconstruction model choices.)
    - Background/overview (no need for header here)
        - (figure) retrieval block diagram. 
        -  In presenting this diagram, make very clear that the simulation block is only used for the validation effort and that the on-orbit retrieval block will take actual L1C images.
        - Note that exosphere not completely understood and don’t have access to ground truth data -- must rely on models from physics simulations and prior retrievals from limited data.  
        - Table 11
    - Temporal binning/image stacking (e.g., 1 hour vs 6 hours)
        - governs SNR
        - Remind about temporal resolution requirements here -- longer than 6 hours would require a moving window.
    - Temporal baseline of stacked image sequence (e.g, 1 month vs 14 days vs 1 hour).  
        - Tradeoff between viewing geometry diversity and effective averaging of temporal variability in the scene.  
        - If you haven't already, describe expectations for timescales of changes -- note here that we only have one full-year model of quiet-time seasonal variability (Pratik -- check to see if this one includes realistic variations in solar activity) and one full-month model of quiet-time season + solar activity (MATE MSIS).   Note that quiet-time variations are so slow as to not affect reasonable image stacking described above.  Storm-time variation precludes a long baseline, point to Chapter 8.
        -  (figures) Error introduced by static assumption on quiet-time data for various observation window durations.
    - Image binning and masking 
        - Have you already motivated curvilinear science pixel geometry?
        - Tradeoffs involving compute time for raytracing
        - Azimuthal resolution, shadow masking (example figure of high resolution albedo)
        - Radial resolution, nonlinear gradients within curvilinear pixel, inner exosphere masking (example figure of radiance profile in radial distance, error plots)
        - Table 13
    - Reconstruction density grid 
        - Refer to Appendix D
        - Definition of cartesian coordinate system (GSE) and spherical coordinate system, sph coordinate convention (e/a instead of phi/theta)
        - nyquist argument - 2x highest frequency of continuous model (gonz thesis pg 52) 
        - Log vs. linear distribution in radial distance, etc.
        - Table 14 (Table 12 belongs in Chaper 8, no?)
    - Spherical harmonic order 
        - Direct fits to ground truth - determine which L works for retrieval - agnostic of viewing geometry 
        - Retrieval of full set of higher order terms on each radial shell is a huge number of coefficients.  Without constraints on radial smoothness, yields "ringing" errors.  Motivates spline approach...
    - Spline control point selection
3. Retrieval Performance 
    - Noiseless vs. noisy.  Precision (in terms of Monte Carlo variances)
    - Unbiased vs biased calibration.  
        - Retrieval algorithm should be able to cope with expected biases on orbit
        - Bias in absolute intensity calibration, photon (IPH) background, and instrument background subtraction effect accuracy of L1C radiance used in the retrieval.  Describe how.  Note that absolute intensity calibration scales scene linearly, as does g-factor.
        - Future work?  Figures on sensitivity to calibration bias in these parameters
    - Sensitivity to specified g-factor
        - Assumed g-factor is expected to have accuracy of 3%, measured at 30 minute cadence independently in pipeline (is a calibration product)
        - g-factor sensitivity tests are equivalently evaluating sensitivity to bias in the absolute calibration.
        - Figures.
    - Conclusion.  Summarize perfomance margin on retrieval requirements (accuracy is the only one relevant here for now, no?).



