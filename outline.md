# Introduction
  - Earth's Exosphere and Hydrogen
      - definition, why do we care about exosphere (baliukin pg 2)
          - H exosphere indicates presence of water (important for planetary/exoplanetary exploration), Mars
          - unknown exosphere Lyman-α contributions pollute astronomical UV studies
      - (figure) sample radial profiles (radius vs density)
      - source of hydrogen (baliukin pg 1)
          - produced in lower atmosphere, dissociation of h2o & methane, upwards diffusion. escaping and nonescaping atoms
          - ballistic, satellite, escaping
      - responds to radiation pressure from sun (baliukin pg 2 → thomas & bohlin)
          - photon absorption provides antisolar momentum, reemission is anisotropic → net antisolar momentum
      - time varying response to storm
  - Prior Exospheric Measurements and Retrievals
      - #link(label("previous_mission"))[(figure) relative positions of missions]
      // - show results here - refer to historical techniques in @static_retrieval
      - OGO-5 - 1968, 24 Re apogee
          - early studies of Hydrogen geocorona, detectors saturated
          - detected H out to 7 Re
          - later proved existence of Lyman-α background due to interstellar hydrogen flowing through solar system (baliukin pg 2)
          - showed existence of geotail (thomas and bohlin)
      - Apollo 16 camera - 1872, 63 Re apogee
      - Galileo
          - took 1/2 day from perigee (1990-12-08) to get to Moon (62Re)
          - picture taken ~1990-12-13 → 5.5 days → 5.5x2x62 → 680 Re
      - hisaki 1D (exceed instrument) @hisaki
          - 1000km apogee
          - discovered exosphere extends out to 38 Re
      - SWAN/SOHO (1996)
          - discovered geocorona extends to 100Re, beyond moon @baliukin
      - PROCYON/LAICA - lyman alpha imaging camera (2015)
      - dynamics explorer 1 (rairden, 1986) - spherically symmetric fit
          - based on numerical solution to RT
      - IMAGE (2000)
          - geocoronal imager (GEO) instrument
          - Østgaard retrieval @ostgaard
      - TWINS
  - Carruthers Mission Overview
      - Why has exosphere not already been studied?
          - past measurements have been limited by either low cadence or poor vantage
      - Advantages of L1: good vantage, low fuel maintenance cost, solar power
      - Direct sensing with mass spec. impractical (too large)

      - Computational tomography in exospheric/atmospheric studies (reference to Gonzalo, Butala)
      - Suitability (or lack thereof of classical CT techniques in atmospheric tomography)
        - e.g. FBP requires regular circular path
  - Contributions of This Thesis
      - spherical raytracer that can be adapted to other missions studying planetary atmospheres
      - reconstruction algorithms built for the Carruthers mission
  - Thesis Organization

# Measurement Constraints
  - Chapter introduction section
      - what is lyman-α - spectral line of hydrogen
          - a useful indicator of hydrogen. sun is strong source of Lyman-α photons that are absorbed and reemitted through resonance scattering at the same wavelength
      - Carruthers carries two UV cameras capable of detecting Lyman-α
          - NFI for fast inner exosphere evolution
          - WFI for global estimation of slowly-changing outer exosphere
  - Carruthers Camera and Orbit Details
      - carruthers will be inserted in a halo orbit around langrange point L1 1.5e6 km from the Earth
          - from distant vantage, should be able to observe full 100Re radial extent of exosphere (soho/swan)
          - full orbit takes approximately 6 months (FIXME: SOHO claims 1 yr orbit, how?) (lara: there are other kinds of orbits, probably ok)
          - about ±7 degrees off sun-earth north/south (FIXME: terminology)
          - about ±28 degrees off sun-earth east/west (FIXME: terminology)
          - orbital ellipse limits: dX 71 Re, dY 221 Re, 54
              - FIXME: dY seems huge.  mistake?
      - Camera details
          - NFI: 3.6°, 30 min for fast inner exosphere evolution
          - WFI: 18°, 60 mins for global estimation of slowly-changing outer exosphere
          - spatial resolution
              - #link(label("los_evolution"))[(figure) LOS evolution over 15 days]
              - NFI: X°/pixel (Y Re when projected to Earth)
              - WFI: X°/pixel (Y Re when projected to Earth)
  - Emission Model
      - definition of optically thin and importance of this assumption on emission model computational complexity
      - definition of optical depth (zoennchen 2015)
      - g-factor
      - definition of albedo (gonzalo thesis appendix c)
      - other sources of signal
          - IPH, moon, stars
          - diagram showing hydrogen density in solar system?
  - Simulation, Calibration and Post-Processing
      - ability to simulate realistic images is critical to validate algorithm performance under noise conditions
      - brief description of MCP, KBr photocathode, etc. (use overview list from Instrument Model section)
      - IPH subtraction and instrument calibration
      - science pixel binning - computational simplificity
      - masking - 3Re, moon,
      - see @appendix_sim for simulation and calibration details
      - g-factor estimation (from Pratik) and conversion of calibrated radiance to column density


# Inverse Problem Formulation

- section introduction
    - purpose of tomography is retrieval of a volumetric structure which produced a set of measurements
    - this process is generally known as an inverse problem
- discretization
    - discretization necessary for numerical approximations to be made when analytic solution is infeasible
    - also rediscretize measurements to reduce computational burden (refer to section XXX. FIXME)
    - spherical grid definition in GSE coordinates
- inverse problem
    - under conditions described in @measurement_constraints (namely, single scattering in optically thin regime), tomographic inversion can be formulated as solution to the linear inverse problem $y = L x$ (ignoring noise)
    - direct inversion of the tomographic operator requires that the forward matrix L be non-singular in order for an inverse L^-1 to exist
    - measurement constraints are limited - less rows than free variables. ill posed
    - not injective (not one-to-one) (measurement might correspond to more than one potential solution)
    - not surjective (not onto) (a measurement may not necessarily even correspond to a feasible solution, e.g. noise)
    - still necessarily enough to get good retrieval under presence of noise
    - changes in density do not have significant impact on measurements.
    - inversely, small changes in measurements have large effects on solution during retrieval - ill conditioned @illposed
        - potentially problematic in the presence of noise
    - ill posedness and conditioning fixed through regularization, low rank models
    - concepts defined mathematically by hadamard @hadamard
    - generalized inverse (moore-penrose pseudoinverse), direct least squares minimization, SVD, RRPE
    - selecting hyperparameter values - trial and error, GCV, many others (gonz thesis pg 49)
    - solve with iterative method
- math notation (model, operator, coeffs, measurements, density), etc.


# Raytracing Optically Thin Gases <raytracer>
- section introduction
    - why FBP not appropriate.  SIRT, ART, TVMIN
    - fast projection are required
    - how iterative algs work
    - my contribution: raytracer which is differentiable, gpu-enabled, machine-learning enabled through pytorch, has visualizations
- gpu utilization - each pixel treated as separate task
    - out of core - definition.  my software is not capable of this
- autodifferentiability
    - important for some iterative algorithms
    - hand-coded gradient based algorithms slow to implement
    - table of capabilities of common tomography libraries
- definition of grid
    - assumptions this library makes about grid
- raytracer algorithm
    - summary of how raytracer is implemented
- api overview
    - this section summarizes how to use the library and gives examples of usage
    - SphericalGrid - how to use this class
    - ViewGeometry - how to use this class
    - combining view geometries to form a multi-vantage view geometry
    - Operator - instantiating an operator from grid and view geometry
    - reconstruction examples
- validation and tests
    - procedures used to verify correctness of code
- operator timing and benchmarking
    - how long operator takes to run and how much memory it uses

# Static Retrieval of Exosphere <static_retrieval>
- intro
    - numerous retrieval approaches to limtied exospheric data
    - some simple (1D estimate, few parameters), some more complex (HDOF)
    - summarizes point of this chapter
    - purpose of chapter is to introduce historical static retrieval approaches (generally increasing in complexity), then finally introducing a new static model contributed to this thesis

- 1D retrievals
    - early retrievals often relied on simple 1D models of exosphere
        - [x] either one-off measurements taken opportunistically (galileo flyby), limited viewing diversity
        - simple parametric formulations
        - in either case, a consequence of these simpler models is that are computationally easier to retrieve
        - avoid problems of underdetermined system having only a handful of free parameters
        - good measurement diversity from a single vantage → good inverse problem conditioning (esp under sph symm assumption)
    - SOHO/SWAN (1996) @baliukin
        - kinetic model of H atoms fit to observed data
        - simple model which ignores radiation pressure (chamberlain) and extension which includes radiation pressure
            - not quite 1D
        - requires knowledge of exobase (e.g. derived from NRLMSIS model)
        - applied onion-peeling technique in which outermost shells of model are fit to the data before proceeding inwards. high TP alt LOS pierce through few shells, low TP alt LOS pierce through more
    - ostgaard - IMAGE (2003) double exponential parametric form
        - sph symmetric
        - $n(r) = n_1 e^(-r/alpha_1) + n_2 e^(-r/alpha_2)$
        - this specific formulation may have been chosen because of ease of computing analytic coldens/radiance
        - density analytically converted to radiance given LOS and fit to data
    - PROCYON/LAICA (2015)
    - gonzalo MAP estimate
        - derived from earlier work from Butala
        - not enough viewing diversity, stereoscopic
    - Zoennchen 2024
        - functional form based on spherical harmonics
        - defined at single shell with exponential decay (and other terms)
    - spherical harmonic model
        - applicability of spherical harmonic bases to modelling exospheres
            - figure: show direct fits for different L and the max error on each
            - conclusion: spherical harmonic basis is a good low rank representation of exospheric models
        - splines enforce smoothness
            - agnostic shape
        - figure: spherical bases
        - single measurement A00 retrieval


# Static Retrieval Validation
- reconstruction requirements
    - spatial resolution requirements and reporting interval
- validation technique
    - [ ] how to choose sph harm order? - sph harm order sweep plots
        - direct fits to ground truth - determine with L works for retrieval - agnostic of viewing geometry
        - full retrievals (i.e. retrieval montecarlo simulation with known ground truth)
    - describe technique for verifying that each discretization is sufficient
        - density grid
        - LOS grid
        - science pixel binning
    - exosphere not completely understood and don't have access to ground truth data
    - must rely on models from physics simulations and prior retrievals from limited data
        - #link(label("datasets"))[(table) ground truth datasets]
    - motivate choice of discretization grid
        - nyquist argument - 2x highest frequency of continuous model (gonz thesis pg 52)
        - #link(label("stormbins"))[(table) storm time discretization]
        - #link(label("quietbins"))[(table) quiet time discretization]
    - definition of cartesian coordinate system (GSE) and spherical coordinate system, sph coordinate convention (e/a instead of phi/theta)
    - calibration bias - g factor, future work: IPH, radiation
    - #link(label("codeoverview"))[(figure) retrieval block diagram]

# Dynamic Retrieval of Exosphere
- work in progress

# Dynamic Retrieval Validation
- work in progress
  - reconstruction requirements
      - spatial resolution requirements and reporting interval - same as static
      - change detection requirement

# Conclusion
- work in progress

# Appendix - Instrument Simulation and Background Subtraction
- this section describes how noise is simulated for the instrument and how to remove the effects of noise
- instrument model
    - complete mathematical description of instrument (camera) model from a perspective of random processes
    - components of the UV camera
        - optical filter wheel
        - kbr photocathode
        - micro channel plate
        - phosphor and CCD
        - ADC
    - simulation of instrument processes and frame stacking
        - central limit theorem approximation
        - numerical validation


# Appendix - Discretization Considerations <discretization>

- 3 sources of discretization error (below)
    + science pixel binning too coarse for col density → binning has implicit assumption that does not hold
        - i.e. values of pixels being binned are all roughly equal
        - column densities are being undersampled at 3Re from large H-density gradients
    + grid is too coarse for data
    + same as #1, but for shadow region instead of 3Re region
- procedure for verifying that grids are sufficiently fine:
    + obtain col. dens. measurements with very high resolution grid
    + obtain col. dens. measurements with target grid
    + check that % err between #1 and #2 is below some threshold
- this is pretty ad-hoc.  how to choose fine grid? threshold?
    - can this be made more formal? ask Farzad
