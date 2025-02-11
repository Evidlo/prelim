#import "template.typ": *
#import "@preview/glossarium:0.5.1": make-glossary, register-glossary, print-glossary, gls, glspl
#import "glossary.typ": glossary
#show: make-glossary
#register-glossary(glossary)

// custom functions for math/text
#import "functions.typ": *

// FIXME
// - make figure line width narrower*/
// - remove heading numbers when introducing autograd, gpu usage, discretization
// - fix tick sizes in plots

#let title = "Exospheric Hydrogen Density Retrieval Using Dynamic Optical Tomography"
#set document(title: title)
#show: ieee.with(
    title: title,
    abstract: [
    ],
    authors: (
        (
            name: "Evan Widloski",
            // department: "Electrical and Computer Engineering",
            // organization: "University of Illinois Urbana-Champaign",
            email: "evanw3@illinois.edu"
        ),
        (
            name: "Lara Waldrop",
            // department: "Electrical and Computer Engineering",
            // organization: "University of Illinois Urbana-Champaign",
            email: "lwaldrop@illinois.edu"
        ),
        (
            name: "Farzad Kamalabadi",
            // department: "Electrical and Computer Engineering",
            // organization: "University of Illinois Urbana-Champaign",
            email: "farzadk@illinois.edu"
        ),
    ),
    index-terms: ("raytracer", "projector", "spherical coordinates", "remote sensing", "tomography", "PyTorch", "GPU"),
    // bibliography-file: "refs.bib",
)
#set math.equation(numbering: "(1)", block: true)

// ----------- Document Styling ------------

#set heading(numbering: "1.")
#show heading: set block(above: 1.4em, below: 1em)
#set text(bottom-edge: "descender")
#set grid(row-gutter: 0.5em,)
// #show figure: set block(inset: (top: 0.5em, bottom: 1em))

#outline(indent: auto, depth: 2)

= Introduction
  #rt([
  - Earth's Exosphere and Hydrogen
      - definition, why do we care about exosphere (baliukin pg 2)
          - H exosphere indicates presence of water (important for planetary/exoplanetary exploration)
          - unknown exosphere Lyman-α contributions pollute astronomical UV studies
      - (figure) sample radial profiles (radius vs density)
      - source of hydrogen (baliukin pg 1)
          - produced in lower atmosphere, dissociation of h2o & methane, upwards diffusion. escaping and nonescaping atoms
          - ballistic, satellite, escaping
      - responds to radiation pressure from sun (baliukin pg 2 → thomas & bohlin)
          - photon absorption provides antisolar momentum, reemission is anisotropic → net antisolar momentum
      - time varying response to storm
  - Carruthers Mission Overview
      - Why has exosphere not already been studied?
          - past measurements have been limited by either low cadence or poor vantage
      - Advantages of L1: good vantage, low fuel maintenance cost, solar power
      - Direct sensing with mass spec. impractical (too large)
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
      #figure(
        grid(columns: 2, column-gutter: 1pt,
            subfigure(image("figures/procyon.png", height: 10em), "meas", "PROCYON/LAICA"),
            subfigure(image("figures/mooncarruthers.png", height: 10em), "meas", "Apollo 16 Carruthers camera")
        ),
          caption: "Previous measurements of exospheric Lyman-α"
      )

  - Contributions of This Thesis
      - spherical raytracer that can be adapted to other missions studying planetary atmospheres
      - reconstruction algorithms built for the Carruthers mission
  - Thesis Organization
  ])

  #figure(
      image("figures/previous_missions.svg", width: 100%),
      caption: "Past observations of exospheric Hydrogen at Lyman-α.\n* not representative of actual spacecraft location"
  ) <previous_mission>

  == Earth's Exosphere and Hydrogen

  == Carruthers Mission Overview

  == Prior Exospheric Measurements and Retrievals

  == Contributions of This Thesis

  == Thesis Organization


= Measurement Constraints <measurement_constraints>

  #rt([
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
      - #link(label("viewgeom1"))[(figure) viewing geometry, orbit around L1]
      - #link(label("viewgeom3"))[(figure) camera details, FOV, int. time, etc.]
      - #link(label("viewgeom2"))[(figure) static 3D picture of multiple view geometries captured for 1month baseline]
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
      - #link(label("scatteringphase"))[(figure) scattering phase function]
      - #link(label("emissionmod"))[(figure) column line integral]
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
  ])

  #image("figures/viewgeom1.svg", height: 5em) <viewgeom1>
  #image("figures/viewgeom2.gif", height:5em) <viewgeom2>
  #image("figures/viewgeom3.png", height:5em) <viewgeom3>
  #image("figures/los_evolution_15d.png", height:5em) <los_evolution>
  #image("figures/emissionmod.png", height:5em) <emissionmod>
  #image("figures/scatteringphase.png", height: 5em) <scatteringphase>

  == Carruthers Cameras and Orbit Details

  == Emission Model

  == Simulation, Calibration and Post-Processing

= Inverse Problem Formulation <inverse_problem>

  #rt([
      WIP outline
      - [x] purpose of tomography is retrieval of a volumetric structure which produced a set of measurements
          - this process is generally known as an inverse problem
      - [x] discretization
          - discretization necessary for numerical approximations to be made when analytic solution is infeasible
          - also rediscretize measurements to reduce computational burden (refer to section XXX. FIXME)
          - [x] spherical grid definition in GSE coordinates
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
  ])


  #figure(
      pad(x: -50%,
          table(
              table.header([Symbol], [Meaning], [Description], [Shape]),
              align: horizon,
              columns: (5em, 10em, 14em, 12em),
              $f$, [Forward operator], [Mapping from 3D H-density to column density], [$f: bb(R)^3 → bb(R)^3$ (static) \ $f: bb(R)^4→bb(R)^3$ (dynamic)],
              $m$, [Model], [Mapping from model \ parameters to 3D H-density], [$m: bb(R)^* → bb(R)^3$ (static) \ $m: bb(R)^* → bb(R)^4$ (dynamic)],
              $bold(c)$, [Model \ params./coeffs.], "Free model variables,\n usually low-dimensional.", [$bold(c) ∈ bb(R)^*$ \ \* model dependent],
              $bold(x)$, [H density], [Spatial distribution of \ exospheric Hydrogen], [$bold(x) ∈ bb(R)^3$ (static) \ $bold(x) ∈ bb(R)^4$ (dynamic)],
              $bold(y)$, [Measurements.], "Column densities measured\n by instrument", [$bold(y) ∈ bb(R)^3$],
          ),
      ),
      caption: [Symbols],
      // FIXME: move to glossary section?
  )

  An inverse problem is the procedure of determining the causative factors of a set of measurements derived from some observation process.  In exospheric tomography, the factor driving the intensity of column density measurements is the distribution of hydrogen in the regions being observed.  In general, a direct analytic solution to tomography or other inverse problems is not feasible, so numerical approximations and discretization become necessary.  In this chapter, I lay out key concepts of linear inverse problems, detail discretization for approaching tomographic inversion numerically, and introduce notation which will be used later in the manuscript to describe tomographic retrieval algorithms.

  == Discretization

    #rt("WIP")

  In general, direct analytic solutions to tomographic estimation problems are infeasible, necessitating a numerical approach.  Instead, numerical methods usually divide the solution space into a finite grid of $N$ non-overlapping regions called #glspl("voxel") where density is assumed to be constant.  There is a large variety of grids to choose from, including regular grids where voxels are uniformly spaced in some domain (e.g. spherical, cylindrical, cartesian) and non-regular grids which may utilize hierarchical structures (e.g. octree) or tetrahedral meshes.  Some of these discretizations schemes have been designed to adaptively update voxel boundaries during retrieval to better fit the object being retrieved @adaptivemesh1.

    - figure: types of grids. cartesian, spherical, octree, etc

  Choosing a discretization grid which is well-suited to the data is critical, as an improper grid can create aliasing and other artifacts which causes the numerical result to deviate from the functions they approximate.  Since the exosphere is well-understood to smoothly vary with density gradients increasing at lower altitudes, a regular spherical grid with logarithmically spaced radial bins is appropriate.  The nature of regular grids allows for a simple multidimensional array as the underlying data structure, and the property of shared boundaries between voxels simplifies tomography calculations, as demonstrated in @raytracer.

    This manuscript uses a regular spherical grid with zero-elevation pole aligned to ecliptic north (+Z Cartesian GSE axis) and azimuth branch point (±180°) pointed away from the sun (-X Cartesian GSE axis), as shown in @grid_details.

    #figure(
        box(width:100pt, height:100pt, stroke:1pt),
        caption: [Spherical grid definition and a spherical voxel\ #rt([FIXME: WIP])]
    ) <grid_details>

    This manuscript uses convention $r$, $e$, $a$ when referring to radial, elevational, and azimuthal dimensions to avoid ambiguity with astrophysical versus mathematical conventions for spherical coordinates as in @sph_convention.


    #figure(
        table(
            columns: (100pt, auto, auto, auto),
            align: center + horizon,
            table.header([Dimension], [This\ Manuscript], [Physics\ Convention], [Mathematical\ Convention]),
            [Radial], $r$, $r$, $r$,
            [Elevational], $e$, $theta$, $phi$,
            [Azimuthal], $a$, $phi$, $theta$
        ),
        caption: [Spherical coordinate conventions]
    ) <sph_convention>

  == Inverse Problem

    #rt("WIP")

  Under the conditions described in @measurement_constraints (single scattering in optically-thin exosphere) and ignoring noise, tomographic inversion can be formulated as the solution to the linear inverse problem

  #math.equation(
      $bold(y) = L bold(x)$
  )

    Direct inversion of the tomographic operator requires that the matrix $L$ be non-singular in order for an inverse $L^(-1)$ to exist.  However, tomography problems generally have fewer measurement constraints (#gls("LOS")) than free variables (#glspl("voxel")), known as an _underdetermined system_, making inversion of $L$ impossible as there are infinitely many solutions for $x$.
    A _generalized inverse_, such as the Moore-Penrose pseudoinverse, can be constructed from $L$ which selects the solution with the smallest norm or that best fits the measurements.  Another approach is to assume that $x$ has some low degree-of-freedom representation on a subspace (linear) or manifold (non-linear) on the space of solutions.  Such a mapping $m$ is referred to as a #gls("model") in this manuscript and is represented as

    #math.equation(
        $bold(x) = m(bold(c))$
    )

    where $bold(c)$ are the low degree-of-freedom coefficients.  If the model is sufficiently low-dimensional, then uniqueness is guaranteed if a solution exists.

    #rt([FIXME: awkward.  terminology: "instability"])

    Even when a unique solution exists, this solution may be sensitive to small perturbations in measurements $y$, which is problematic in the presence of noise.  This is especially true of computational imaging, where smoothing effects of integration dampen hight frequency information about $bold(x)$ in the measurements $bold(y)$.

    - define regularization

    The issues of solution existence and sensitivity in inverse problems were formally defined by Hadamard in the early 20th century and are known as _ill-posed_ and _ill-conditioned_ inverse problems @hadamard @illposed.



= Raytracing Optically Thin Gases <raytracer>

  // Tomography is a method for determining the internal structure of objects from a set of measurements which penetrate into the object being measured.
Fast tomographic reconstruction algorithms that implement explicit inversion formulas typically work only for specific view geometries (such as circular or helical path) and are referred to as _filtered back projection_ (FBP) algorithms @fbp.
  However, spaceborne sensors have view geometries that are determined by orbital parameters of a spacecraft.
  For these situations requiring more flexible view geometries where an exact inverse solution is not available, _iterative reconstruction_ (IR) algorithms prevail at the expense of increased computational complexity.  Examples include SIRT @sirt, TV-MIN @tvmin, ART @art, CGLS @cgls, Plug-and-play @plugandplay and many others.
  These algorithms depend on the generation of synthetic projections of a candidate object using an operator (sometimes called a _raytracer_) which simulates waves traveling through the object medium.
  They produce a reconstruction by repeatedly tweaking the candidate object to minimize discrepancy between synthetic and actual projections, and they stand to benefit the most from a fast operator implementation.

  In this section, we describe the design and implementation of a tomographic projector which is designed for use in iterative reconstruction algorithms on spherical grids.  The algorithm is automatically differentiable, GPU-enabled, and easily integrable into machine learning methods through PyTorch.
  Additionally, the implementation provides a suite of visualization functions described in @api_overview which can be used for visual verification of view geometry and discretization or for generating publishable figures.
  The raytracer has been released as a open-source Python package *TomoSphero*.

=== GPU Utilization

  Most tomographic operators treat each pixel on the detector as an independent computational task, which has led to the development of tomography libraries that are capable of simultaneously utilizing multiple cores on a CPU or hardware accelerator.  TomoSphero is parallelized and GPU-enabled, and its speed has been benchmarked in @benchmarking.

  In cases where a simultaneous computation for every pixel of every measurement would consume more memory than is available, some algorithms operate _out-of-core_, where they parallelize as many tasks as will fit into available memory, then serially queue the remaining tasks for processing after current tasks are complete.  TomoSphero is not capable of out-of-core operation, so a characterization of its memory usage is also available in @benchmarking.

=== Object Discretization

  #rt([FIXME: make this mesh with prevoius section about discretization, or delete])

  Another consideration in tomographic reconstruction is the choice of grid type for discretization of the reconstructed object.  Most publications consider a regular rectilinear grid, which is a reasonable choice when the underlying structure of the object is completely unknown or the scale of features is uniform throughout the object.  However, in cases where some prior information is known about the location of high-detail regions within the object, or when symmetries in the view geometry exist @thibaudeau, a good choice of coordinate system can sample the object more efficiently for decreased computational requirements and lower quantization errors.
  The primary focus of TomoSphero is in the domain of atmospheric tomography, where regular spherical grids are well-suited for modeling solar and planetary atmospheres which exhibit spherical symmetries @solartomography1 @solartomography2.
=== Autodifferentiability

  Computing gradients is an important part of many reconstruction algorithms.
  Many iterative reconstruction algorithms rely on gradient-based optimization to solve for an object whose structure corresponds to measurement data.  Hand-coded gradients give an exact solution but can be time-consuming and error-prone to implement correctly, especially in a language designed for a hardware accelerator or if the object lies in a transform domain (e.g. wavelets, truncated spherical harmonics).  Finite differentiation and symbolic differentiation exist to solve for gradients for any generic problem, but can be extremely slow for high-dimensional problems or can yield complicated expressions for the symbolic derivative (a.k.a expression swell).  Automatic differentiation (_autograd_) is a class of techniques which convert an arbitrary expression into a computational graph of simpler functions, then compute the overall derivative by applying chain rule at each node.  Modern machine learning libraries such as PyTorch @pytorch and Jax @jax provide such capabilities for building this computational graph.  TomoSphero is implemented on top of PyTorch and inherits its autograd capabilities as a result.

  A comparison of TomoSphero's capabilities against other popular libraries is shown in @other_libraries.

  #figure(
      table(
          columns: (130pt, 65pt, auto, auto, auto, auto, auto),
          inset: 4pt,
          align: horizon,

          // vert([Number of rays]),
          // table.cell(rowspan:2)[Num. Rays &\ Volume Shape],
          table.header(
              // [Name], [Volume Grid], [GPU\ Support], [Autograd], [Recon.\ Algs.], [Vis.], [Out-of-core]
              [Library],
              [Grid Type],
              // rotate([Grid Type], -90deg, reflow: true),
              rotate([GPU Support], -90deg, reflow: true),
              rotate([Autograd], -90deg, reflow: true),
              rotate([Reconstruction\ Algs. Provided], -90deg, reflow: true),
              rotate([Visualization], -90deg, reflow: true),
              rotate([Out-of-core], -90deg, reflow: true),
          ),
          table.hline(stroke: 2pt),
          [TIGRE @tigre], [Cartesian], [Yes], [No], [Yes], [No], [Yes],
          [LEAP @leap], [Cartesian,\ analytic], [Yes], [Yes], [Yes], [No], [Yes],
          [ASTRA @astra1 @astra2 @astra3], [Cartesian], [Yes], [Yes], [Yes], [No], [Yes],
          [mbirjax @mbirjax], [Cartesian], [Yes], [Yes], [Yes], [No], [No],
          [ToMoBAR @tomobar], [Cartesian], [Yes], [No], [Yes], [No], [Yes],
          [CIL @cil], [Cartesian], [Yes], [No], [Yes], [Yes], [Yes],
          [Tomosipo @tomosipo], [Cartesian], [Yes], [Yes], [Yes], [Yes], [Yes],
          // FIXME - exospy not a generic raytracer
          // [EXOSPy @exospy], [Spherical], [No], [No], [No], [No], [No],
          [TomoSphero], [Spherical], [Yes], [Yes], [No], [Yes], [No],
      ),
      caption: [Non-exhaustive overview of other tomography libraries]
  ) <other_libraries>

== Grid Definition

    The grid determines the location of voxels in 3D space.  In spherical coordinates the grid lines, or #emph[grid boundaries] are spheres, cones, or planes which correspond to the radial, elevational, and azimuthal dimensions as in @boundaries.  We refer to the space between two adjacent grid boundaries as a #emph[grid region].  For a spherical grid of shape ($N_r$, $N_e$, $N_a$), there should be $N_r+1$, $N_e+1$, and $N_a+1$ boundaries of each type, respectively.

    The user-defined locations of the grid boundaries need not be uniform nor must they cover the extent of the entire sphere.  This can be useful when spatially varying resolution is needed or only a small wedge needs to be raytraced, as in @wedge. @gridexample shows a 2D representation of a grid of shape (1, 2, 7) which does not cover the entire sphere and will be used as an example throughout the next section.


  #figure(
      image("figures/tomosphero/grid.svg", height: 10em),
      caption: [Radial, elevational, azimuthal grid boundaries for a grid with shape (2, 4, 8)]
      // FIXME: clip azimuthal planes at sphere, make +Z arrow black
  ) <boundaries>


  #figure(
      image("figures/tomosphero/alg_wedge.svg", width: 30%),
      caption: [A spherical grid with shape (2, 2, 2) defined over a small wedge]
  ) <wedge>

  #figure(
      image("figures/tomosphero/alg_grid.svg", height: 10em),
      caption: [Sample grid of shape (1, 2, 7) broken into radial, elevational, and azimuthal boundaries and shown in 2D.  An index is assigned to each boundary (in black) and each region (in gray).  The shaded area is within the grid and white area outside.]
  ) <gridexample>

== Raytracer Algorithm <alg_outline>

  The steps for this raytracer can be broken into 4 parts, which are described in detail in this section:

  + For a given ray, compute intersection points with all boundaries and their distances from ray starting point
  + Determine crossing direction for all intersection points and compute region index
  + Sort intersection points and regions by distance from ray starting point and find lengths between points
  + Raytrace ray by computing inner product of ray lengths with values from voxels along ray

  Since the view geometries in the tomography problem concerning this thesis are fixed and known, steps 1-3 can be precomputed and stored before reconstruction begins.  This saves computation time at the expense of increased memory usage.


=== Step 1 - Ray-Boundary Intersection Points

  The first step is to compute the intersection points of all rays with all boundaries and their distances from ray starting points.  An intersection point array and distance array is computed for each type of boundary separately.  We consider only a single ray in this section without loss of generality.

  It is desirable for the arrays allocated for ray-boundary intersection point coordinates and distances to be equal size for all rays to facilitate array programming and parallelization.  A ray may intersect a sphere or a cone 2 times and a plane 1 time.  Therefore, we assume an upper bound where each ray has exactly $2(N_r + 1)$, $2(N_e + 1)$, and $N_a + 1$ intersections with each type of boundary.

  For the cases where a boundary has no intersection points with a ray, coordinates are still computed and stored for these points, but the associated value in the distance array $t$ is set to #emph("inf"), as illustrated in @alg_intersection.

  #figure(
      image("figures/tomosphero/alg_intersection.svg", height: 10em),
      caption: [Intersection points along a #gls("LOS") labelled with index of associated boundary.  Points along bottom of figure are associated with non-intersecting boundaries.]
  ) <alg_intersection>


  The expressions for determining intersection points  are given in @int_raycone @int_raycone2 @int_raysphere
  .
=== Step 2 - Boundary Crossing Direction <section_crossdir>

  Since the overall goal is to determine the indices of particular voxels the ray is crossing as it passes through the grid, we must convert each boundary index computed in the previous section to a region index.  To do this, we must determine for each intersection whether the ray crosses the boundary in a "positive" direction (increasing boundary index) or a "negative" direction (decreasing boundary index) as shown in @alg_crossdir0.
  The algorithm for determining crossing direction is specific to each boundary type and depends only on the crossing point coordinate and ray direction.  It is described in @appendix_crossdir in the appendix.
  Finally, boundary indices are converted to region indices by decrementing indices of crossing points with a negative direction, as shown in @alg_crossdir and @alg_crossdir2.

  #figure(
      image("figures/tomosphero/alg_crossdir0.svg", height: 10em),
      caption: [Positive and negative boundary crossing directions]
  ) <alg_crossdir0>
  #figure(
      image("figures/tomosphero/alg_crossdir.svg", height: 10em),
      caption: [Intersection points along a LOS annotated with crossing direction (↑:postive, ↓:negative)]
  ) <alg_crossdir>
  #figure(
      image("figures/tomosphero/alg_crossdir2.svg", height: 10em),
      caption: [Region indices computed from crossing direction]
  ) <alg_crossdir2>

=== Step 3 - Sorting Intersection Points
    In this step we collect the region indices for all intersection points into a single list and sort the points by their distance $t$ from the ray start location, as illustrated in @alg_sort.  Next we compute the difference between adjacent points $Delta t$ which represents the intersection length of the ray with a particular voxel.  Note that this results in $Delta t$ equal to NaN or #sym.infinity for intersection points with $t=infinity$.
  #figure(
      image("figures/tomosphero/alg_sort.svg"),
      caption: [Intersection points along a LOS for all boundary types collected into a list and sorted by distance $t$ from ray start position, marked with x]
  ) <alg_sort>

    Next, we convert our list of region indices into a list of voxel indices by computing the voxel index of the ray starting location.  For each boundary crossed by the ray, the index of the next voxel changes in only a single dimension.

    This is illustrated in @alg_table (a), where the sorted region indices are inserted into an empty table starting with the ray starting point, then in @alg_table (b) missing indices are filled in from the previous row to form complete 3D voxel indices.

    We also account here for invalid region indices generated in step 2 where the ray passes outside the area defined by the grid or when boundaries have no ray intersection.  We achieve this by setting $Delta t=0$ for any rows where the region index is invalid (less than 0, or greater than or equal to the corresponding $N_r$, $N_e$, $N_a$), or where $Delta t$ is #sym.infinity or NaN.  Setting $Delta t=0$ negates any contribution of this row to the raytracing step.  This is shown in @alg_table (c).

  #figure(
  grid(
      columns: 3, align: bottom, gutter: 12pt,
      subfigure(grid(
              columns: 3, align: bottom,
              table(
                  // columns: (auto, auto, auto, auto),
                  columns: 4,
                  stroke: none,
                  // inset: 10pt,
                  align: horizon,
                  inset: 3.5pt,

                  // vert([Number of rays]),
                  // table.cell(rowspan:2)[Num. Rays &\ Grid Shape],
                  table.header("r", "e", "a", $Delta t$),
                  table.hline(),
                  table.vline(x: 3),
                  [1], [0], [0], [.3],
                  [0], [ ], [ ], [.8],
                  [ ], [-1], [ ], [.7],
                  [ ], [ ], [1], [.5],
                  [ ], [ ], [2], [.5],
                  [ ], [ ], [3], [.2],
                  [ ], [0], [ ], [.7],
                  [1], [ ], [ ], [#sym.infinity],
                  [0], [ ], [ ], [NaN],
                  [0], [ ], [ ], [NaN],
                  [], text(size:2em)[#sym.dots.h], [], [],
                  [ ], [ ], [7], [NaN],
              ),
          ),
          "sorting",
          "Sorted"
      ),
      subfigure(grid(
              columns: 3, align: bottom,
              table(
                  // columns: (auto, auto, auto, auto),
                  columns: 4,
                  stroke: none,
                  // inset: 10pt,
                  align: horizon,
                  inset: 3.5pt,

                  // vert([Number of rays]),
                  // table.cell(rowspan:2)[Num. Rays &\ Grid Shape],
                  table.header("r", "e", "a", $Delta t$),
                  table.hline(),
                  table.vline(x: 3),
                  [1], [0], [0], [.3],
                  [0], gt(0), gt(0), [.8],
                  gt(0), [-1], gt(0), [.7],
                  gt(0), gt(-1), [1], [.5],
                  gt(0), gt(-1), [2], [.5],
                  gt(0), gt(-1), [3], [.2],
                  gt(0), [0], gt(3), [.7],
                  [1], gt(0), gt(3), [#sym.infinity],
                  [0], gt(0), gt(3), [NaN],
                  [0], gt(0), gt(3), [NaN],
                  [], text(size:2em)[#sym.dots.h], [], [],
                  gt(0), gt(0), [7], [NaN],
              ),
          ),
          "sorting",
          "Forward-filled"
      ),
      subfigure(grid(
              columns: 3, align: bottom,
              table(
                  // columns: (auto, auto, auto, auto),
                  columns: 4,
                  stroke: none,
                  // inset: 10pt,
                  align: horizon,
                  inset: 3.5pt,

                  // vert([Number of rays]),
                  // table.cell(rowspan:2)[Num. Rays &\ Grid Shape],
                  table.header("r", "e", "a", $Delta t$),
                  table.hline(),
                  table.vline(x: 3),
                  ga("1"), [0], [0], ga("0"),
                  [0], gt(0), gt(0), [.8],
                  gt(0), ga("-1"), [ ], ga("0"),
                  gt(0), ga("-1"), [1], ga("0"),
                  gt(0), ga("-1"), [2], ga("0"),
                  gt(0), ga("-1"), [3], ga("0"),
                  gt(0), [0], gt(3), [.7],
                  [1], gt(0), gt(3), ga("0"),
                  [0], gt(0), gt(3), ga("0"),
                  [0], gt(0), gt(3), ga("0"),
                  [], text(size:2em)[#sym.dots.h], [], [],
                  gt(0), gt(0), [7], ga("0"),
              ),
          ),
          "sorting",
          "Zeroed"
      ),
  ),
      caption: [a) Sorted region indices are inserted into empty table along with voxel intersection lengths. b) Missing table entries are filled in from the previous row to form a complete 3D voxel index.  c) Rows with invalid entries are set to zero so the row has no effect during raytracing.]
  ) <alg_table>

=== Step 4 - Line Integration

    The final step of the raytracer is to perform a weighted sum of the object values at the voxel indices given in the table in the previous section with the $Delta t$ column as weights.

  // FIXME - use x instead of rho

    #math.equation(
        $sum_(i=0)^(2(N_r + 1) + \ 2(N_e + 1) + \ (N_a + 1) - 1) bold(rho)[r_i, e_i, a_i] * Delta t_i = innerproduct(rho, Delta t)$
    )

    where $bold(rho)$ is the object to be integrated and $r_i$, $e_i$, and $a_i$ are indices from the $i$th row of the table.


== API Overview <api_overview>
  We have designed this library to allow the user to easily construct tomography problems in just a few lines by instantiating classes that represent the grid, view geometries and forward raytrace operator.  These classes also support visualization via the ```python .plot()``` method to aid in validation and presentation.  In this section, we provide some some examples of setting up these classes for typical tomography problems as well as visualizations generated from the library.

=== Grid
    // FIXME: don't use theta/phi, switch to e/a
    The grid defines the physical extent and shape of the object that will be traced.  It is created using the `SphericalGrid` class and may be defined by providing `shape` and `size` arguments when a uniform spherical grid is desired, as in @api_grid.  Alternatively, one may manually specify the grid boundary locations via `r_b`, `e_b`, `a_b` arguments for a non-uniform grid.

    #figure(
        grid(columns: 2, align: horizon,
            rbox(
                ```python
                grid = SphericalGrid(
                  shape=(30, 30, 30), # voxels
                  size_r=(0, 10)
                )
                grid.plot()
                ```
            ),
            image("figures/tomosphero/api_grid.png", height: 8em)
        ),
        caption: [An origin-centered grid defined out to radius 10],
    ) <api_grid>

=== View Geometries
    A _view geometry_ is a logical collection of lines of sight.  This library provides a few built-in view geometries for common tomography paradigms shown in @api_conegeom and @api_parallelgeom.  In these examples, `pos` and `shape` control the detector position and shape.  Detector orientation is inferred automatically, but this behavior may be overriden.

    #figure(
        grid(columns: 2, align: horizon,
            rbox(
                // FIXME: need `fov` in this example
                ```python
                geom = ConeRectGeom(
                  pos=(10, 0, 0),
                  shape=(64, 64) # pixels
                )
                geom.plot()
                ```
            ),
            image("figures/tomosphero/api_conerectgeom.png", height: 10em)
        ),
        caption: [Cone beam view geometry with rectangular detector located at vertex],
    ) <api_conegeom>

    #figure(
        grid(columns: 2, align: horizon,
            rbox(
                // FIXME: need `size` in this example
                ```python
                geom = ParallelGeom(
                  pos=(10, 0, 0),
                  shape=(4, 4) # pixels
                )
                geom.plot()
                ```
            ),
            image("figures/tomosphero/api_parallelgeom.png", height: 10em)
        ),
        caption: [Parallel beam view geometry with rectangular detector],
    ) <api_parallelgeom>

    For more exotic detector geometries, one may define an entirely custom view geometry by specifying each ray manually via the `ViewGeom` class as in @api_customgeom or by subclassing `ViewGeom`.
    #figure(
        grid(columns: 2, align: horizon,
            rbox(
                ```python
                import torch as t
                rays = t.rand((4, 4, 3)) / 5
                rays[:, :, 0] = -1
                ray_starts = t.tensor((10, 0, 0)).broadcast_to(rays.shape)
                geom = ViewGeom(
                  rays=rays,
                  ray_starts=ray_starts
                )
                geom.plot()
                ```
            ),
            image("figures/tomosphero/api_customgeom.png", height: 10em)
        ),
        caption: [Custom view geometry for a 4-by-4 detector with random rays],
    ) <api_customgeom>
=== Composing View Geometries
    The view geometries given above can be composed into a `ViewGeomCollection` to create more complex geometries.  A `ViewGeomCollection` may represent a set of sensors which are taking measurements at the same instant or a single moving sensor making sequential measurements.  Programmatically, two view geometries can be combined by Python's addition operator, as in @api_composition.  In this case, `.plot()` returns an animated visualization showing the sequence of view geometries.

    #figure(
        // raw(lang: "python", "hello()"),
        // code2(raw(lang: "python", "hello()")),
        grid(columns: 2, align: horizon,
            rbox(
                ```python
                geom = None
                # sweep 360° around origin
                for w in t.linspace(0, 2*t.pi, 50):
                  pos=(5*t.cos(w), 5*t.sin(w), 1)
                  # combine geoms by addition
                  geom += ConeRectGeom(
                    pos=pos,
                    shape=(100, 100), # pixels
                    fov=(25, 25) # degrees
                  )
                anim = geom.plot()
                ```
            ),
            image("figures/tomosphero/api_collection.png", height: 9.5em),
        ),
        caption: [Circular orbit constructed from individual view geometries and a single frame from the resulting animation],
    ) <api_composition>

=== Operator

    When both view geometry and grid have been defined, the user may instantiate an `Operator` class, which carries out the precomputation steps described in @alg_outline and stores results in memory.  When called with a 3D (static) or 4D (dynamic) object array, this instantiated `Operator` returns raytraced measurements.  It is automatically differentiable on account of PyTorch's autograd capability, which is useful in reconstruction algorithms that require access to gradients.

    #figure(
        // raw(lang: "python", "hello()"),
        // code2(raw(lang: "python", "hello()")),
        grid(columns: 2, align: horizon,
            rbox(
                ```python
                op = Operator(grid, geom)
                # example object - broken torus
                rho = t.ones(grid.shape)
                rho[12:15, 12:18, 0:26] = 1
                # raytrace and get measurements
                y = op(rho)
                # y.shape matches geom.shape
                ```
            ),
            // image("figures/api_operator.png", height: 10em)
            image("figures/tomosphero/api_measurements.png", height: 10em)
        ),
        caption: [Operator constructed from previously defined view geometry and grid with associated measurements],
    ) <api_operator>


=== Reconstruction Examples

    The primary purpose of this library is as a building block for implementing experimental reconstruction algorithms. @api_reconstruction shows a simple example of a reconstruction algorithm which makes use of gradient descent to reconstruct the example object in the previous section.  The example minimizes the $L_2$ error between the measurements `y` and synthetic measurements of candidate reconstructed object `rho_hat`, but this loss may be combined with regularizers or any differentiable function of `rho_hat`.

    #figure(
        grid(columns: 2, align: horizon,
            rbox(
                ```python
                # choose initial guess for object
                rho_hat = t.zeros(
                  grid.shape, device='cuda',
                  requires_grad=True
                )
                # choose optimizer vars and method
                optim = t.optim.Adam(
                  [rho_hat], lr=1e-3,
                  weight_decay=1e2
                )
                # iterative reconstruction loop.
                # minimize L2 loss w.r.t. rho_hat
                for i in range(500):
                  loss = t.sum((y - op(rho_hat))**2)
                  loss.backward()
                  optim.step()

                y_hat = op(rho_hat)
                ```
            ),
            image("figures/tomosphero/api_reconstruction.png", height: 9em),

        ),
        caption: [GPU-enabled reconstruction algorithm using off the shelf optimizer from PyTorch and a squared error loss function],
    ) <api_reconstruction>

    If the reconstructed object lies in some low-dimensional transform domain, the user can provide a differentiable transformation function and instruct PyTorch to optimize the latent representation of the object `c_hat` instead of the object directly, as shown in @api_reconstruction_latent.

    #figure(
        rbox(
            ```python
            def f(c):
              # map from low-dimensional latent space to R³.
              # use PyTorch-differentiable functions here
              ...
              return rho
            # choose initial guess for latent variables
            c_hat = t.zeros(...)
            # choose optimizer vars and method
            optim = t.optim.Adam([c_hat], ...)
            # iterative reconstruction loop.
            # minimize L2 loss w.r.t. c_hat
            for i in range(500):
              loss = t.sum((y - op(f(c_hat)))**2)
              ...
            ```
        ),
        caption: [Pseudocode for reconstructing an object which lies in some transform domain defined by `f`],
    ) <api_reconstruction_latent>

    Finally, we remark that the differentiability of the operator allows it to be exploited in machine learning models, such as _physics-informed neural networks_ (PINNs), where the operator forms layers of the network itself @pinn.
=== Assumptions and Limitations

    Unlike most medical imaging paradigms, the detector is located at the vertex of `ConeRectGeom`.  We also assume the detector lies outside of the grid, rays extend to infinity, and that the object is zero outside of the grid.

    When combining view geometries into a collection, the geometries must have the same shape so that the corresponding stack of measurements returned during raytracing is a rectangular array.
    In the case of a dynamic object, the grid of the object being raytraced cannot move nor can its shape or extent change.

== Validation and Tests

    The correctness of the raytracer has been thoroughly tested both by visual verification of various test objects from a variety of positions as in @visual_validation as well as comparing individually raytraced lines of sight against known analytic results.  Additionally, we have tested each phase of the raytracer algorithm outline described in @alg_outline, including the steps for deriving boundary intersection points, boundary crossing direction, and intersection point sorting.  These tests are integrated into the raytracer package's automated testing framework for evaluation when new software releases are made.

    #figure(
      grid(columns: 2, column-gutter: 1pt,
          subfigure(image("figures/limb_darkening.png", height: 10em), "samplert", "Solid sphere with limb darkening"),
          subfigure(image("figures/checkerboard.png", height: 10em), "samplert", "Hollow checkerboard shell"),
      ),
      caption: [Raytraced sample objects],
    ) <visual_validation>

  #rt([
  // FIXME: WIP
  - FIXME: more detail about validation? - I ran many more tests comparing column densities computed analytically to raytracer output which I could describe here, (figure) testing different LOS cases

  ])

== Operator Memory and Timing Benchmarking <benchmarking>

    Since precomputed arrays are held in memory, it is important to know the peak memory usage of the algorithm to determine whether a particular raytracer forward operator will fit in memory.  If running the raytracer on a GPU, memory is even more limited, with most consumer graphics cards constrained to 12GB or 24GB of VRAM as of 2024.

    The peak memory required for steps 1-4 can be computed in gigabytes approximately with @memusage.

    #math.equation($#text[Peak GB] approx 6N#sub[rays] (2N#sub[r] + 2N#sub[e] + N#sub[a]) (8) \/ 10#super[9]$) <memusage>

    A breakdown of the terms of this expression is below:

    - $2N#sub[r] + 2N#sub[e] + N#sub[a]$ - array length (maximum number of voxels intersecting a ray)
    - $6N#sub[rays]$ - 6 arrays for every ray containing
      - voxel $[r, e, a]$ indices
      - voxel intersection lengths
      - voxel values
      - temporary array for intermediate computation
    - 8 - element size of int64 and float64 arrays, in bytes
    - 10#super[9] - conversion from bytes to gigabytes

    @memtable gives an overview of peak memory use and compute time for several configurations of grid and view geometry shapes.
    We conducted these tests on Ubuntu 22.04 with AMD Threadripper 5995WX (128GB RAM) and NVIDIA RTX 4070 (12GB VRAM) with CUDA 12.4 and PyTorch 2.2.2.  Note that additionally memory will be consumed by autograd and the runtime and may be platform-dependent.

    #figure(
        table(
            columns: (auto, 40pt, 40pt, 40pt, 40pt, auto),
            // columns: 6,
            // inset: 10pt,
            align: horizon,

            // vert([Number of rays]),
            // table.cell(rowspan:2)[Num. Rays &\ Grid Shape],
            table.header(
                // table.cell(rowspan:2)[#rotatex([Num. Rays &\ Grid Shape], -90deg)],
                table.cell(rowspan:2)[Num. Rays &\ Grid Shape],
                table.cell(colspan:2, y:0, x:1)[Precompute Time\ (steps 1-3)],
                table.cell(colspan:2, y:0, x:3)[Raytrace Time\ (step 4)],
                table.cell(rowspan:2)[#vert([Peak Memory])],
                vert([CPU]),
                vert([GPU]),
                vert([CPU]),
                vert([GPU]),
            ),
            table.hline(stroke: 2pt),
            [1x256x256\ (100, 100, 100)], [0.88s], [0.28s], [42ms], [1.5ms], [1.5GB],
            [1x512x512\ (50, 50, 50)], [1.5s], [0.17s], [82ms], [3.1ms], [2.9GB],
            [30x64x64\ (50, 50, 50)], [0.75s], [0.093s], [37ms], [1.4ms], [1.5GB],
        ),
        caption: [Resource usage for various grid and view geometry shapes.]
    ) <memtable>


= Static Retrieval of Exosphere <static_retrieval>
  #rt([
    - [x] intro
        - numerous retrieval approaches to limtied exospheric data
        - some simple (1D estimate, few parameters), some more complex (HDOF)
        - summarizes point of this chapter
        - purpose of chapter is to introduce historical static retrieval approaches (generally increasing in complexity), then finally introducing a new static model contributed to this thesis



    - [x] 1D retrievals
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
    - [x] gonzalo MAP estimate
        - derived from earlier work from Butala
        - not enough viewing diversity, stereoscopic
    - [x] Zoennchen 2024
        - functional form based on spherical harmonics
        - defined at single shell with exponential decay (and other terms)
    - [x] spherical harmonic model
        - applicability of spherical harmonic bases to modelling exospheres
            - figure: show direct fits for different L and the max error on each
            - conclusion: spherical harmonic basis is a good low rank representation of exospheric models
        - splines enforce smoothness
            - agnostic shape
        - figure: spherical bases
        - single measurement A00 retrieval
    - [ ] longterm: summary table?
  ])

  This chapter presents several historical and recent methods of static retrieval of hydrogen density.  The retrieval methods are organized roughly in increasing complexity, ending with a novel method in @spline_model.  Mathematical notation differs from the original publications to conform with notation used in this manuscript (see @inverse_problem).

  == 1D Retrievals <1d_retrieval>

    Early retrievals often relied on simple spherically symmetric 1D models of the exosphere.  In cases where geocoronal studies were taken oppooortunistically (e.g. Galileo Earth flyby), measurements are often only available from a single vantage which has limited viewing geometry diversity.
    The assumption of spherically symmetry naturally produces a well-conditioned inverse problem without requiring significant measurement diversity and avoids an underdetermined system by keeping model dimensionality low, both important for single snapshot retrievals.

    A fundamental contribution to the field is the Chamberlain model, which is a spherically symmetric model derived from knowledge of motion of H atoms in the upper atmosphere (section XXX), (#rt([FIXME: section reference])).  Under Liouville's theorem, the model makes some simple assumptions about the trajectories followed by exospheric H atoms and derives a density distribution from the resulting PDE.
    Baliukin et al. @baliukin applied the Chamberlain model to measurements from the SOHO/SWAN mission, extending it to take into account solar radiation pressure, loss of Lyman-α-scattering H atoms over time due to ionization, and prior knowledge of H density at lower altitudes from empirical models (e.g. NRLMSIS). (#rt([FIXME: citation, does Chamberlain also depend on exobase knowledge?])).  The retrieval utilizes an onion-peeling technique where layers of the exosphere are retrieved one-at-a-time starting at the outermost layer and moving inwards.

    Other techniques include parametric 1D models found empirically and fit to data.  Østgaard et al. @ostgaard use a double-exponential model following the form

    #math.equation(
        $m(bold(c)) = n_H (r) = n_1 e^(-r/alpha_1) + n_2 e^(-r/alpha_2)$
    )

    where $bold(c) = {n_1, n_2, alpha_1, alpha_2}$ are unknown model coefficients.  The two exponentials express a belief that the H density distribution consists of two populations of H density atoms.  Additionally, this simple exponential formulation permits an analytic line integration that Østgaard et al. utilize to perform model fits directly in the measurement domain with no tomographic operator.  The exact method of fitting is not described, but most likely a linear regression similar to

    #math.equation(
        $
            bold(hat(c)) = arg min_bold(c) ||bold(y) - bold(y)_bold(c)||_2^2 \
            bold(hat(c)) = m(bold(hat(c)))
        $
    )

    where $bold(y)_bold(c)$ are column densities derived by analytically integrating $m(bold(c))$.

    These 1D solutions can be expanded to more sophisticated models shown in the next few sections when more measurement data is available from orbits designed to provide better view geometry diversity.

  == Spherical Harmonic Representation <sph_power>

      Zoennchen et al. @zoennchen_new present a retrieval method (derived from earlier work @zoennchen2011, @zoennchen2013, @zoennchen_old) featuring a partially separable parametric model based on a spherical harmonic representation (SHR).
      Spherical harmonics are a family of complex, two-dimensional, continuous, orthonormal functions defined over a spherical domain.  These functions are commonly used in planetary atmospheric sciences @hodges, atomic physics (electron orbitals), computer graphics and more.  They serve as an analog to the Fourier basis in a spherical domain and approximate smoothly-varying spherical functions by linear combination of the basis.

      Spherical harmonics are typically denoted $Y_(l m)(e, a)$ where $l$ and $m$ are known as the function _degree_ and _order_.  In planetary physics, the convention is to denote coefficients of this linear combination as $A_(l m)$ (when $m ≥ 0$) and $B_(l m)$ (when $m < 0$).  @sph_harm gives a graphical overview of the real/imaginary parts of the first few degrees of spherical harmonic functions.

      #figure(
          image("figures/harm.svg", height: 15em),
          caption: "First 3 degrees of spherical harmonic functions"
      ) <sph_harm>


      The model density distribution is given as

      #math.equation(
          $n_H (r, e, a) = c r^(-k) d^(1/r) "SHR"(r, e, a)$
        )

      with

      #math.equation($
          "SHR"(r, e, a) = sqrt(4 pi) sum_(l=0)^3 sum_(m=0)^l
           A_(l m)(r) dot.op "Re"[Y_(l m)(e, a)] + B_(l m)(r) dot.op "Imag"[Y_(l m)(e, a)] \
          A_(l m) = a_(l m) + b_(l m) ln(r) \
          B_(l m) = p_(l m) + q_(l m) ln(r)
      $
      ) <scipy_sph>

      where $bold(c) = {{a_(0 0), ...}, {b_(0 0), ...}, {p_(0 0), ...}, {q_(0 0), ...}, c, k}$ are model parameters.  Note that $∀l : B_(l 0)(r) = 0$ by definition.

      The formulation of $"SHR"(r, e, a)$ in @scipy_sph differs slightly from @zoennchen_new so that $Y_(l m)$ corresponds directly to Scipy's `sph_harm_y` function.

        Rewriting the model using this manuscript's notation and taking $m(bold(c)) = n_H (r, e, a)$ gives

      #math.equation(
          $
              hat(bold(c)) = arg min_bold(c) ||bold(y) - f(m(bold(c)))||_2^2 + lambda ||D_e m(bold(c))||_2^2 + lambda ||D_a m(bold(c))||_2^2 \
              hat(bold(x))_"HDOF" = m(hat(bold(c)))
          $
      )

      where $D_e$ and $D_a$ are finite difference operators for Tikhonov regularization that enforces elevational and azimuthal smoothness.

    == High Degree-of-Freedom

      An alternative approach presented in Zoennchen et al. @zoennchen_new  @gonzalolaica and based on @solartomography1 is the #gls("HDOF") method which uses a non-parametric model where each voxel in the underyling density discretization is a separately optimizable parameter.  This method utilizes #gls("MAP") estimation which is optimal in a probabilistic sense but requires assumptions about the density distribution and noise statistics.  Specifically this method assumes that hydrogen density is a Gaussian Markov Random Field (GMRF) distributed according to a known prior distribution $cal(N)(bold(x)_"pr", Σ_"pr")$ with mean and covariance $bold(x)_"pr"$ and $Σ_"pr"$ provided externally (e.g. derived from other datasets #rt("FIXME: which dataset?")) and measurements follow a noise distribution $Σ_"e"$.

      Under these assumptions, the #gls("MAP") solution has closed form

      #math.equation(
          $hat(bold(x))_"MAP" = (L^T Σ_e^(-1) L + Σ_"pr"^(-1))^(-1)(L^T Σ_e^(-1) bold(y) + Σ_"pr"^(-1) x_"pr")$
      ) <gonzalo_map>

      where $L$ is a matrix representation of forward operator $f$.

      An alternative formulation of @gonzalo_map is

      #math.equation(
          $
              hat(bold(x))_"MAP" = bold(x)_"pr" + underbrace(Σ_"pr" L^T (L Σ_"pr"L^T + Σ_e)^(-1), K)(bold(y) - L bold(x)_"pr")
          $
      ) <gonzalo_map_alt>

      // #math.equation(
      //     $
      //         hat(bold(x))_"MAP" = bold(x)_"pr" + K(bold(y) - L bold(x)_"pr") \
      //         K = (L^T Σ_e^(-1) L + Σ_"pr"^(-1))^(-1) L^T Σ_e^(-1)
      //     $
      // ) <gonzalo_map_alt>

      which illustrates that the #gls("MAP") solution is the density prior which has been shifted by residual measurement error $bold(y) - L bold(x)_"pr"$ projected back on to the density domain.

      The matrix $K$ serves as a weighted pseudoinverse which balances between using prior information from $bold(x)_"pr"$ versus measurement data $bold(y)$ based on uncertainty in the prior ($Σ_"pr"$) and measurement noise ($Σ_e$).

      Unfortunately, the formulation given in @gonzalo_map_alt is not always computationally tractible, as simply storing the result of $L Σ_"pr" L^T$ for a relatively small forward operator requires excessive memory (e.g. 10·64² LOS, 50³ voxels → 125 GB).  Zoennchen solve this by observing that $Σ_"pr"^(-1)$ is sparse given locality assumptions of a GRMF (distant voxels are independent), and full computation of a dense matrix is unnecessary.  Another approach is to write the solution as a minimization problem @icon_inversion like

      #math.equation(
          $
              hat(bold(x))_"MAP" = arg min_(bold(x)) ||bold(y) - L bold(x)||_(Σ_e^(-1))^2 + ||bold(x) - bold(x)_"pr"||_(Σ_"pr"^(-1))^2
          $
      )

      or using this manuscript's notation

      #math.equation(
          $
              hat(bold(c)) = arg min_(bold(c)) ||bold(y) - f(m(bold(c)))||_(Σ_e^(-1))^2 + ||m(bold(c)) - bold(x)_"pr"||_(Σ_"pr"^(-1))^2 \
              hat(bold(x))_"MAP" = m(bold(hat(c)))
          $
      )

      where model $m$ is the identity operator for this non-parametric case.

      This minimization problem can be solved using a number of iterative methods, such as conjugate gradient or automatic differentiation (utilizing the differentiable forward operator presented in @raytracer).

  == Robust, Regularized, Positive Estimation

    Cucho-Padin et al., present another non-parametric method called #gls("RRPE") @rrpe @gonzalorrpe where voxels are allowed to vary independently. This model uses several 1D Tikhonov regularizers to enforce smoothness in each dimension.  The optimization problem is written

    #math.equation(
        $
            hat(bold(c)) = arg min_bold(c) ||bold(y) - f(m(bold(c)))||_2^2
            + lambda_r ||D_r m(bold(c))||_2^2
            + lambda_e ||D_e m(bold(c))||_2^2
            + lambda_a ||D_a m(bold(c))||_2^2 \
            hat(bold(x))_"RRPE" = m(hat(bold(c)))
        $
    )

    where $D_r$, $D_e$ and $D_a$ are finite difference operators for radial, elevational and azimuthal directions, and $m$ is identity operator.  Cucho-Padin et al. use generalized cross validation to derive values for regularizer hyperparameters.

  == Spherical Harmonic Spline Model <spline_model>

    This section describes a partially separable model which uses spherical harmonics as basis functions like in @sph_power.  However, instead of an inverse power law for controlling radial variation of spherical harmonic coefficients, this method samples the coefficients from cubic spline functions.  This formulation is less strict about the shape of radial density profiles but still enforces smoothness of spherical harmonic coefficients.

    The model can be written as

    #math.equation(
        $
            m(bold(c)) &= sum_(l=0)^L sum_(m=-l)^l S_(l m)(r) X_(l m)(e, a) \
            S_(l m)(r) &= sum_(k=1)^K c_(l m k) B_(k 2)(r)
        $
    )

    #math.equation($
        X_(l m)(e, a) = cases(
        "Re"[Y_(l m)(e, a)] "if" m ≥ 0,
        "Im"[Y_(l m)(e, a)] "if" m < 0,
        )
      $)

    where $S_(l m)(r)$ is a cubic spline function composed of $K$ third order B-splines ${B_(1 2), ...}$ that defines spherical harmonic coefficients for any radial shell from a small number of control points.  This formulation of basis function $X_(l m)$ is equivalent to @sph_power but avoids notational awkwardness when $m ≥ 0$ versus $m < 0$.

    The number and location of the $K$ control points used in the spline functions and the maximum spherical harmonic degree $L$ are fixed and selected before optimization.
    $bold(c) = {c_(0 0 0), ...}$ are the model parameters and the minimization problem is

    #math.equation(
        $
            hat(bold(c)) = arg min_(bold(c)) 1 / (|y|) ||y - f(m(bold(c)))||_1 + lambda  1 / (|m(bold(c))|) ||"clip"_(-infinity, 0)(m(bold(c)))||_1 \
            bold(hat(x)) = m(bold(hat(c)))
        $
    )

    #rt("FIXME: need to justify choice of L1 instead of L2 here.")

    Note that a choice of $L=0$ enforces spherical symmetry, which can be useful for single snapshot retrievals like those presented in @1d_retrieval.

    === Implementation Notes

      - Basis functions ${X_(l m), ...}$ may be computed once during initialization and used for all grid radial shells
      - Dimensions $l$ and $m$ should be flattened and merged for ${X_(0 0), ...}$ and ${c_(0 0 0), ...}$ to avoid an awkward pyramidal array structure
      - affine map for log-spaced control points

      #rt([FIXME: summary table])

      #table(
          columns: (auto, auto),
      table.header([*Method*], [*Parametric*]),
          [Chamberlain], [?],
          [Østgaard], [Yes],
          [SHR], [Yes],
          [HDOF], [No],
          [RRPE], [No],
          [Spline], [Yes],
      )


= Static Retrieval Validation
  #rt([
  - reconstruction requirements
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
  ])

  #figure(
      image("figures/retrieval_overview.png", width: 90%),
      caption: "Diagram of simulator and retrieval loop used during validation"

  ) <codeoverview>

  #rt("needs citations for each row")

  #figure(
      table(
          columns: 5,
          inset: 4pt,
          align: horizon,
          table.header(
              "Dataset", [Dynamic/\ Static], "Derived from", "Storm Condition", "Source",
          ),
          table.hline(stroke: 2pt),
          [Zoennchen], [Static], [TWINS], [Quiet], [Zoennchen 2024],
          [Cucho-Padin], [Static], [TWINS], [Quiet], [Zoennchen 2024],
          [Pratik], [Dynamic], [NRLMSIS], [Quiet], [Pratik Joshi],
          [Connor], [Dynamic], [Physics-based], [Storm/Quiet], [],
          [Cucho-Padin\ Dynamic], [Dynamic], [TWINS], [], []
      ),
      caption: "Storm Time"
  ) <datasets>

  #figure(
      table(
          columns: 4,
          table.header(
              "Dimension", "Range", "# of Bins", "Spacing",
          ),
          table.hline(stroke: 2pt),
          [time], [NA], [NA], [linear, 1 hr],
          [radius], [3-25 Re], [200], [logarithmic],
          [elevation], [0-180 deg.], [60], [linear, 3 deg.],
          [azimuth], [0-360 deg.], [80], [linear, 4.5 deg.],
      ),
      caption: "Storm Time"
  ) #label("stormbins")



  #figure(
      table(
          columns: 4,
          table.header(
              "Dimension", "Range", "# of Bins", "Spacing",
          ),
          table.hline(stroke: 2pt),
          [time], [NA], [NA], [linear, 6 hr],
          [radius], [3-25 Re], [200], [logarithmic],
          [elevation], [0-180 deg.], [60], [linear, 3 deg.],
          [azimuth], [0-360 deg.], [80], [linear, 4.5 deg.],
      ),
      caption: "Quiet Time"
  ) <quietbins>

= Dynamic Retrieval of Exosphere

= Dynamic Retrieval Validation

= Conclusion

  #set heading(numbering: "A.1.1", supplement: [Appendix])
  #counter(heading).update(0)

= Appendix - Instrument Simulation and Background Subtraction <appendix_sim>

  This section describes a statistical model for the instrument noise and background signals present in the NFI and WFI cameras during measurement of exospheric Lyman-α.  Modelling these processes is important for converting raw sensor measurements in digital numbers (DN) as telemetered by the spacecraft to corresponding radiances that can be used for tomographic reconstruction, summarized in @calibration.   A statistical model is also important for generating synthetic noisy measurements to validate the performance of retrieval algorithms.  As a result, the Carruthers cameras have undergone extensive laboratory characterization to determine instrument model parameters and periodic on-orbit calibration is planned to account for parameter drift due to exposure to the space environment.

  == Instrument Model

  The Carruthers spacecraft contains two #gls("UV") cameras designed for detecting exospheric Lyman-α emission, an indicator of the amount of atomic hydrogen present along a #gls("LOS").  These cameras utilize a design which has heritage with other #gls("UV") instruments such as ICON @icon and GUVI @guvi and contain the following stages @rider:

  // TODO: WIP

  - *optical filter wheel* - provides multispectral measurement capability for out-of-band signal measurement used in thick exosphere retrievals
  - *KBr photocathode* - potassium bromide photocathode for converting #gls("UV") photons to electrons (known as _photoelectrons_)
  - *#gls("MCP")* - amplifies individual photoelectrons to a detectable shower of electrons
  - *phosphor and #gls("CCD")* - phosphor produces light pulses from electrons which are detected by a #gls("CCD") and converted to an electrical charge
  - *ADC* - #gls("ADC") for reading out #gls("CCD") charge.  (together with the #gls("CCD") this is sometimes referred to as an #gls("APS"))



  The rest of this section will consist of a derivation of a single pixel noisy measurement in DN given a photon spectral radiance and other quantities in @knownvariables and @randomvariables.

  Let $L_("exo")(lambda) + L_("bkg")(lambda)$ represent exospheric photon spectral radiance and background.
  Using the etendue of a single pixel $Omega dot.op A$, it is possible to convert photon spectral radiance to a photon spectral flux
  #math.equation(numbering: none,
      $p_("phot")(lambda) = (L_("exo")(lambda) + L_("bkg")(lambda)) dot.op Omega dot.op A gt("(phot/s/nm)")$
  )
  After entering the instrument, the photons pass through an #rt([optical filter and photocathode where they are converted to a stream of photoelectrons (also known as _events_) with rate $e_"phot"$.  This conversion happens with efficiency $f_"opt" epsilon(lambda)$ (events/phot)]) where unitless scaling factor $f_"opt"$ (known as a _flat-field_) accounts for spatial variation in the optics, filters and photocathode.

  #math.equation(numbering: none,
      $e_"phot" = f_"opt" integral_lambda p_("phot")(lambda) epsilon(lambda) dif lambda gt("(events/s)")$
  )

  Alternatively, if $L_"exo"$ and $L_"bkg"$ are monochromatic with intensity $I_"exo"$ and $I_"bkg"$ and wavelength $lambda_0$

  #math.equation(numbering:none,
      $e_"phot" &= f_"opt" integral_lambda p_("phot")(lambda) epsilon(lambda) dif lambda \
          &= f_"opt" (integral_lambda L_("exo")(lambda) epsilon(lambda) dif lambda + integral L_("bkg")(lambda)  epsilon(lambda) dif lambda) dot.op Omega dot.op A \
          &= f_"opt" epsilon(lambda_0) (I_"exo" + I_"bkg") dot.op Omega dot.op A gt("(events/s)")$

  )

  Assuming Lyman-α photons are emitted in a Poisson process, the actual number of photoelectrons created in a single camera frame of duration $t_"fr"$ is given by random variable

  #math.equation(numbering:none,
      $E_"phot" tilde.op "Pois"(t_"fr" e_"phot") gt("(events)")$
  )

  In general, single photoelectrons are difficult to detect, so the Carruthers cameras employ an #gls("MCP") for turning a single particle into a detectable shower of particles.  #glspl("MCP") consist of an array of small glass tubes (channels) which are electrically charged so that a photoelectron striking the wall of one of these tubes will cause a cascade of particles via secondary emission @microchannelplate.  The small size of these channels ensures that the subsequent shower of particles exits the #gls("MCP") in the same location as the photoelectron.  Due to the nature of secondary emission, the number of particles created by the #gls("MCP") from a photoelectron is given by the discrete random variable $G_"mcp"$, whose distribution has been measured in the laboratory.  The number of electrons leaving the MCP #rt([(known as _counts_)]) due to Lyman-α photoelectrons is given by

  #math.equation(numbering: none,
      $f_"mcp" sum_(l=1)^(E_"phot") G_"mcp,l" gt("(counts)")$
  )

  where $f_"mcp"$ is a flat-field representing spatial gain variation on the #gls("MCP").

  Energetic particles from space which directly penetrate the spacecraft body and impinge on the MCP are modeled by $E_"mcp" tilde.op "Pois"(t_"fr" e_"mcp")$ and the number of counts is given by

  #math.equation(numbering: none,
      $f_"mcp" sum_(l=1)^(E_"mcp") G_"mcp,l" gt("(counts)")$
  )

  Additional sources of electron counts are energetic particles impinging on the #gls("APS") ($C_"aps"$), as well as thermally generated electrons from the instrument electronics known as _dark current_ ($C_"dark"$).  The charge from these electrons accumulates in the #gls("CCD") until it is amplified and read out by the #gls("ADC") as #gls("DN") given by

  #math.equation(numbering: none,
      $Y &= (f_"mcp" sum_(l=1)^(E_"phot" + E_"mcp") G_"mcp,l" + C_"aps" + C_"dark") g_"aps" + R gt("(DN)")$
  )


  where $R$ is normally distributed read noise with bias introduced during readout and $g_"aps"$ is configurable gain in the #gls("APS") electronics.


  @bd_annotated shows the Lyman-α signal at each stage in the camera, starting with photon spectral flux and ending with a final measurement in DN.
A summary of all variables and sources of randomness is given in @knownvariables and @randomvariables.

  #figure(
      image("figures/bd_annotated2.svg"),
      caption: [Statistical model of instrument, showing instrument stages and sources of background noise]
  ) <bd_annotated>

  #figure(
      table(
          columns: 3,
          inset: 0.5em,
          table.header(
              "Variable", "Description", "Units",
          ),
          table.hline(stroke: 2pt),
          [$L_("exo")(lambda)$],  [exosphere photon spectral radiance], [phot/cm²/sr/s/nm],
          [$L_("bkg")(lambda)$], [background photon spectral radiance], [phot/cm²/sr/s/nm],
          [$A$], [pixel area], [cm²],
  [$Omega$], [pixel solid angle], [sr],
  [$epsilon(lambda)$], [optical efficiency], [events/phot],
  [$e_("mcp")$], [MCP radiation rate], [events/s],
  [$f_("mcp")$], [MCP gain flat-field], [unitless],
  [$f_("opt")$], [optical efficiency flat-field], [unitless],
  [$c_("aps")$], [APS radiation rate], [counts/s],
  [$c_("dark")$], [dark count rate], [counts/s],
  [$b$], [ADC bias], [DN],
  [$t_("fr")$], [frame integration time], [s],
  [$g_("aps")$], [APS gain], [DN/count],
  [$N$], [APS binning factor], [unitless],
  [$M$], [number of frames to stack], [frames],
      ),
      caption: "Known or derived quantities",
  ) <knownvariables>


  #figure(
      table(
          columns: 3,
          inset: 0.5em,
          table.header(
              "Variable", "Description", "Units",
          ),
          table.hline(stroke: 2pt),
  [$E_"phot" tilde.op "Pois"(t_"fr" e_"phot")$], [Lyman-α events], [events],
          [$E_"mcp" tilde.op "Pois"(t_"fr" e_"mcp")$], [MCP radiation], [events],
          [$G_"mcp" tilde.op cal(Z)("mean"=mu_G, "var"=sigma^2_G)$], [MCP gain distribution], [counts/event],
          [$C_"aps" tilde.op "Pois"(t_"fr" c_"aps")$], [APS radiation], [counts],
          [$C_"dark" tilde.op "Pois"(t_"fr" c_"dark")$], [Dark current], [counts],
          [$R tilde.op cal(N)(b, sigma^2)$], [Read noise and bias], [DN],
          ),
      caption: "Random Variables",
  ) <randomvariables>

    == Frame Stacking, Binning, and Simulation Approximation

    A common technique to increase measurement #gls("SNR") in spaceborne telescopes is to increase exposure time at the expense of temporal resolution.  As frames are read out by the #gls("APS"), they are accumulated digitally on the spacecraft in a process known as _frame stacking_.  This effectively achieves the #gls("SNR") of a single, long exposure measurement while avoiding the problem of charge well saturation in the #gls("CCD") and has the added benefit of reducing bandwidth requirements by only downlinking a single frame-stacked measurement.  The loss of temporal resolution is negligible relative to the rate of evolution of the exosphere (#rt("FIXME: citation?")).

    Likewise, spatial resolution can be sacrificed to boost measurement #gls("SNR") by _pixel binning_ adjacent pixels on the detector together.  Binning may be performed digitally, but additional read noise reduction can be achieved in the #gls("APS") if pixels can be merged before being read out by the #gls("ADC"). @binningtype demonstrates the effect of binning on read noise variance for different detector types, with #gls("CCD") detectors achieving a factor $N^2$ reduction in read noise variance against purely digital binning.

    #figure(
      grid(columns: 2,
          subfigure(image("figures/nfi_1.png", height: 12em), "fs", "NFI, single frame"),
          subfigure(image("figures/wfi_1.png", height: 12em), "fs", "WFI, single frame"),
          subfigure(image("figures/nfi_14400.png", height: 12em), "fs", "NFI, stacked"),
          subfigure(image("figures/wfi_28800.png", height: 12em), "fs", "WFI, stacked"),
      ),
      caption: [Comparison of uncalibrated single frame vs uncalibrated framestacked measurement],
    ) <framestacking>

    @framestacking gives shows examples of simulated measurements which utilize both frame stacking and pixel binning and compares them to a single, unbinned frame.  More features of the underlying radiance distribution become visible when many frames are stacked.
    @integrationtime gives an overview of the stacking and binning parameters selected in order to meet mission #gls("SNR") requirements and maintain an acceptable loss of spatial and temporal resolution given prior knowledge of the exosphere @gonzalostorm. #rt("(clarify this is not a contribution of the thesis?)")

    #figure(
        table(columns: 4, align: horizon,
            table.header("Binning\nType", "Derivation\n(2x2 binning example)", "", "Read Noise\n(NxN binning)"),
            table.hline(stroke: 2pt),

            "Digital",
            math.equation(numbering:none, $
                x_1, &..., x_4 tilde.op cal(N)([mu_1, ..., mu_4], sigma^2) \
                y &= x_1 + ... + x_4 \
                    &tilde.op cal(N)(mu_1 + ... + mu_4, 4sigma^2)
            $),
            image("figures/bin_digital.png", height: 10em),
            $R tilde.op cal(N)(..., N^2 sigma^2)$,

            "CMOS",
            math.equation(numbering:none, $
                x_a &tilde.op cal(N)(mu_1 + mu_3, sigma^2) \
                x_b &tilde.op cal(N)(mu_2 + mu_4, sigma^2) \
                y &= x_a + x_b \
                    &tilde.op cal(N)(mu_1 + ... + mu_4, 2sigma^2)
            $),
            image("figures/bin_cmos.png", height: 7.5em),
            $R tilde.op cal(N)(..., N sigma^2)$,

            "CCD",
            math.equation(numbering:none, $
                x &tilde.op cal(N)(mu_1 + ... +  mu_4, sigma^2) \
                y &= x
                    tilde.op cal(N)(mu_1 + ... + mu_4, sigma^2)
            $),
            image("figures/bin_ccd.png", height: 7.5em),
            $R tilde.op cal(N)(..., sigma^2)$,
        ),
        caption: "Binned read noise for different detector binning types"
    ) <binningtype>

  #figure(
      table(
          columns: 3,
          inset: 4pt,
          align: horizon,
          table.header(
              [Parameter], [NFI], [WFI],
          ),
          table.hline(stroke: 2pt),
          [Total Integration Time (minutes)], [30], [60],
          [Stacked Frames (frames)], [14400], [28800],
          table.hline(stroke: 2pt),
          [Binning Factor (unitless)], [2x2], [4x4],
          [Downlinked Resolution (pixels)], [1024x1024], [512x512],
      ),
      caption: "Camera Integration Time and Frame Stacking"
  ) <integrationtime>

    Let subscripts $m$, $n$, and $l$ denote i.i.d. copies of the random variables defined in the previous section.

    - $m$ - frame
    - $n$ - APS pixel
    - $l$ - photon

    Then we can write a frame-stacked and binned measurement of $M$ individual frames and $N^2$ pixels as

    #math.equation(
    $Z &= sum_(m=1)^M sum_(n=1)^(N^2) Y_(m,n) \
        &= sum_(m=1)^M sum_(n=1)^(N^2) (f_("mcp") sum_(l=1)^(E_("phot",m,n) + \ E_("mcp",m,n)) G_("mcp",m,n,l) + C_("aps",m,n) + C_("dark",m,n) ) g_"aps" + R_(m,n) gt("(DN)")$
    ) <stackingequationdirect>

    Generating noisy measurements directly using the equation given above is computationally expensive due to the triply-nested summation, large number of stacked frames and large number of photoelectrons.  A more practical approach is to make use of the #gls("CLT") to approximate this summation of i.i.d. random variables with a Normal distribution.

    The central limit theorem states that as $M$ and $N$ go to infinity, random variable $Z$ is well-approximated by

    #math.equation(numbering:none,
        $Z tilde.op cal(N)(M N^2mu_Y, M N^2sigma_Y^2)$
    )

    where

    #math.equation(numbering:none,
        $mu_Y &= ex(Y) \
        &= (f_"mcp" ex(sum_(l=1)^(E_"phot" + E_"mcp") G_("mcp",l)) + ex(C_"aps") + ex(C_"dark") ) g_"aps" + ex(R) \
        &= (f_"mcp" ex(sum_(l=1)^(E_"phot" + E_"mcp") G_("mcp",l)) + t_"fr" c_"aps" + t_"fr" c_"dark" ) g_"aps" + b \
        & gt("See" #ref(label("appendix_rvsum")) "for derivation of expectation of sum") \
        &= (f_"mcp" mu_G t_"fr" (e_"phot" + e_"mcp") + t_"fr" c_"aps" + t_"fr" c_"dark" ) g_"aps" + b \
        &= t_"fr" (f_"mcp" mu_G (e_"phot" + e_"mcp") + c_"aps" + c_"dark" ) g_"aps" + b \
        sigma_Y^2 &= var(Y) \
        &= var( f_"mcp" sum_(l=1)^(E_"phot" + E_"mcp") G_("mcp",l) + C_"aps" + C_"dark" ) g_"aps"^2 + var(R) \
        &= ( f_"mcp"^2 var(sum_(l=1)^(E_"phot" + E_"mcp") G_("mcp",l)) + var(C_"aps") + var(C_"dark") ) g_"aps"^2 + var(R) \
        &= ( f_"mcp"^2 var(sum_(l=1)^(E_"phot" + E_"mcp") G_("mcp",l)) + t_"fr" c_"aps" + t_"fr" c_"dark" ) g_"aps"^2 + sigma^2 \
        & gt("See" #ref(label("appendix_rvsum")) "for derivation of variance of sum") \
        &= ( f_"mcp"^2 ex(G_"mcp"^2) t_"fr" (e_"phot" + e_"mcp") + t_"fr" c_"aps" + t_"fr" c_"dark" ) g_"aps"^2 + sigma^2 \
        &= t_"fr" ( f_"mcp"^2 ex(G_"mcp"^2) (e_"phot" + e_"mcp") + c_"aps" + c_"dark" ) g_"aps"^2 + sigma^2
        $
    )

    To validate the accuracy of the approximation, we compare a Monte Carlo simulation of @stackingequationdirect to a Normal distribution with parameters given above by the #gls("CLT"), shown in @clt. (#rt("FIXME: give specific parameters used in simulation below?  at least mention if photon rate is worst case"))

    #figure(
      grid(columns: 2, column-gutter: 1pt,
          subfigure(image("figures/montecarlo_1_NFI.png"), "mcclt", [NFI, 2x2 binning, 1 frame]),
          subfigure(image("figures/montecarlo_1_WFI.png"), "mcclt", [WFI, 4x4 binning, 1 frame]),
          subfigure(image("figures/montecarlo_14400_NFI.png"), "mcclt", [NFI, 2x2 binning, 14400 frames]),
          subfigure(image("figures/montecarlo_28800_WFI.png"), "mcclt", [WFI, 4x4 binning, 28800 frames]),
      ),
        caption: [Monte Carlo simulation of actual distribution of $Z$ vs CLT approximation]
    ) <clt>

    == Berry-Esseen Bound

    #rt([Section incomplete.  Derivation of $ex(|Y^3|)$ is very complicated but necessary for Berry-Essen bound which gives theoretically guarantees about the convergence of the blue→red lines in the plots above.  May remove this section])

    #math.equation(numbering:none,
        $ex(Y^2) = var(Y) + ex(Y)^2$
    )

  - #rt("WIP: CLT verification experiments")
  - #rt("WIP: Berry-Esseen Theorem")

    #image("figures/clt_wfi.png", height: 15em)
    #image("figures/clt_nfi.png", height: 15em)

    === Calibration and Subtraction <calibration>


  - #rt("WIP")



= Appendix - Boundary Crossing Direction <appendix_crossdir>

  This section describes how to compute the boundary crossing direction across the three types of grid boundaries as described in @section_crossdir.  All vectors are in Cartesian coordinates and we assume the principal axis of the spherical grid to be aligned with the Cartesian $+Z$ axis.


  #figure(
      // grid(
      //     columns: 3, align: bottom,
      // subfigure(image("figures/appendix_crossdir_r.jpg"), "cro", "radial"),
      // subfigure(image("figures/appendix_crossdir_e.jpg"), "cro", "elevational"),
      // subfigure(image("figures/appendix_crossdir_a.jpg"), "cro", "azimuthal"),
      // ),
      image("figures/tomosphero/appendix_crossdir.svg"),
      caption: [Three types of boundary crossings]
  ) <app_crossdir>

  Let $vc(p)$ be the boundary crossing point and $vc(n)$ be the direction of the ray, as illustrated in @app_crossdir.
  To compute the boundary crossing direction over a radial boundary, check for the condition

  #math.equation(
      $#text[negative] = (vc(p) dot vc(n)) > 0$
  )

  For an elevational boundary, check the condition

  #math.equation(
      $vc(m) = vc(p) times (-p_y, p_x, 0) \
          #text[negative] = (vc(m) dot vc(n)) > 0$
  )

  where $p_x$ and $p_y$ are the $x$ and $z$ components of vector $vc(p)$.

  For an azimuthal boundary, check the condition

  #math.equation(
      $vc(m) = vc(p) times vc(n) \
          #text[negative] = m_z < 0$
  )

  where $n_z$ is the $z$ component of vector $vc(n)$.

  == Moments of Random Sum of Random Variables <appendix_rvsum>

    For brevity, denote $E = E_"phot" + E_"mcp"$ and $e = e_"phot" + e_"mcp"$ where $E tilde.op "Pois"(e)$.  Also let $G = G_"mcp"$, first raw moment $mu = ex(G)$, second raw moment $nu = ex(G^2)$, third raw moment $xi = ex(G^3)$ and

    #math.equation(numbering:none,
        $sum_(l=1)^(E) G_l = sum_(l=1)^(E_"phot" + E_"mcp") G_("mcp",l)$
    )

    $G_1, G_2, ...$ are assumed to be i.i.d. random variables.

  Then, by definition of expectation

  #math.equation(numbering:none,
      $ex(sum_(l=1)^(E) G_l)
      &= sum_(L=1)^infinity P(E = L) sum_(l=1)^L ex(G_l) \
      &= sum_(L=1)^infinity P(E = L) L ex(G) \
      &= ex(E) ex(G) = mu e
  $)

  #math.equation(numbering:none,
      $ex(( sum_(l=1)^(E) G_l )^2)
      &= sum_(L=1)^infinity P(E = L)  ex( ( sum_(l=1)^L G_l )^2) \
      &= sum_(L=1)^infinity P(E = L)  ( sum_(l=1)^L ex(G^2) + sum_(l=1)^(L^2 - L) ex(G)^2 )^2 \
      &= sum_(L=1)^infinity P(E = L)  ( L ex(G^2) + (L^2 - L) ex(G)^2 )^2 \
      &= ex(G^2) sum_(L=1)^infinity P(E = L) L + ex(G)^2 sum_(L=1)^infinity P(E = L) L^2 - L  \
      &= ex(G^2) ex(E) + ex(G)^2 (ex(E^2) - ex(E))  \
          & #gt([Substituting moments of $G$ and $E$] + " " + [@poissonmoments]) \
          &= nu e + mu^2 ((e + e^2) - e) = nu e + mu^2 e^2 \
  var(sum_(l=1)^(E) G_l) &= ex(( sum_(l=1)^(E) G_l )^2) - ex(sum_(l=1)^(E) G_l)^2 \
      &= nu e + mu^2 e^2 - mu^2 e^2 = nu e \
  $)

  #rt([FIXME: remove this if not using berry-esseen bound])

  #math.equation(numbering:none,
      $ex(( sum_(l=1)^(E) G_l )^3)
          &= sum_(L=1)^infinity P(E = L) ex( (sum_(l=1)^L G_l)^3 ) \
          &= sum_(L=1)^infinity P(E = L) (sum_(l=1)^L ex(G_l^3) + sum_(l=1)^(6L) ex(G^2) ex(G) + sum_(l=1)^(L^3 - L - 6L) ex(G)^3) \
          &= sum_(L=1)^infinity P(E = L) (L ex(G^3) + 6L ex(G^2) ex(G) + (L^3 - L - 6L) ex(G)^3) \
          &= ex(G^3) sum_(L=1)^infinity P(E = L) L + ex(G^2) ex(G) sum_(L=1)^infinity P(E=L) 6L +  ex(G)^3 sum_(L=1)^infinity P(E = L) (L^3 - L - 6L)\
          &= ex(G^3) ex(E) + ex(G^2) ex(G) 6 ex(E) +  ex(G)^3 (ex(E^3) - 7ex(E))\
          & #gt([Substituting moments of $G$ and $E$] + " " + [@poissonmoments]) \
          &= xi e + 6 nu mu e + mu^3 ((e + 3e^2 + e^3) - 7e)\
          &= xi e + 6 nu mu e + mu^3 (-6e + 3e^2 + e^3)\
  $)

= Appendix - Discretization Considerations <discretization>

  #rt([
      WIP

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
  ])

  #figure(
      image("figures/discretization_considerations.jpg"),
      caption: "Potential discretization Issues"
  )

== Glossary <symbols>



  // - *#acr("LOS")* - #acrfull("LOS")

  #print-glossary(glossary)

  - *MCP* - microchannel plate 23, 24
  - *SNR* - signal-to-noise ratio 26, 27
  - *UV* - ultraviolet: description goes here 22, 23


= Figures


#bibliography("refs.bib")