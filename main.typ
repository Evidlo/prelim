#import "template.typ": *
#import "./diagbox.typ": *
#import "@preview/glossarium:0.5.1": make-glossary, register-glossary, print-glossary, gls, glspl
#import "glossary.typ": glossary
#show: make-glossary
#register-glossary(glossary)

#show link: set text(style: "italic")

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


// ----------- Document Styling ------------

#show heading: set block(above: 1.4em, below: 1em)
#set text(bottom-edge: "descender")
#set grid(row-gutter: 0.5em,)

#outline(indent: auto, depth: 2)

= Introduction

  #rt([
  - Earth's Exosphere and Hydrogen
      - Definition, why do we care about exosphere (Baliukin pg 2)
          - H exosphere indicates presence of water (important for planetary/exoplanetary exploration)
          - Unknown exosphere Lyman-α contributions pollute astronomical UV studies
      - #link(label("radial_profile"))[(figure) Example H density profiles (radius vs density)]
      - Source of hydrogen (Baliukin pg 1)
          - Produced in lower atmosphere, dissociation of h2o & methane, upwards diffusion. escaping and non-escaping atoms
          - Ballistic, satellite, escaping
      - Responds to radiation pressure from sun (Baliukin pg 2 → Thomas & Bohlin)
          - Photon absorption provides antisolar momentum, reemission is anisotropic → net antisolar momentum
      - Time varying response to storm
  - Carruthers Mission Overview
      - Why has exosphere not already been studied?
          - Past measurements have been limited by either low cadence or poor vantage
      - Advantages of L1: good vantage, low fuel maintenance cost, solar power
      - Direct sensing with mass spec. impractical (too large)
  - Prior Exospheric Measurements and Retrievals
      - #link(label("previous_mission"))[(figure) Relative positions of missions]
      // - show results here - refer to historical techniques in @static_retrieval
      - OGO-5 - 1968, 24 Re apogee
          - Early studies of Hydrogen geocorona, detectors saturated
          - Detected H out to 7 Re
          - Later proved existence of Lyman-α background due to interstellar hydrogen flowing through solar system (Baliukin pg 2)
          - Showed existence of geotail (Thomas and Bohlin)
      - Apollo 16 camera - 1872, 63 Re apogee
      - Galileo
          - Took 1/2 day from perigee (1990-12-08) to get to Moon (62Re)
          - Picture taken ~1990-12-13 → 5.5 days → 5.5x2x62 → 680 Re
      - Hisaki 1D (exceed instrument) @hisaki
          - 1000km apogee
          - Discovered exosphere extends out to 38 Re
      - SWAN/SOHO (1996)
          - Discovered geocorona extends to 100Re, beyond moon @baliukin
      - PROCYON/LAICA - Lyman-alpha imaging camera (2015)
      - Dynamics explorer 1 (Rairden, 1986) - spherically symmetric fit
          - Based on numerical solution to RT
      - IMAGE (2000)
          - Geocoronal imager (GEO) instrument
          - Østgaard retrieval @ostgaard
      - TWINS

  - Contributions of This Thesis
      - Spherical raytracer that can be adapted to other missions studying planetary atmospheres
      - Statistical noise model for UV instruments used by Carruthers
      - Reconstruction algorithms built for the Carruthers mission
  - Thesis Organization
  ])

  #figure(
    grid(columns: 2, column-gutter: 1pt,
        subfigure(image("figures/procyon.png", height: 10em), "meas", "PROCYON/LAICA"),
        subfigure(image("figures/mooncarruthers.png", height: 10em), "meas", "Apollo 16 Carruthers camera")
    ),
      caption: "Previous measurements of exospheric Lyman-α"
  )
  #figure(
      image("figures/scratch_radial_profile.png", height: 10em),
      caption: "Radial profile under quiet and storm conditions.  a) Subsolar point b) Geographic North pole"
  ) <radial_profile>


  == Earth's Exosphere and Hydrogen <earth_exosphere>

  == Carruthers Mission Overview

  == Prior Exospheric Measurements and Retrievals

  #figure(
      image("figures/previous_missions.svg", width: 100%),
      caption: "Past observations of exospheric Hydrogen at Lyman-α.\n* not representative of actual spacecraft location"
  ) <previous_mission>

  == Contributions of This Thesis

  == Thesis Organization


= Measurement Constraints <measurement_constraints>

  #rt([
  - [x] Chapter introduction section
      - What is Lyman-α - Spectral line of hydrogen
          - Only practical indicator of hydrogen distribution
      - Carruthers carries two UV-capable cameras to provide constraints on hydrogen distribution
      - Chapter summary: Carruthers orbit, camera specs, emission model, post-processing
  - [x] Carruthers Orbit and Camera Details
      - Definition of #gls("GSE") coordinate system
      - Carruthers will be inserted in a halo orbit around Lagrange point L1 1.5e6 km from the Earth
          - From distant vantage, should be able to observe full extent of geocorona (soho/swan)
          - Orbit duration information and lateral deviation (important for tomographic meas. diversity)
      - #link(label("viewgeom1"))[(figure) view geometry, orbit around L1]
      - #link(label("viewgeom2"))[(figure) static 3D picture of multiple view geometries captured for 1month baseline]
      - [x] Camera details
          - NFI: 3.6°, 30 min for fast inner exosphere evolution
          - WFI: 18°, 60 minutes for global estimation of slowly-changing outer exosphere
          - #link(label("viewgeom3"))[(figure) camera details, FOV, int. time, etc.]
          - #link(label("camera_specs"))[(table) camera geometry specifications]
          - #link(label("carruthers_orbit"))[(figure) LOS evolution over 15 days]
  - Emission Model
      - [x] Definition of emission model - necessary for retrieval
      - [x] Definition of optically thin and importance of this assumption on emission model computational complexity (Zoennchen 2015)
      - [x] A correction term applied by Gonzalo, but left out of this manuscript
          - Note that raytracer is capable of implementing these correction terms.  See appendix XXX for an overview of implementing correction terms using raytracer
      - #link(label("emission_model_physics"))[(figure) Physics of emission model]
      - g-factor - FIXME: varies with radial distance but assumed constant, see @static_validation for bias
      - Definition of albedo (Gonzalo thesis appendix c)
      - Background sources
          - IPH, moon, stars
          - Diagram showing hydrogen density in solar system?
          - One component of background radiance is IPH
              - tenuous distribution of hydrogen between solar-system planets are illuminated and contribute to radiance measurement in a LOS-dependent way
              - estimating IPH is involved, carruthers will measure in annulus around earth where exospheric hydrogen radiance is nearly zero and interpolate over the exosphere FOV
              - this manuscript assumes IPH is known, bias introduced from IPH estimation process are considered in @static_validation
  - Instrument, Calibration
      - #rt([less detailed version of @calibration])
      - equation for calibrated column density

      - Ability to simulate realistic images is critical to validate algorithm performance under noise conditions
      - Previous section described model for interaction between photons and hydrogen
      - This section will describe model for interaction of photons and instrument
      - Brief description of MCP, KBr photocathode, etc. (use overview list from Instrument Model section)
      - See @appendix_sim for simulation and calibration details
          - Science pixel binning - computational simplicity
          - Masking - optically thick region, moon, stars
  - Post-Processing
  ])

  #rt([FIXME: Lara has better version?]) #image("figures/viewgeom1.svg", height:10em) <viewgeom1>
  #rt([delete]) #image("figures/viewgeom2.gif", height:10em) <viewgeom2>
  #rt([delete]) #image("figures/viewgeom3.png", height:10em) <viewgeom3>
  // #image("figures/emissionmod.png", height:10em) <emissionmod>

  The Sun is a strong source of Lyman-α photons, which enter the Earth's atmosphere and interact with neutral hydrogen atoms before being reemitted in a process known as resonant scattering.  Aside from in-situ spectrometric measurements of hydrogen (impractical for global atmospheric measurements), this Lyman-α emission is the only available indicator of hydrogen in the exosphere.
  Carruthers is equipped with UV-capable cameras and will observe the entirety of the exosphere at Lyman-α from a distant vantage, providing constraints on the global distribution of hydrogen.

  This chapter will describe the Carruthers observational orbit, camera geometry, an emission model of Lyman-α interaction with and propagation through exospheric hydrogen, instrument model, and end with an overview of the post-processing procedure for converting raw sensor measurements to usable science constraints.


  == Carruthers Orbit and Camera Geometry

    The #gls("GSE") coordinate system (shown in @gse_coordinates) is a natural fit when describing spacecraft position and exospheric hydrogen distribution, since properties of the exosphere such as the hypothesized nightside (#rt("FIXME: first time this is mentioned in thesis")) tail are aligned to an Earth-Sun frame.  The rest of this manuscript uses #gls("GSE") in units of Earth radii (1 Re = 6371 km) either in cartesian or spherical coordinates.

    #figure(
        image("figures/gse_coordinates_placeholder.jpg", height: 10em),
        caption: [#rt([FIXME: placeholder.]) GSE coordinate system.  +X is Sunwards, +Y is Earth wake (approximately dusk), +Z is ecliptic North]

    ) <gse_coordinates>

    Carruthers will be inserted into a halo orbit around the L1 Lagrange point (about 1.5 million km from Earth, 235 Re), which is distant enough to observe the entirety of the geocorona from an outside vantage @baliukin over a long period.

    Shown in @carruthers_orbit, this orbit deviates above/below the ecliptic plane by ±7° (xxx Re #rt([FIXME])) and ±28° (xxx Re #rt([FIXME])) in the dawn/dusk direction, providing important angular measurement diversity for tomographic analysis.
    While multiple spacecraft would be ideal for improving spatiotemporal measurement resolution from the relatively slow 6 month orbital period, analysis in @static_validation shows that an observation window of just 2 weeks is sufficient to meet mission retrieval requirements during quiet exospheric conditions.

    #figure(
        grid(
            columns: 2, gutter: 5pt,
            subfigure(image("figures/carruthers_orbit_placeholder.jpg", width: 15em), "orbitfig", "6 month Carruthers orbit around L1"),
            subfigure(image("figures/los_evolution_15d.png"), "orbitfig", "Angular diversity of orbit in ecliptic plane")
        ),
        caption: [#rt([FIXME: placeholder.]) Carruthers orbit details]
    ) <carruthers_orbit>


    Two cameras within Carruthers make up the #gls("GCI"), the primary instrument on board the spacecraft.  The #gls("NFI") observes the inner exosphere at high resolution on a 30 minute cadence and the #gls("WFI") observes global distribution of exospheric hydrogen out to 25 Re on a 60 minute cadence.  Together, these cameras are capable of measuring both the fast-evolving dynamics of the inner exosphere and global distribution of hydrogen along with interplanetary background Lyman-α signal (discussed in the next section).  @earth_fov shows the #gls("FOV") of the cameras relative to the inner exosphere boundary and outer exospheric limit (#rt[justifying 25 Re here?  exosphere goes out to 100 Re and beyond]).  A summary of camera geometry specifications is given in @camera_specs.

    #figure(
        image("figures/scratch_fov.jpg"),
        caption: [#rt([FIXME: placeholder.])  (a) Carruthers camera FOV projected onto Earth tangent plane. (b) FOV relative to outer exosphere boundaries.]
    ) <earth_fov>

    #figure(
        table(
            table.header([Camera], [FOV\ (degrees)], [Resolution\ (pixels)], [Angular Res.\ (degrees)], [Spatial Res.\ (Re, projected)]),
            align: horizon,
            [WFI], [18°], [512²], $35·10^(-3)$, [0.14],
            [NFI], [3.6°], [1024²], $3.5·10^(-3)$, [0.014],
            columns: (auto, auto, auto, auto, auto),
        ),
        caption: [Camera geometry specifications.\ Spatial resolution is projected onto Earth tangent plane as in @earth_fov(a)]
    ) <camera_specs>


  == Emission Model <emission_model>

    In order to recover a hydrogen density distribution from measurements, it is necessary to mathematically model the process by which Lyman-α photons propagate through the exosphere and enter the camera.
    This is known as an emission model and is a central component of tomographic retrieval algorithms.

    Numerically modelling the physics of radiative transfer is a computationally complex task, as photons entering the atmosphere usually scatter multiple times in several locations, creating complicated interdependencies between distant portions of the exosphere.  However, in regions of the atmosphere where hydrogen is sparse (known as the #gls("optically thin") regime #rt([FIXME: citation])), it is possible to assume photons scatter only once without significant loss of accuracy, simplifying computational requirements and implementation complexity of the emission model.

    Anderson and Hord Jr. @opticaldepththin, define the optically thin regime as starting when $tau <= 0.1$, where optical depth $tau$ is a measure of #rt([FIXME: optical depth description]).

    With these assumptions, the measurement by the spacecraft is proportional to a line integral of the total mass of hydrogen present in the atmosphere along the #gls("LOS") of a particular pixel, given by photon radiance

    // #math.equation(
    //     $y' = g phi.alt(beta_(i,j)) integral_(l=0)^(infinity) bold(a)(vc(x) + vc(n)l) bold(rho)(vc(x) + vc(n)l) dif l$
    // ) <integral2>
    #math.equation(
        $I_"exo,t" = g^*_t phi.alt(beta) integral_l bold(a)_t (vc(r))  bold(rho)_t (vc(r)) dif s #gt([(phot/s/cm²/sr)])$
    ) <integral2>

    #rt([FIXME: notation inconsistent with other sections])


    where $l$ is the pixel #gls("LOS"), illustrated in @emission_model_physics.
    Through Beer's law and a logarithmic transform, this integral equation is also valid for many modalities featuring absorptive rather than emissive media. #rt([FIXME: citation])


    // knlown as a #gls("LOS"), an example of the Fredholm integral of the first kind.


    // #math.equation(
    //     $y prop integral_V bold(rho)(v) dif v$
    // )


    where $g^*_t$ is angular g-factor, $t$ is time, $phi.alt_t$ is scattering phase function, $bold(a)_t$ is albedo, $vc(r)$ is position vector, $bold(rho)$ is hydrogen density and $V$ is the exosphere volume in the #gls("FOV").  $phi.alt$, known as the #gls("scattering phase function"), represents the directional distribution of resonantly-scattered photons relative to the direction of the sun, shown in @scatteringphase.

    #figure(
        image("figures/physics_scattering.svg", height: 15em),
        caption: [Scattering phase function]
    ) <scatteringphase>

    Angular solar g-factor $g^*_t$ (phot/s/atom/sr) relates hydrogen quantity (atoms) to photon flux emitted per unit solid angle (phot/s/sr).  This quantity is directly related to solar activity and is assumed to be a known, external input to the emission model.
    @gfactor_difference illustrates the important $4 pi$ difference between angular g-factor (per steradian) and _isotropic_ g-factor (per whole sphere), which are sometimes not explicitly distinguished in literature.

    #figure(
        grid(
            columns: 2, gutter: 12pt,
            subfigure(image("figures/g_factor_angular.svg", height:8em), "gfact", "Angular"),
            subfigure(image("figures/g_factor.svg", height:8em), "gfact", "Isotropic"),
        ),
        caption: [Angular vs. Isotropic g-factor],
    ) <gfactor_difference>

    Finally, $bold(a)_t (vc(r))$ is a unitless multiplicative correction factor (assumed known) to account for high-density regions of the inner optically thick exosphere which act as a secondary source of Lyman-α photons illuminating the outer exosphere from below.

    - #rt([FIXME: derivation of line integral from first principles?  Make use of etendue, pixel area, paraxial approximation to convert volume integral to line integral])


    #figure(
        image("figures/physics.svg", height:20em),
        caption: [Emission model overview for a #gls("LOS")]
    ) <emission_model_physics>

    @emission_units gives an overview of different measurement quantities and their units.

    #figure(
        table(
            columns: 3,
            align: horizon,
            table.header([Quantity], [Symbol], [Units]),
            "Density", $bold(rho)$, "atom/cm³",
            "Non-spectral radiance", $I_"exo"$, "phot/s/atom/cm²/sr",
            "Isotropic g-factor", $g^*$, "phot/s/atom/sr",
            "Albedo", $bold(a)$, "unitless",
            [Scattering phase function], $phi.alt$, "unitless"
            // "Emissivity", "", "",
            // "Pixel solid angle", "Ω", "steradian (sr)",
            // "Pixel area", "A", "cm²",
        )
    ) <emission_units>

    Zoennchen et al. #rt([(FIXME: cite Gonzalo)]) have found that extending the basic optically thin assumptions to additionally model extinction along the LOS #rt([(FIXME: and the other term)]) reduced discrepancy between tomographic retrievals and physics-based simulations.  The scope of this manuscript does not include these correction terms, but @appendix_extra_physics describes the procedure for implementing these terms with the raytracer in @raytracer.

    A tenuous distribution of hydrogen throughout the solar system, known as #gls("IPH"), also contributes to the radiance detected by the spacecraft.  Estimation of #gls("IPH") is an involved process, and Carruthers will dedicate a portion of on-orbit operations to making observations of an annulus around the Earth where exospheric hydrogen is not present.  @iph shows a typical #gls("IPH") distribution expected to be observed during the mission.  #gls("IPH") radiance contribution from behind the exosphere envelope is interpolated from the measured annulus.  The impact of biases introduced by #gls("IPH") estimation are considered in @static_validation.
    Other unwanted sources of Lyman-α signal which violate emission model assumptions include the moon and stars and optically thick exosphere, as shown in @moon_stars, but these sources will be masked out and ignored during retrieval instead of estimated.  Together, these radiance sources are referred to as $I_"bkg"$.


    #figure(
        image("figures/iph_placeholder.jpg", height: 10em),
        caption: [#rt([FIXME: placeholder]) IPH distribution and observation annulus]
    ) <iph>

    // The quantity in @integral1 is known as _column density_, but on orbit the spacecraft is actually measuring _photon flux_ shown in @integral2, the total photon rate entering the instrument.

  == Instrument Model <instrument_model>

    #rt([
        // FIXME:

    // While the volume integral in @integral1 is accurate for any pixel geometry, the expression can be simplified to a line integral when the pixel #gls("FOV") is small.  Using conservation of etendue and the paraxial approximation (pixel #gls("FOV") lies close to the optical axis), @integral1 can be written

    // #math.equation(
    //     $I_"exo,t" = g^*_t phi.alt(beta) integral_V bold(a)_t (vc(r)) bold(rho)_t (vc(r)) dif v #gt([(phot/s/cm²/sr)])$
    // ) <integral1>

    ])

    #rt([FIXME: question for Lara: past section uses standard convention $I_"exo"$, but this is incompatible with convention of capitalized letters being random variables.  Would $i_"exo"$ be OK?])

    As mentioned in the previous section, accurately modelling the physical processes involved in obtaining tomographic measurements of a density distribution is critical for retrieval accuracy.  While the emission model describes the physics of Lyman-α photons interacting with the exosphere, an instrument model explains how photon radiance received at the front of the instrument passes through the stages of the #gls("GCI") and is converted to digital measurements.

  This section describes a statistical model for the instrument noise and background signals present in the NFI and WFI cameras during measurement of exospheric Lyman-α.  Modelling these processes is important for converting raw sensor measurements in digital numbers (DN) as telemetered by the spacecraft to corresponding radiances that can be used for tomographic reconstruction.   A statistical model is also important for generating synthetic noisy measurements to validate the performance of retrieval algorithms.  As a result, the Carruthers cameras have undergone extensive laboratory characterization to determine instrument model parameters and periodic on-orbit calibration is planned to account for parameter drift due to exposure to the space environment.

    The Carruthers spacecraft contains two #gls("UV") cameras designed for detecting exospheric Lyman-α emission, an indicator of the amount of atomic hydrogen present along a #gls("LOS").  These cameras utilize a design which has heritage with other #gls("UV") instruments such as ICON @icon and GUVI @guvi and contain the following stages @rider (shown in @instrument_stages):

  // TODO: WIP

  - *optical filter wheel* - provides multispectral capability for out-of-band signal measurement used in thick exosphere retrievals
  - *KBr photocathode* - potassium bromide photocathode for converting #gls("UV") photons to electrons (known as _photoelectrons_)
  - *#gls("MCP")* - amplifies individual photoelectrons to a detectable shower of electrons
  - *phosphor and #gls("CCD")* - phosphor produces light pulses from electrons which are detected by a #gls("CCD") and converted to an electrical charge
  - *ADC* - #gls("ADC") for reading out #gls("CCD") charge.  (together with the #gls("CCD") this is sometimes referred to as an #gls("APS"))

  #figure(
      box(width:100pt, height:100pt, stroke:1pt),
      caption: [#rt([FIXME: obtain figure from Lara]) Optical stages of instrument.]
  ) <instrument_stages>

    The rest of this section will consist of a derivation of a single pixel noisy measurement in #gls("DN") given a photon spectral radiance and other quantities in @knownvariables and @randomvariables.

  Let $L_("exo")(lambda) + L_("bkg")(lambda)$ represent exospheric photon spectral radiance and background.
    After entering the instrument and passing through several optical stages, photon spectral radiance is converted to a photon flux
  #math.equation(numbering: none,
      $p_("phot")(lambda) = lr(\[D(L_("exo")(lambda) + L_("bkg")(lambda)) * h)\]_j dot.op A dot.op Omega gt("(phot/s/nm)")$
  )
    where $A dot.op Omega$ is pixel etendue, $D$ is a non-linear spatial distortion, and $h$ is a convolutional kernel representing optical blur.  Subscript $j$ denotes that only a single pixel is considered of the spatially-distorted and blurred signal.
    Note that the above expression implicitly assumes that radiance is constant over the entire pixel, an approximation that holds true so long as hydrogen density does not vary much in the #gls("FOV") in the direction perpendicular to the #gls("LOS").
  After entering the instrument, the photons pass through an optical filter and photocathode where they are converted to a stream of photoelectrons (also known as _events_) with rate $e_"phot"$.  This conversion happens with efficiency $f_"opt" epsilon(lambda)$ (events/phot) where unitless scaling factor $f_"opt"$ (known as a _flat-field_) accounts for spatial variation in the optics, filters and photocathode.

  #math.equation(numbering: none,
      $e_"phot" = f_"opt" integral_lambda p_("phot")(lambda) epsilon(lambda) dif lambda gt("(events/s)")$
  )

  Alternatively, if $L_"exo"$ and $L_"bkg"$ are monochromatic with intensity $I_"exo"$ and $I_"bkg"$ and wavelength $lambda_0$

  #math.equation(numbering:none,
      $e_"phot" &= f_"opt" integral_lambda p_("phot")(lambda) epsilon(lambda) dif lambda \
          &= f_"opt" lr(\[D(integral_lambda L_("exo")(lambda) epsilon(lambda) dif lambda + integral_lambda L_("bkg")(lambda)  epsilon(lambda) dif lambda) * h\])_j dot.op A dot.op Omega \
          &= f_"opt" epsilon(lambda_0) \[D(I_"exo" + I_"bkg") * h\]_j dot.op A dot.op Omega gt("(events/s)")$

  )

  Assuming Lyman-α photons are emitted in a Poisson process, the actual number of photoelectrons created in a single camera frame of duration $t_"fr"$ is given by random variable

  #math.equation(numbering:none,
      $E_"phot" tilde.op "Pois"(t_"fr" e_"phot") gt("(events)")$
  )

  In general, single photoelectrons are difficult to detect, so the Carruthers cameras employ an #gls("MCP") for turning a single particle into a detectable shower of particles.  #glspl("MCP") consist of an array of small glass tubes (channels) which are electrically charged so that a photoelectron striking the wall of one of these tubes will cause a cascade of particles via secondary emission @microchannelplate.  The small size of these channels ensures that the subsequent shower of particles exits the #gls("MCP") in the same location as the photoelectron, preserving image spatial resolution.  Due to the nature of secondary emission, the number of particles created by the #gls("MCP") from a photoelectron is given by the discrete random variable $G_"mcp"$, whose distribution has been measured in the laboratory.  The number of electrons leaving the MCP #rt([(known as _counts_)]) due to Lyman-α photoelectrons is given by

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
      image("figures/bd_annotated2.svg", width: 120%),
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

    == Frame Stacking, Time Binning, and Instrument Model Approximation

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

    @framestacking shows examples of simulated measurements which utilize both frame stacking and pixel binning and compares them to a single, unbinned frame.  More features of the underlying radiance distribution become visible when many frames are stacked.
    @integrationtime gives an overview of the stacking and binning parameters selected in order to meet mission #gls("SNR") requirements and maintain an acceptable loss of spatial and temporal resolution given prior knowledge of the exosphere @gonzalostorm.

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

      Naively stacking and binning thousands of noisy frames generated using the model described in the last section can be computationally, especially when generating stacked images from many vantages.  Instead, the #gls("CLT") can be applied to generate a noisy stacked image $Z$ efficiently, assuming all individual frames are identically distributed.  #gls("CLT") states

    #math.equation(numbering:none,
        $Z tilde.op cal(N)(M N^2mu_Y, M N^2sigma_Y^2)$
    )

      where $M$ is the number of stacked frames, $N$ is the number of binned pixels, and $u_Y$ and $sigma^2_Y$ are the mean and variance of individual frames.  See @instrument_clt for a derivation of these values and numerical validation of the #gls("CLT") approximation.

  == Post-Processing <post_processing>

    As will be shown in @inverse_problem, many retrieval algorithms amount to applying the emission and instrument models to a candidate density and comparing the resulting candidate measurements to real measurements obtained on orbit to update the density in an iterative fashion.  However, this can be computationally expensive, especially when a retrieval algorithm requires hundreds or thousands of iterations to converge to a solution.   An alternative to applying the emission and instrument models every iteration is to reverse the effects of these models on the real measurements in a process known as _subtraction_ or _calibration_.  Subtraction need only occur once after images are downlinked, greatly accelerating the iterative retrieval process as in @subtraction_efficiency.

    #figure(
        image("figures/subtraction_efficiency_placeholder.jpg", height:25em),
        caption: [Subtraction eliminates the need for most of emission and instrument model from the retrieval loop.],
    ) <subtraction_efficiency>

    // Subtraction greatly accelerates retrieval and involves reversing the effects of all steps described in the past two sections excluding the integral given in @integral2.
    While the details of subtraction are out of scope of this manuscript, an overview of the steps is below.  Refer to the expression for $mu_Y$ given in @clt_mean and $I_"exo"$ given in @integral2.

    #figure(
        align(left, box(
            stroke: 1pt,
            inset: 1em,
            [
                #gt([Instrument model subtraction])
                1. Divide out effects of stacking by $M$ frames and binning by $N$ pixels
                2. Subtract off read noise bias
                3. Divide out #gls("APS") gain and camera integration time
                4. Subtract external signal sources like dark current, #gls("APS") noise
                5. Divide out #gls("MCP") gain and #gls("MCP") flat-field
                6. Subtract MCP noise contribution
                7. Divide out pixel etendue, optical flat-field, and efficiency
                8. De-distort and de-blur total radiance
                9. Subtract background radiance, including estimated #gls("IPH")
                #gt([Emission model subtraction])
                10. Divide out g-factor and scattering phase function
            ])),
        caption: [Subtraction process overview]
    )

    Notably, subtraction excludes the integral given in @integral2, as this is the job of the retrieval algorithm.
    The resulting estimate is known as a _column density_ and is independent of any instrument parameters or solar conditions

    #math.equation(
        $y_(t,j) approx integral_(l_(t,j)) bold(a)_t (vc(r)) bold(rho)_t (vc(r)) dif l gt("(atom/cm²)")$
    )

    where #gls("LOS") $l$ is annotated with index $t,j$ to emphasize time and pixel dependence.

    As mentioned previously, some #gls("LOS") contain Lyman-α signals which are unknown or violate the emission model assumptions and must be marked so they are ignored during retrieval.  These include the moon, stars, optically thick exosphere, and Earth shadow, shown in @moon_stars. (#rt([FIXME: cite old Zoennchen]))

    - #rt([FIXME: refer to Earth shadow mask in albedo intro section])

    #figure(
        grid(
            columns: 3, gutter: 10pt,
            subfigure(image("figures/moon_mask_placeholder.jpg"), "moonstars", "Moon in WFI FOV"),
            subfigure(image("figures/star_mask_placeholder.jpg"), "moonstars", "Stars in WFI FOV"),
            subfigure(image("figures/thick_mask_placeholder.jpg"), "moonstars", "Optically thick exosphere and shadow"),
        ),
        caption: [#rt([FIXME: placeholder]) IPH sources that must be masked. (#rt([combine these into a single figure that illustrates and labels all 3?]))]
    ) <moon_stars>

    #figure(
        image("figures/scratch_fov.jpg"),
        caption: [#rt([FIXME: placeholder.])  (a) Polar-binned camera FOV projected onto Earth tangent plane. (b) Polar-binned FOV relative to outer exosphere boundaries.]
    ) <earth_fov_science>

    A final step of post-processing is to bin the 512² pixel WFI and 1024² pixel NFI images to a lower resolution.  The dense #gls("LOS") spatial sampling and fast temporal cadence of these cameras creates significant memory and computational requirements for retrieval algorithm, which can be avoided by downsampling images in a way that preserves as much detail as possible.  This is similar to the stacking and binning described in @instrument_model, but occurs after downlink.  A good binning scheme is a log-polar grid, due to the roughly spherically-symmetric, radially exponential distribution of the exosphere.  Within the Carruthers mission, this is referred to as _science pixel binning_ and shown in @earth_fov_science. Similarly, temporal averaging can reduce the volume of data ingested by the retrieval algorithms without much loss of information on exosphere dynamic evolution, especially during quiet conditions.  These data reduction techniques help improve the memory resources required during retrieval, which is analyzed in @benchmarking.

    This chapter has focused on treating a single detector pixel $j$ and associated #gls("LOS") at a single time index $t$ for the sake of notational simplicity, but tomographic retrieval algorithms rely on many measurements from different angles to fully constrain a reconstruction.
    This manuscript uses the convention $bold(y) = {y_(t,j) | forall t, forall j}$ to denote a stack of subtracted 3D measurements, which will be used extensively in the next section.


= Inverse Problem Formulation <inverse_problem>

  #rt([
      - Purpose of tomography is retrieval of a volumetric structure which produced a set of measurements
          - This process is generally known as an inverse problem
      - Discretization
          - Discretization necessary for numerical approximations to be made when analytic solution is infeasible
          - Also rediscretize measurements to reduce computational burden (refer to section XXX. FIXME)
          - Spherical grid definition in GSE coordinates
      - Inverse problem
          - Tomography can be formulated as linear inverse problem y = Lx
          - Usually _underdetermined system_. Direction inversion not possible
          - Solution often formulated as minimizing some objective function or _cost function_ such as $hat(rho) = min_rho ||bold(y) - L rho||_2^2$ (least squares)
          - Generalized least squares can find a solution to underdetermined and overdetermined systems (pseudoinverse)
          - Alternatively, formulate solution space as low-dimensional subspace or manifold
          - Solution may be unstable under noise due to smoothing effects in common inverse problems
          - Regularization impacts specific method for solution
              - e.g. Tikhonov - $||bold(y) - L bold(rho) ||_2^2 + ||T x||_2^2$
                  - hash explicit closed form solution (Regularized least squares)
          - Hadamard defined these concepts as _ill-conditioned_ and _ill-posed_
          // - Under conditions described in @measurement_constraints (namely, single scattering in optically thin regime), tomographic inversion can be formulated as solution to the linear inverse problem $y = L x$ (ignoring noise)
          // - Direct inversion of the tomographic operator requires that the forward matrix L be non-singular in order for an inverse L^-1 to exist
          // - Measurement constraints are limited - less rows than free variables. ill posed
          // - Not injective (not one-to-one) (measurement might correspond to more than one potential solution)
          // - Not surjective (not onto) (a measurement may not necessarily even correspond to a feasible solution, e.g. noise)
          // - Still necessarily enough to get good retrieval under presence of noise
          // - Changes in density do not have significant impact on measurements.
          // - Inversely, small changes in measurements have large effects on solution during retrieval - ill conditioned @illposed
          //     - Potentially problematic in the presence of noise
          // - Ill posedness and conditioning fixed through regularization, low rank models
          // - Concepts defined mathematically by hadamard @hadamard
          // - Generalized inverse (moore-penrose pseudoinverse), direct least squares minimization, SVD, RRPE
          // - Selecting hyperparameter values - trial and error, GCV, many others (gonz thesis pg 49)
          - Solve with iterative method
      - Math notation (model, operator, coefficients, measurements, density), etc.
  ])



  An inverse problem is the procedure of determining the causative factors of a set of measurements derived from an observation process.  In exospheric tomography, the factor driving the intensity of column density measurements is the distribution of hydrogen in the regions being observed.  Direct analytic solutions to tomographic or other inverse problems are not always possible, so numerical approximations and discretization become necessary.  In this chapter, I lay out key concepts of linear inverse problems, detail a discretization scheme for approaching tomographic inversion numerically, and introduce notation which will be used later in the manuscript to describe tomographic retrieval algorithms.

  == Discretization

    #rt("FIXME: convert line integral to discrete sum")

    In general, direct analytic solutions to tomographic estimation problems are infeasible, necessitating a numerical approach where the solution space is divided into a finite grid of $N$ non-overlapping regions called #glspl("voxel") where density is assumed to be constant.  There are a large variety of grid types to choose from, including regular grids (e.g. spherical, cylindrical, cartesian), non-regular grids which may utilize hierarchical structures (e.g. octree) and tetrahedral meshes. // as shown in @grid_examples.
    Some of these discretizations schemes have been designed to adaptively update voxel boundaries during retrieval to better fit the object being retrieved @adaptivemesh1.

    // #figure(
    //     box(width:100pt, height:100pt, stroke:1pt),
    //     caption: [Discretizations of a spherical domain\ #rt([FIXME: WIP])]
    // ) <grid_examples>

  Choosing a discretization grid which is well-suited to the data is critical, as an improper grid can create aliasing and other artifacts that cause the numerical result to deviate from the functions they approximate.  Since the exosphere is well-understood to smoothly vary with larger density gradients at lower altitudes, a regular spherical grid with logarithmically spaced radial bins is appropriate.  The nature of regular grids allows for a simple multidimensional array as the underlying data structure, and the property of shared boundaries between voxels simplifies tomography calculations, as demonstrated in @alg_outline.

    The regular spherical grid used in this manuscript has its 0° elevation pole aligned to ecliptic north (+Z Cartesian GSE axis) and ±180° azimuth branch point pointed away from the sun (-X Cartesian GSE axis), as shown in @grid_details.

    This manuscript uses convention $r$, $e$, $a$ when referring to radial, elevational, and azimuthal dimensions to avoid ambiguity with astrophysical versus mathematical conventions for spherical coordinates as in @sph_convention that denote the elevational and azimuthal dimensions differently.


    #figure(
        image("figures/sph_convention.svg", height: 10em),
        caption: [Spherical grid definition and a spherical voxel]
    ) <grid_details>

    #figure(
        table(
            columns: (100pt, auto, auto, auto, auto),
            align: center + horizon,
            table.header([Dimension], [Physics\ Convention], [Math\ Convention], [Manuscript\ Convention], [Manuscript\ Range]),
            table.hline(stroke: 2pt),
            [Radial], $r$, $r$, $r$, $[3 "Re", 25 "Re"]$,
            [Elevational], $theta$, $phi$, $e$, [$[0°, 180°]$],
            [Azimuthal], $phi$, $theta$, $a$, [$[-180°, 180°]$]
        ),
        caption: [Spherical coordinate conventions.]
    ) <sph_convention>

    The grid resolution necessary to avoid aliasing and sample errors is data-dependent and is covered in @discretization.

  == Inverse Problem

  Under the conditions described in @measurement_constraints (single scattering in optically-thin exosphere) and ignoring noise, tomographic inversion can be formulated as the solution of the linear inverse problem

  #math.equation(
      $bold(y) = F bold(rho)$
  )


  #math.equation(
      $y_"tij" = sum_(r, e, a) F_("ij","trea") rho_"trea"$
  )

    with measurements $bold(y)$, forward operator $F$, solution $bold(rho)$, time t, pixel ij, and spatial voxel rea.

    Direct inversion of the tomographic operator requires that the matrix $F$ be non-singular in order for an inverse $F^(-1)$ to exist.  However, tomography problems generally have fewer measurement constraints (#gls("LOS", display:"lines of sight")) than free variables (#glspl("voxel")), making them _ill-posed_ (lacking existence, uniqueness or stability of solutions).

    While inversion of $F$ is impossible, the problem of obtaining estimate $hat(bold(rho))$ is sometimes reformulated as a minimization of some objective function, such as the common least squares

    #math.equation(
        $hat(bold(rho)) = min_bold(rho) ||bold(y) - F bold(rho)||_2^2$
    ) <cost_function>

    A _generalized inverse_ such as the Moore-Penrose pseudoinverse, can be constructed from $F$ which selects the solution with the smallest norm or which best fits the measurements.  Another approach is to assume that $bold(rho)$ has some low degree-of-freedom representation on a subspace (linear) or manifold (non-linear) within the space of all solutions.  Such a mapping $M$ is referred to as a #gls("model") in this manuscript and is represented as

    #math.equation(
        $bold(rho) = M(bold(c))$
    )

    where $bold(c)$ are the low degree-of-freedom coefficients.  If the model is sufficiently low-dimensional and linear, then uniqueness is guaranteed if a solution exists.

    However, this solution may still be _ill-conditioned_, or unstable in the presence of noise.  In these cases, $F$ is usually rank-deficient or has small singular values meaning that small perturbations in the measurements $bold(y)$ are amplified greatly in the inversion process.
    This is especially true of computational imaging problems like tomography, where smoothing effects of integration dampen the high frequency information about $bold(rho)$ available in the measurements.

    A common technique is to add additional terms to @cost_function which encode prior knowledge about the distribution of $bold(rho)$ which may not be captured by the model.  This is known as _regularization_ and can be used to combat ill-posedness.
    An alternative approach is to view inverse problems as a probabilistic inference problem which uses Bayes theorem to model uncertainty in reconstruction $bold(rho)$ given measurements $bold(y)$.

    The issues of solution existence, uniqueness and sensitivity in inverse problems were formally described by Hadamard in the early 20th century and used to define ill-posed and ill-conditioned inverse problems @hadamard @illposed.

    The notation introduced in this section is summarized in @symbols and will be used throughout the manuscript.

  #figure(
      pad(x: -50%,
          table(
              table.header([Symbol], [Meaning], [Description], [Shape]),
              align: horizon,
              columns: (5em, 10em, 14em, 12em),
              $F$, [Forward operator], [Mapping from 3D H-density to column density], [$F: bb(R)^3 → bb(R)^3$ (static) \ $F: bb(R)^4→bb(R)^3$ (dynamic)],
              $M$, [Model], [Mapping from model \ parameters to 3D H-density], [$M: bb(R)^* → bb(R)^3$ (static) \ $M: bb(R)^* → bb(R)^4$ (dynamic)],
              $bold(c)$, [Model \ params./coeffs.], "Free model variables,\n usually low-dimensional.", [$bold(c) ∈ bb(R)^*$ \ \* model dependent],
              $bold(rho)$, [H density], [Spatial distribution of \ exospheric Hydrogen], [$bold(rho) ∈ bb(R)^3$ (static) \ $bold(rho) ∈ bb(R)^4$ (dynamic)],
              $bold(y)$, [Measurements.], "Column densities measured\n by instrument", [$bold(y) ∈ bb(R)^3$],
          ),
      ),
      caption: [Symbols],
      // FIXME: move to glossary section?
  ) <symbols>



= Raytracing Optically Thin Gases <raytracer>

  // Tomography is a method for determining the internal structure of objects from a set of measurements which penetrate into the object being measured.


// Fast tomographic reconstruction algorithms that implement explicit inversion formulas typically work only for specific view geometries (such as circular or helical path) and are referred to as _filtered back projection_ (FBP) algorithms @fbp.
//   However, spaceborne sensors have view geometries that are determined by orbital parameters of a spacecraft.
//   For these situations requiring more flexible view geometries where an exact inverse solution is not available, _iterative reconstruction_ (IR) algorithms prevail at the expense of increased computational complexity.  Examples include SIRT @sirt, TV-MIN @tvmin, ART @art, CGLS @cgls, Plug-and-play @plugandplay and many others.
//   These algorithms depend on generating synthetic projections of a candidate object using an operator (sometimes called a _raytracer_) that simulates waves traveling through the object medium.
//   They produce a reconstruction by repeatedly tweaking the candidate object to minimize discrepancy between synthetic and actual projections, and they stand to benefit the most from a fast operator implementation.

  In this section, we describe the design and implementation of a tomographic projector which is designed for use in iterative reconstruction algorithms on spherical grids.  The algorithm is automatically differentiable, GPU-enabled, and easily integrable into machine learning methods through PyTorch.
  Additionally, the implementation provides a suite of visualization functions described in @api_overview which can be used for visual verification of view geometry and discretization or for generating publishable figures.
  The raytracer has been released as a open-source Python package *TomoSphero*.@tomosphero.

=== GPU Utilization

  Most tomographic operators treat each pixel on the detector as an independent computational task, which has led to the development of tomography libraries that are capable of simultaneously utilizing multiple cores on a CPU or hardware accelerator.  TomoSphero is parallelized and GPU-enabled, and its speed has been benchmarked as described in @benchmarking.

  In cases where a simultaneous computation for every pixel of every measurement would consume more memory than is available, some algorithms operate _out-of-core_, where they parallelize as many tasks as will fit into available memory, then serially queue the remaining tasks for processing after current tasks are complete.  TomoSphero is not capable of out-of-core operation, so a characterization of its memory usage is also presented in @benchmarking.

=== Object Discretization

  #rt([FIXME: make this mesh with previous section about discretization, or delete])

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

  For the cases where a boundary has no intersection points with a ray, coordinates are still computed and stored for these points, but the associated value in the distance array $l$ is set to #emph("inf"), as illustrated in @alg_intersection.

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
      caption: [Intersection points along a LOS annotated with crossing direction (↑:positive, ↓:negative)]
  ) <alg_crossdir>
  #figure(
      image("figures/tomosphero/alg_crossdir2.svg", height: 10em),
      caption: [Region indices computed from crossing direction]
  ) <alg_crossdir2>

=== Step 3 - Sorting Intersection Points
    In this step we collect the region indices for all intersection points into a single list and sort the points by their distance $l$ from the ray start location, as illustrated in @alg_sort.  Next we compute the difference between adjacent points $Delta l$ which represents the intersection length of the ray with a particular voxel.  Note that this results in $Delta l$ equal to NaN or #sym.infinity for intersection points with $l=infinity$.
  #figure(
      image("figures/tomosphero/alg_sort.svg"),
      caption: [Intersection points along a LOS for all boundary types collected into a list and sorted by distance $l$ from ray start position, marked with x]
  ) <alg_sort>

    Next, we convert our list of region indices into a list of voxel indices by computing the voxel index of the ray starting location.  For each boundary crossed by the ray, the index of the next voxel changes in only a single dimension.

    This is illustrated in @alg_table (a), where the sorted region indices are inserted into an empty table starting with the ray starting point, then in @alg_table (b) missing indices are filled in from the previous row to form complete 3D voxel indices.

    We also account here for invalid region indices generated in step 2 where the ray passes outside the area defined by the grid or when boundaries have no ray intersection.  We achieve this by setting $Delta l=0$ for any rows where the region index is invalid (less than 0, or greater than or equal to the corresponding $N_r$, $N_e$, $N_a$), or where $Delta l$ is #sym.infinity or NaN.  Setting $Delta l=0$ negates any contribution of this row to the raytracing step.  This is shown in @alg_table (c).

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
                  table.header("r", "e", "a", $Delta l$),
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
                  table.header("r", "e", "a", $Delta l$),
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
                  table.header("r", "e", "a", $Delta l$),
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

    The final step of the raytracer is to perform a weighted sum of the object values at the voxel indices given in the table in the previous section with the $Delta l$ column as weights.

  // FIXME - use x instead of rho

    #math.equation(
        $sum_(i=0)^(2(N_r + 1) + \ 2(N_e + 1) + \ (N_a + 1) - 1) bold(rho)[r_i, e_i, a_i] * Delta l_i = innerproduct(rho, Delta l)$
    )

    where $bold(rho)$ is the object to be integrated and $r_i$, $e_i$, and $a_i$ are indices from the $i$th row of the table.


== API Overview <api_overview>
  We have designed this library to allow the user to easily construct tomography problems in just a few lines by instantiating classes that represent the grid, view geometries and forward raytrace operator.  These classes also support visualization via the ```python .plot()``` method to aid in validation and presentation.  In this section, we provide some some examples of setting up these classes for typical tomography problems as well as visualizations generated from the library.

=== Grid
    // FIXME: don't use theta/phi, switch to e/a
    The grid, which defines the physical extent and shape of the object to be traced, is created using the `SphericalGrid` class and may be defined by providing `shape` and `size` arguments when a uniform spherical grid is desired, as in @api_grid.  Alternatively, one may manually specify the grid boundary locations via `r_b`, `e_b`, `a_b` arguments for a non-uniform grid.

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
    A _view geometry_ is a logical collection of lines of sight.  This library provides a few built-in view geometries for common tomography paradigms shown in @api_conegeom and @api_parallelgeom.  In these examples, `pos` and `shape` control the detector position and shape.  Detector orientation is inferred automatically, but this behavior may be overridden.

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
            def model(c):
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
              loss = t.sum((y - op(model(c_hat)))**2)
              ...
            ```
        ),
        caption: [Pseudocode for reconstructing an object which lies in some transform domain defined by `model`],
    ) <api_reconstruction_latent>

  Finally, we remark that the differentiability of the operator allows it to be exploited in machine learning models, such as #glspl("PINN"), where the operator forms layers of the network itself @pinn.
=== Assumptions and Limitations

    Unlike most medical imaging paradigms, the detector is located at the vertex of `ConeRectGeom`.  We also assume that the detector lies outside of the grid, that rays extend to infinity, and that the object is zero outside of the grid.

    When combining view geometries into a collection, the geometries must have the same shape so that the corresponding stack of measurements returned during raytracing is a rectangular array.
    In the case of a dynamic object, the grid of the object being raytraced cannot move nor can its shape or extent change.

== Validation and Tests

    The raytracer's correctness has been thoroughly tested through both visual verification of various test objects from different positions (as shown in @visual_validation) and comparison against known analytic results.  Additionally, we have tested each phase of the raytracer algorithm outline described in @alg_outline, including the steps for deriving boundary intersection points, boundary crossing direction, and intersection point sorting.  These tests are integrated into the raytracer package's automated testing framework for evaluation when new software releases are made.

    #figure(
      grid(columns: 2, column-gutter: 1pt,
          subfigure(image("figures/limb_darkening.png", height: 10em), "samplert", "Solid sphere with limb darkening"),
          subfigure(image("figures/checkerboard.png", height: 10em), "samplert", "Hollow checkerboard shell"),
      ),
      caption: [Raytraced sample objects],
    ) <visual_validation>

  #rt([
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


= Static Retrieval of Exospheric Density <static_retrieval>
  #rt([
    - Chapter introduction section
        - This thesis makes a distinction between two types of retrieval algorithms: static, where hydrogen density is assumed to have no temporal component, and dynamic, where the density distribution is allowed to vary between vantages.
        - Static algorithms were introduced historically, but more recent algorithms have been developed that can handle a changing hydrogen distribution.
        - Static algorithms are important because they are conceptually and computationally simpler, can serve as a basis for more sophisticated dynamic algorithms, and still perform well for dynamic densities under quiet, non-storm conditions (see @static_dynamic for analysis)
        - This chapter introduces historical static retrieval approaches in order of increasing complexity, ending with a new static model contributed by this thesis.

    - 1D Retrievals
        - Early retrievals often relied on simple 1D models of exosphere
            - Either one-off measurements taken opportunistically (galileo flyby), limited viewing diversity
            - Avoids underdetermined system by limiting number of free parameters
            - Spherically symmetric assumption means good problem conditioning from any single vantage
            - Simpler models are computationally easier to retrieve
        - SOHO/SWAN (1996) @baliukin
            - Kinetic model of H atoms fit to observed data
            - Simple model which ignores radiation pressure (chamberlain) and extension which includes radiation pressure
                - Not quite 1D
            - Requires knowledge of exobase (e.g. derived from NRLMSIS model)
            - Applied onion-peeling technique in which outermost shells of model are fit to the data before proceeding inwards. high TP alt LOS pierce through few shells, low TP alt LOS pierce through more
        - Østgaard - IMAGE (2003) double exponential parametric form
            - Sph symmetric
            - $n(r) = n_1 e^(-r/alpha_1) + n_2 e^(-r/alpha_2)$
            - This specific formulation may have been chosen because of ease of computing analytic coldens/radiance
            - Density analytically converted to radiance given LOS and fit to data
        - PROCYON/LAICA (2015)
    - Gonzalo MAP estimate
        - Derived from earlier work from Butala
        - Not enough viewing diversity, stereoscopic
    - Zoennchen 2024
        - Functional form based on spherical harmonics
        - Defined at single shell with exponential decay (and other terms)
    - Spherical harmonic model
        - Applicability of spherical harmonic bases to modelling exospheres
            - Figure: show direct fits for different L and the max error on each
            - Conclusion: spherical harmonic basis is a good low rank representation of exospheric models
        - Splines enforce smoothness
            - Agnostic shape
        - Figure: spherical bases
        - Single measurement A00 retrieval
    - [ ] Longterm: summary table?
  ])


  This thesis makes a distinction between two types of retrieval algorithms: static, where hydrogen density is assumed to have no temporal component, and dynamic, where the density distribution is allowed to vary between vantages.
  Static algorithms were introduced historically, but more recent algorithms have been developed (see @dynamic_retrieval) that can accommodate a changing hydrogen distribution.
  Static algorithms are important because they are conceptually and computationally simpler, can serve as a basis for more sophisticated dynamic algorithms, and still perform well for dynamic densities under quiet, non-storm conditions (see @static_dynamic for analysis).
  This chapter introduces historical static retrieval approaches used in previous retrievals of exosphere density in order of increasing complexity.  However, many of these methods are either too restrictive in their assumptions about exospheric hydrogen distribution (e.g. assuming a specific functional form or even spherical symmetry) or require extensive hyperparameter tuning.  The chapter ends with a new static method in @spline_model contributed by this thesis, which attempts to address some of the issues of these other methods.
  Mathematical notation differs from the original publications to conform with notation used in this manuscript (see @inverse_problem).

  // This chapter presents several historical and recent methods of static retrieval of hydrogen density.  The retrieval methods are organized roughly in increasing complexity, ending with a novel method in @spline_model.


  == 1D Retrievals <1d_retrieval>


    Early retrievals often relied on simple spherically symmetric 1D models of the exosphere.  In cases where geocoronal studies were taken opportunistically (e.g. Galileo Earth flyby), measurements are often only available from a single vantage which has limited view geometry diversity.
    The assumption of spherical symmetry naturally produces a well-conditioned inverse problem from a single measurement taken at any vantage and avoids an underdetermined system (ill-posedness) by keeping model dimensionality low.

    A fundamental contribution to the field is the Chamberlain model, which is a spherically symmetric model derived from knowledge of motion of H atoms in the upper atmosphere (@earth_exosphere).  Under Liouville's theorem, the model makes some simple assumptions about the trajectories followed by exospheric H atoms and derives a density distribution from the resulting PDE.
    Baliukin et al. @baliukin applied the Chamberlain model to measurements from the SOHO/SWAN mission, extending it to take into account solar radiation pressure, loss of Lyman-α-scattering H atoms over time due to ionization, and prior knowledge of H density at lower altitudes from empirical models (e.g. NRLMSIS). The retrieval utilizes an onion-peeling technique where layers of the exosphere are retrieved one-at-a-time starting at the outermost layer and moving inwards.

    Other techniques include parametric 1D models found empirically and fit to data.  Østgaard et al. @ostgaard use a double-exponential model following the form

    #math.equation(
        $n_H (r) = n_1 "exp"(-r/alpha_1) + n_2 "exp"(-r/alpha_2)$
    )

    where $n_1, n_2, alpha_1, alpha_2$ are unknown model coefficients.  The two exponentials express the expectation that the H density distribution consists of two populations of H density atoms having distinct radial decay rates (known as _scale heights_).  Additionally, this simple exponential formulation permits an analytic line integration that Østgaard et al. utilize to perform model fitting directly in the measurement domain with no tomographic operator.  The exact method of fitting is not described, but most likely a linear regression similar to

    #math.equation(
        $
            bold(hat(c)) = arg min_bold(c) ||bold(y) - bold(y)_bold(c)||_2^2 \
            bold(hat(c)) = M(bold(hat(c)))
        $
    )

    where $bold(c) = {n_1, n_2, alpha_1, alpha_2}$ are the 4 free model parameters and $bold(y)_bold(c)$ is column densities derived by analytically integrating $M(bold(c))$.

    While computationally simple, these 1D solutions are not capable of capturing spherical asymmetries such as the geotail which are hypothesized to exist in the exosphere and are a major science target of Carruthers.  The next few sections build upon these ideas to create more sophisticated models that are useful when more measurement data is available from orbits designed to provide better view geometry diversity.

  == Spherical Harmonic Representation <sph_power>

    Zoennchen et al. @zoennchen_new present a retrieval method (derived from earlier work @zoennchen2011, @zoennchen2013, @zoennchen_old) featuring a partially separable parametric model based on a spherical harmonic representation (SHR).  Partially separable means the model can be expressed as a product of lower-dimensional functions which can decrease computational requirements.
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

    where $bold(c) = {{a_(0 0), ...}, {b_(0 0), ...}, {p_(0 0), ...}, {q_(0 0), ...}, c, k}$ are model parameters.  Zoennchen et al. truncate to $l ≤ 3$ for a total of $4 × 16 + 2 = 66$ free parameters. Note that $∀l : B_(l 0)(r) := 0$ by definition.

      The formulation of $"SHR"(r, e, a)$ in @scipy_sph differs slightly from the original in @zoennchen_new so that $Y_(l m)$ corresponds directly to Scipy's `sph_harm_y` function.

        Rewriting the model using this manuscript's notation and taking $M(bold(c)) := n_H (r, e, a)$ gives

      #math.equation(
          $
              hat(bold(c)) = arg min_bold(c) ||bold(y) - F(M(bold(c)))||_2^2 + lambda ||D_e M(bold(c))||_2^2 + lambda ||D_a M(bold(c))||_2^2 \
              hat(bold(rho))_"SHR" = M(hat(bold(c)))
          $
      )

      where $D_e$ and $D_a$ are finite difference operators for Tikhonov regularization that enforces elevational and azimuthal smoothness.

    == High Degree-of-Freedom <hdof>

      An alternative approach presented in Zoennchen et al. @zoennchen_new  @gonzalolaica and based on @solartomography1 is the #gls("HDOF") method which uses a non-parametric model where each voxel in the underyling density discretization is a separately optimizable parameter.  This method utilizes #gls("MAP") estimation which is optimal in a probabilistic sense but requires assumptions about the density distribution and noise statistics.  Specifically this method assumes that hydrogen density is a #gls("GMRF") distributed according to a known prior distribution $cal(N)(bold(rho)_"pr", Σ_"pr")$ with mean and covariance $bold(rho)_"pr"$ and $Σ_"pr"$ provided externally (e.g. a spherically-symmetric Chamberlain model) and measurements follow a noise distribution $Σ_"e"$.

      Under these assumptions, the #gls("MAP") solution has closed form

      #math.equation(
          $hat(bold(rho))_"MAP" = (L^T Σ_e^(-1) L + Σ_"pr"^(-1))^(-1)(L^T Σ_e^(-1) bold(y) + Σ_"pr"^(-1) bold(rho)_"pr")$
      ) <gonzalo_map>

      where $L$ is a matrix representation of forward operator $f$.

      An alternative formulation of @gonzalo_map is

      #math.equation(
          $
              hat(bold(rho))_"MAP" = bold(rho)_"pr" + underbrace(Σ_"pr" L^T (L Σ_"pr"L^T + Σ_e)^(-1), K)(bold(y) - L bold(rho)_"pr")
          $
      ) <gonzalo_map_alt>

      // #math.equation(
      //     $
      //         hat(bold(rho))_"MAP" = bold(rho)_"pr" + K(bold(y) - L bold(rho)_"pr") \
      //         K = (L^T Σ_e^(-1) L + Σ_"pr"^(-1))^(-1) L^T Σ_e^(-1)
      //     $
      // ) <gonzalo_map_alt>

      which illustrates that the #gls("MAP") solution is the density prior which has been shifted by residual measurement error $bold(y) - L bold(rho)_"pr"$ projected back on to the density domain.

      The matrix $K$ serves as a weighted pseudoinverse which balances between using prior information from $bold(rho)_"pr"$ versus measurement data $bold(y)$ based on uncertainty in the prior ($Σ_"pr"$) and measurement noise ($Σ_e$).

      Unfortunately, the formulation given in @gonzalo_map_alt is not always computationally tractible, as simply storing the result of $L Σ_"pr" L^T$ for a relatively small forward operator requires excessive memory (e.g. 10·64² LOS, 50³ voxels → 125 GB).  Zoennchen et al. solve this by observing that $Σ_"pr"^(-1)$ is sparse given locality assumptions of a GRMF (distant voxels are independent), and full computation of a dense matrix is unnecessary.  Another approach is to write the solution as a minimization problem @icon_inversion like

      #math.equation(
          $
              hat(bold(x))_"MAP" = arg min_(bold(x)) ||bold(y) - L bold(x)||_(Σ_e^(-1))^2 + ||bold(x) - bold(x)_"pr"||_(Σ_"pr"^(-1))^2
          $
      )

      or using this manuscript's notation

      #math.equation(
          $
              hat(bold(c)) = arg min_(bold(c)) ||bold(y) - F(M(bold(c)))||_(Σ_e^(-1))^2 + ||M(bold(c)) - bold(rho)_"pr"||_(Σ_"pr"^(-1))^2 \
              hat(bold(rho))_"MAP" = M(bold(hat(c)))
          $
      )

      where model $m$ is the identity operator for this non-parametric case.

      This minimization problem can be solved using a number of iterative methods, such as conjugate gradient descent or optimizers which utilize automatic differentiation (taking advantage of the differentiable forward operator presented in @raytracer).

  == Robust, Regularized, Positive Estimation

    Cucho-Padin et al., present a non-parametric method called #gls("RRPE") @rrpe @gonzalorrpe where voxels are allowed to vary independently. This model uses several 1D Tikhonov regularizers to enforce smoothness in each dimension.  The optimization problem is written

    #math.equation(
        $
            hat(bold(c)) = arg min_bold(c) ||bold(y) - F(M(bold(c)))||_2^2
            + lambda_r ||D_r M(bold(c))||_2^2
            + lambda_e ||D_e M(bold(c))||_2^2
            + lambda_a ||D_a M(bold(c))||_2^2 \
            hat(bold(rho))_"RRPE" = M(hat(bold(c)))
        $
    )

    where $D_r$, $D_e$ and $D_a$ are finite difference operators for radial, elevational and azimuthal directions, and $m$ is identity operator.  Cucho-Padin et al. use generalized cross validation to derive values for regularizer hyperparameters.

  == Spherical Harmonic Spline Model <spline_model>

    This section describes a new partially separable model which uses spherical harmonics as basis functions like in @sph_power.

    Previous models are either incapable of representing spherical asymmetries, require a density prior (problematic for an understudied atmospheric regime), or lack the interpretability of a spherical harmonic basis.
    This method replaces the relatively strong assumption of an inverse power law governing radial density variation with a more adaptable cubic spline model.
      This formulation is less strict about the shape of radial density profiles but still enforces smoothness of spherical harmonic coefficients.

    The model is written as

    #math.equation(
        $
            M(bold(c)) &= H(S(bold(c))) := sum_(l=0)^L sum_(m=-l)^l H_(l m)(e, a) S_(l m)(r) \
            S_(l m)(r) &= sum_(k=1)^K c_(l m k) B_(k 2)(r)
        $
    )

    #math.equation($
        H_(l m)(e, a) = cases(
        "Re"[Y_(l m)(e, a)] "if" m ≥ 0,
        "Im"[Y_(l m)(e, a)] "if" m < 0,
        )
      $)

    where $S_(l m)(r)$ is a cubic spline function composed of $K$ third order B-splines ${B_(1 2), ...}$ that define spherical harmonic coefficients for any radial shell from a small number of control points, illustrated in @sph_harm_spline.  This basis function formulation for $H_(l m)$ is equivalent to the one presented in @sph_power but avoids notational awkwardness when $m ≥ 0$ versus $m < 0$.


      #figure(
          pad(x:-50%,
              image("figures/harm_spline.svg", height: 15em),
          ),
          caption: "Smooth variation of coefficients in radial direction enforced by B-splines"
      ) <sph_harm_spline>

    The number and location of the $K$ control points used in the spline functions and the maximum spherical harmonic degree $L$ are selected at model initialization time.


    $bold(c) = {c_(0 0 0), ...}$ are the model parameters and the minimization problem is

    #math.equation(
        $
            hat(bold(c)) = arg min_(bold(c)) ||bold(y) - F(M(bold(c)))||_1 + lambda_1 ||S(bold(c))||_1 + lambda_2 ||"clip"_(-infinity, 0)(M(bold(c)))||_1 \
            bold(hat(rho)) = M(bold(hat(c)))
        $
    ) <sphharmspline_min>

    @sphharmspline_min uses an $L_1$ fidelity loss for robustness against outliers in measurements $bold(y)$, an $L_1$ regularization (LASSO) on spline-interpolated coefficients for model interpretability, and an $L_1$ penalization term on negative densities to enforce model correctness.
    Note that a choice of $L=0$ enforces spherical symmetry like the methods presented in @1d_retrieval, which can be useful for single snapshot retrievals.


    Because of the generalizability of spline density profiles and interpetability of spherical harmonics by the atmosphere science community, this is the primary method that will be used by the Carruthers mission for retrieving exospheric density during exospheric quiet conditions.

    // === Implementation Notes

    //   - Basis functions ${H_(l m), ...}$ may be computed once during initialization and used for all grid radial shells
    //   - Dimensions $l$ and $m$ should be flattened and merged for ${H_(0 0), ...}$ and ${c_(0 0 0), ...}$ to avoid an awkward pyramidal array structure
    //   - affine map for log-spaced control points

      #rt([FIXME: include summary table? TBD])

      #table(
          columns: (auto, auto, auto, auto),
          table.header([*Method*], [*Parametric*], [*\# Free Parameters*], [*Cost Function*]),
          [Chamberlain], [?], [], [],
          [Østgaard], [Yes], [], [],
          [SHR], [Yes], [66], [],
          [HDOF], [No], [N/A], [],
          [RRPE], [No], [N/A], [],
          [Spline], [Yes], [256], [],
      )


= Static Retrieval Validation <static_validation>
  #rt([
  - Intro paragraph
  - Reconstruction requirements
      - Contractual spatial resolution requirements and reporting interval (these are not precisely defined in some ways)
      - Precise mathematical interpretation of requirements
      - Exosphere not completely understood
          - must rely on models from physics simulations and prior retrievals from limited data
          - #link(label("datasets"))[(table) ground truth datasets]
  - Implementation Approach Justification
      - Temporal binning and #strike[Image Stacking] (reserve "image stacking" for on-orbit ops)
          - (not sure about the need for this section)
      - Temporal Baseline of Images
          - Static Algorithms on Dynamic Data
              - Static algorithms naturally induce an averaging effect on dynamic data.
              - #link("static_assumption", [(figure) Error introduced by static assumption on quiet-time data for various observation window durations])
              - Lara: which paper to cite?  Gonzalo storm time?
      - Science Pixel Binning
          - as mentioned previously, SPB reduces computational burden
          - at expense of some spatial resolution.
          - especially important is radial resolution - direction of largest gradients
          - figure: 1D error plot(s) of binned vs unbinned radiance profile
          - "to limit binned radiance error to 1%, we choose XXX radial bins "
      - Spherical Harmonic Spline Model Parameter Selection
          - figures of direct fits for different L and control points
              - note: this is just a guideline for determining minimum number of parameters to represent H distribution
              - error incurred during retrieval depends on well-posedness of the problem, which depends on model dimensionality
      - Avoiding aliasing and other sampling errors is a serious problem
          - Potential issues with density grid, LOS grid, science pixel binning
          - Refer to @discretization_considerations
          - Motivate choice of discretization grid
              - nyquist argument - 2x highest frequency of continuous model (Gonzalo thesis pg 52)
              - #link(label("stormbins"))[(table) storm time discretization]
              - #link(label("quietbins"))[(table) quiet time discretization]
  - Retrieval Performance
      - #link(label("codeoverview"))[(figure) retrieval block diagram]
      - Simulation block is used for validation prior to launch
      - Montecarlo simulation of reconstruction under noise
      - Performance Under Calibration Bias
          - Bias in g-factor, IPH, radiation, all affect accuracy of measurements
          - Retrieval algorithm should be able to cope with expected biases on orbit
  ])


  #rt("needs citations for each row")




  The Carruthers mission is required to prove that it is capable of meeting mission requirements set during its proposal.  With physics of the forward model, a retrieval algorithm implementation and knowledge of a hypothetical exosphere distribution, it is possible to justify that these requirements will be met through simulation.
  This chapter describes the retrieval algorithm performance requirements on the Carruthers mission, introduces synthetic ground truth datasets created by the exospheric research community, analyzes retrieval results of a few algorithms from @static_retrieval on a few of these datasets, and provides rationale for tunable settings used in the algorithms including model parameters and discretizations.

  == Reconstruction Requirements <recon_requirements>

    The Carruthers proposal stage set forth minimum requirements for the performance of the the thin exospheric density retrieval algorithms to ensure that reconstructions meet mission science objectives regarding the shape of the hydrogen distribution during geomagnetically quiet conditions (objective 1) and its transient response of the exosphere to impulsive geomagnetic storm events (objective 2).

    As the purpose of Carruthers is to make discoveries about an atmospheric regime which is not well-known, algorithm performance validation relies on testing against datasets which are derived from physics simulations or prior retrievals made from limited data.  @datasets goes into detail of available datasets and their origins.

    To meet science objective 1, Carruthers defines an "accuracy" requirement on densities retrieved from the datasets, constraining absolute error of every voxel in the retrieval to within ±50% of the ground truth.  Similarly, a "precision" requirement which is insensitive to bias in the retrieval ensures that temporal enhancements and depletions are measurable to within 20% of their actual value in the dynamic datasets.  These requirements also specify a minimum spatial reporting resolution of the retrievals, and either a 6 hour temporal resolution for geomagnetic quiet conditions or 1 hour for storm conditions.
    Requirements are specified for a single arbitrary voxel at time $t$ and position $r e a$ and are limited to exospheric regions above 3 Re in altitude and where densities exceed 25 atoms/cm³ to avoid problems of ill-posedness during inversion.
    // #rt([FIXME: The proposal does not specifically define a constraint on the confidence of the above requirements, but this will be considered in the next section.])
    The proposal specifies a confidence level requiring 84% success rate of requirements over an ensemble of trials.

    Precise definitions of requirements defined in the proposal are given in @requirements that ensure the reconstructed density distribution $bold(hat(rho))$ is representative of ground truth distribution $bold(rho)$ and that exospheric response to storms is detectable.

  #figure(
      table(
          columns: 2, align: horizon,
          table.header(
              "Requirement", "Criterion",
          ),
          table.hline(stroke: 2pt),
          // [Accuracy],
          // [#math.equation(
          //     $cases(
          //         (|hat(x)_"t,p" - x_"t,p"|) / x_"t,p" ≤ 50% &"if" x_"t,p" ≥ 25 "atom/cm³",
          //         "true" &"else"
          //     )
          //     $)],
          [Accuracy],
          [#math.equation(
              $(|hat(rho)_"trea" - rho_"trea"|) / rho_"trea" ≤ 50% "with 84% confidence"$,
          )],

          [Change],
          [#math.equation(
              $abs((hat(rho)_"(t+1)rea" - hat(rho)_"trea")/hat(rho)_"trea" - (rho_"(t+1)rea" - rho_"trea")/rho_"trea") ≤ 20% "with 84% conf."$
          )],

          [Temporal Resolution], [Storm: $Delta t = 1$ hr\ Quiet: $Delta t = 6$ hr],
          [Radial Resolution], [$Delta r ≤ 0.5$ Re],
          [Elevational Resolution], [$Delta e ≤ 180°$],
          [Azimuthal Resolution], [$Delta a ≤ 120°$],
          // [Additional constraints], [$r ≥ 3$ Re and $x_"trea" ≥ 25$ atom/cm³]

      ),
      caption: [Requirements for every reconstructed voxel $hat(rho)_"trea"$ where ground truth density $x_"trea" ≥ 25$ atom/cm³ and $r > 3$ Re.]
  ) <requirements>

    Using the forward model from @measurement_constraints and ground truth datasets, it is possible to simulate measurements that will be taken by Carruthers on-orbit.  With accurate simulated measurements, performance of retrieval algorithms can be analyzed in a Monte Carlo fashion against the mission requirements given in @requirements. The block diagram in @codeoverview shows an overview of an end-to-end validation of an arbitrary iterative retrieval algorithm against a ground truth dataset.
    Validation consists of four stages: a simulator for generating realistic noisy radiance measurements in #gls("DN"); a calibrator to convert sensor readings to interpretable column densities while removing instrument effects; a retrieval algorithm which iteratively fits a density to the column densities; a validator which checks that retrieved density $hat(rho)$ meets the requirements defined above.  On orbit, the simulator output will be replaced with actual downlinked measurements from Carruthers.

    #figure(
        // image("figures/code_overview.drawio.png", width: 90%),
        image("figures/thin_retrieval_flowchart.svg", width: 140%),
        caption: "Flowchart of simulator, calibration, retrieval loop, and validation test"

    ) <codeoverview>

    @quietbins and @stormbins give the grid used for retrieving densities, derived from the requirements in @requirements.

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
      caption: "Quiet Time Reconstruction Grid"
  ) <quietbins>

  #figure(
      table(
          columns: 4,
          table.header(
              "Dimension", "Range", "# of Bins", "Spacing",
          ),
          table.hline(stroke: 2pt),
          [time], [NA], [NA], [linear, 1 hr],
          [radius], [3-25 Re], [200], [logarithmic],
          [elevation], [0-180 deg.], [45], [linear, 3 deg.],
          [azimuth], [0-360 deg.], [60], [linear, 4.5 deg.],
      ),
      caption: "Storm Time Reconstruction Grid"
  ) #label("stormbins")

  == Datasets <datasets>

    The retrieval algorithm implicitly encodes assumptions about hydrogen density distribution (smoothness, spherical harmonic structures, etc.), and verification that these assumption hold requires testing against plausible ground truth densities.  This section overviews some of the ground truth datasets that are available for validation, their origins, spatiotemporal coverage, and model assumptions

  #rt([FIXME: Lara will give more descriptions about datasets])

    ==== Zoennchen & TWINS

    Using the parametric model described in @sph_power, Zoennchen et. al fit several density distributions to Lyman-α radiance measurements observed by TWINS @twins from specific days in both solar minimum and maximum conditions during geomagnetic quiet conditions.
    Notably, this dataset (referred to as the *Zoennchen* dataset) was generated from a spherical harmonic model with order L=3 and stronger assumptions about radial decay of hydrogen distribution.  It is included here for completeness, but validation against models with dissimilar assumptions is important to avoid an _inverse crime_ @inversecrime.

    Zoennchen et. al also present a second dataset (known as the *TWINS* dataset) derived from the model discussed in @hdof using the same measurement constraints as the *Zoennchen* dataset and a prior of a spherically-symmetric Chamberlain model with exobase parameters extracted from NRLMSIS 2.0 @msis.  This dataset is radially non-monotic in some high-altitude regions, which is not physically realistic based on current understanding of the exosphere.

  ==== Vidal-Madjar and Bertaux

  Another validation dataset which will be referred to as the *VMB* (Vidal-Madjar & Bertaux) model @vmb is a physics-based simulation of the exosphere.  It takes as input a 2D temperature distribution at the exobase (derived from NRLMSIS 2.0 @msis) and propagates this upwards and laterally.  By choosing different dates of the MSIS inputs, the model supports generating both geomagnetically quiet and storm-condition densities.  This simulation utilizes Louisville equations to propagate constraints from the exobase, ignoring physical processes which are known to influence exospheric hydrogen distribution such as solar radiation pressure and charge exchange. #rt([FIXME: wording?])

  @datasets_table gives an overview of these datasets.

  #figure(
      table(
          columns: 5,
          inset: 4pt,
          align: horizon,
          table.header(
              "Dataset", [Dynamic/\ Static], "Derived from", [Storm\ Condition], [Coverage],
          ),
          table.hline(stroke: 2pt),
          [*Zoennchen* @zoennchen_new], [Static], [Data\ (TWINS)], [Quiet], [3-25 Re],
          [Cucho-Padin \ *TWINS* @zoennchen_new], [Static], [Data\ (TWINS)], [Quiet], [3-25 Re],
          [Vidal-Madjar\ & Bertaux (*VMB*)], [Dynamic], [Simulation], [Quiet/Storm], [3-15 Re\ 3 months],
          [Connor\ *MSIS*], [Dynamic], [Simulation\ (MSIS)], [Quiet], [3-25 Re\ 3 weeks],
          [Connor\ *TIMEGCM*], [Dynamic], [Simulation\ (TIMEGCM)], [Storm], [3-25 Re\ 10 days],
          // [Cucho-Padin\ Dynamic], [Dynamic], [TWINS], [], []
      ),
      caption: "Synthetic ground-truth datasets"
  ) <datasets_table>

  #rt([FIXME: add section about connor datasets, or remove the rows])

  == Retrieval Performance <retrieval_validation>

    The Spherical Harmonic Spline method introduced in @spline_model will be used by the Carruthers mission for static retrieval of exosphere during geomagnetically quiet conditions.  @quiet_recon gives accuracy error for a single static density retrieval which assumes the hydrogen density is stationary during a 2 week observation window with evolving vantage.  This test was repeated on the Zoennchen and VMB datasets for #rt("50 trials (FIXME)") to establish that accuracy requirements are met at the necessary confidence level.

    Notably, the retrieval error is the smallest in the transverse (YZ) plane which is especially noticeable in the VMB retrieval.  This plane is oriented roughly perpendicular to measured #gls("LOS"), which intuitively should give good conditioning, and thus high-quality retrievals, to voxels in this region.

    Additionally, errors are larger in the low-altitude region near the midnight (-X) side of Earth, where no #gls("LOS") intersect due to occlusion by the Earth and optically-thick exosphere.  Voxels in this region are constrained only by spatial smoothing in the reconstruction model.

    The complete model parameters (e.g. control points, spherical harmonic order, etc.) and experimental settings (e.g. grid and view geometry resolution) used for reconstructing both of these cases is given in @quiet_recon_params.
    Experimental justification for model parameter selection is given in the next section.

  #let h = ("height": 12em);
  #figure(
      grid(
          columns: 1, column-gutter: 1em,
          subfigure(image("figures/zoennchen_err.png", ..h), "quietrecon", "Zoennchen reconstruction error"),
          subfigure(image("figures/vmb_err.png", ..h), "quietrecon", "VMB reconstruction error"),
      ),
      caption: [Static retrieval error during geomagnetically quiet conditions shown over 3 cardinal Cartesian planes.]

  ) <quiet_recon>

  #let h = ("height": 15em);
  #figure(
      grid(
          columns: 2, column-gutter: 1em,
          subfigure(image("figures/zoennchen_coeffs.png", ..h), "quietcoeffs", "Zoennchen retrieved coefficients"),
          subfigure(image("figures/zoennchen_direct.png", ..h), "quietcoeffs", "Zoennchen ground truth coefficients"),
          subfigure(image("figures/vmb_coeffs.png", ..h), "quietcoeffs", "VMB retrieved coefficients"),
          subfigure(image("figures/vmb_direct.png", ..h), "quietcoeffs", "VMB ground truth coefficients"),
      ),
      caption: [Static retrieval error during geomagnetically quiet conditions shown over 3 cardinal Cartesian planes.  Monte-Carlo trials show reconstructions meet accuracy requirement with sufficiently high confidence.]

  ) <quiet_recon_coeffs>

  #figure(
      table(
          columns: 2,
          inset: 4pt,
          align: horizon,
          "Problem Parameter", "Value",
          table.hline(stroke: 2pt),
          "Simulation Grid", [500x45x60],
          "Reconstruction Grid", [see @quietbins],
          "View Geometry", [see @camera_specs_sci],
          "Observation window", "14 days",
          "Harmonic Truncation Order (L)", "3",
          "Harmonic Spline C. Points", "16"
      ),
      caption: "Reconstruction parameters"
  ) <quiet_recon_params>

  #figure(
      table(
          table.header([Camera], [FOV\ (degrees)], [Resolution\ (pixels)], [Angular Res.\ (degrees)], [Spatial Res.\ (Re, projected)]),
          align: horizon,
          [Polar WFI], [18°], [50x100], [FIXME xxx radial\ yyy azimuthal], [xxx radial],
          [Polar NFI], [3.6°], [50x100], [xxx radial\ yyy azimuthal], [xxx radial],
          columns: (auto, auto, auto, auto, auto),
      ),
      caption: [Circular camera geometry specifications.\ Spatial resolution is projected onto Earth tangent plane as in @earth_fov(a)]
  ) <camera_specs_sci>



  == Implementation Approach Justification <static_justification>

  The purpose of this section is to justify configurable experimental parameters (i.e. excluding algorithm hyperparameters) used to evaluate retrieval performance.  This includes observation window, observation orbit location, science binning resolution, and density spherical grid resolution.  Improper choice of these parameters can introduce numerical errors that significantly impact retrieval performance.

    === Temporal Baseline of Stacked Images <static_dynamic>

    The algorithms described in @static_retrieval implicitly assume that densities being retrieved do not vary temporally within the observation window.  While this holds approximately true for quiet conditions when exospheric dynamics are slowly evolving, an averaging effect from this assumption contributes non-neglible error to the retrieval.  To upper-bound the error from this effect, it is possible to compute the worst case of the error between mean of the density within the window to an instantaneous density within the window.  An explicit expression is given below for voxel $r e a$,

    #math.equation(
        $e_(r e a)(t) &= max_(s in [t, t + Delta t])
            abs((limits("mean")_(tau in [t, t + Delta t]) x_(tau r e a) - x_(s r e a)) / x_(s r e a)) \
            e(t) &= max_(r, e, a) e_(r e a)(t)
        $
    ) <static_equation>

    // where

    // #math.equation(
    //     $"mean"_(tau in [t, t + Delta t]) = sum_(tau=t)^(t + Delta t) x_(tau r e a)$
    // )

        Where $e(t)$ considers the worst-case error across all voxels.  @static_assumption shows @static_equation above evaluated at all times in the VMB dataset for various $Delta t$.   Based on this data, Carruthers has chosen an observation window size of $Delta t = 14 "days"$ to constrain averaging effects to no more than 10% error.

  #figure(
      image("figures/scratch_static_assumption.png", width: 20em),
      caption: [Upper-bound of errors induced by temporal averaging effects of assuming a static scene. #rt([FIXME: rerun this for VMB dataset (and other datasets?)])]
  ) <static_assumption>

    === View Geometry and Grid Discretization <discretization>

    @post_processing introduced the notion of _science pixel binning_, a log-polar binning scheme which reduces the number of data constraints that need to be ingested by the retrieval algorithm while attempting to minimize loss of information in the measurement data.

    Choosing a science pixel binning scheme that has sufficient radial and azimuthal resolution is important to retrieval accuracy, and depends on the 3D discretization grid, the density being observed, and view geometry.  A choice of resolution which is too low introduces binning error into the retrieval algorithm which can dominate noise-induced error and lead to poor retrieval performance even under high #gls("SNR").

    The retrieval algorithm reduces computational complexity by substituting a set of native pixel measurements that fall within a science pixel with the mean and their associated #glspl("LOS") with a single #gls("LOS") at the science pixel centroid.  This operation implicitly assumes that column density varies linearly within a science pixel, but this assumption is violated by the radially-exponential nature of the hydrogen density distribution.
    @spb_bad_discretization illustrates a case where binned radial resolution is too low to capture large hydrogen density gradients at low altitudes, leading to overestimated column densities during retrieval.

    #figure(
        image("figures/scratch_spb_discretization.jpg", height: 15em),
        caption: [#rt([FIXME: placeholder]).  Large measurement quantization error incurred in regions with high gradients when science pixel binning is too coarse.  X axis is tangent point (TP) altitude.]
    ) <spb_bad_discretization>

      // FIXME: delete this.  after changing SPB to 100x50 and increasing grid, things are working OK now
    // A few potential solutions to this problem are summarized below:

    // \
    // 1. Correct column density overestimation analytically.
    // 2. Increase bin resolution to decrease overestimation effect.
    // 3. Mask out LOS which exhibit non-linear change over the science binned pixel.
    // 4. Redefine science pixel binning to avoid problematic #gls("LOS") directions entirely.
    // \

    // Solution 1 requires assuming a specific functional form to the column density distribution across a pixel, which is undesirable.  Solution 2 is also not practical as science pixel bin resolution is presumably already selected to use all available computational resources.  Solution 3 is feasible, but requires changes to how science pixel view geometries are generated.  Post-processing for Carruthers uses solution 4, as it eliminates the overestimation issue (at the expense of some lost measurements) and is the most straightforward to implement.

    // #figure(
    //     grid(
    //         columns: 2, column-gutter: 1em,
    //         subfigure(box(width:100pt, height:100pt, stroke:1pt), "spbbad", "Before mask"),
    //         subfigure(box(width:100pt, height:100pt, stroke:1pt), "spbbad", "After mask"),

    //     ),
    //     caption: [#rt([FIXME: placeholder]).  Error in column densities with problematic LOS eliminated with mask.]
    // ) <spb_bad_discretization_fix>

    // Based on experimentation like the one shown in @spb_bad_discretization_fix, Carruthers has selected a lower limit of 5 Re for science pixel #gls("LOS") tangent points to ensure column densities are sufficiently linear within a pixel and limit error incurred due to binning to no more than #rt([(FIXME: rerun experiment)])%.

    In a similar vein, an undersampled 3D grid that is coarse relative to the number of constraining #gls("LOS") (particularly in regions of high density gradients) can yield discretization artifacts in the raytraced images that significantly impact performance of the retrieval algorithm.  This is especially prevalent in the simulator, which uses native #gls("NFI") view geometry that has both high resolution and narrow #gls("FOV") resulting in a particularly dense #gls("LOS") sampling of the scene as illustrated in @los_discretization.

    #figure(
        image("figures/scratch_grid_discretization.jpg", height: 12em),
        caption: [#rt([FIXME: placeholder]).  Undersampled grid becomes apparent in discontinuities visible in measurements when LOS are sufficiently dense. X axis is tangent point altitude.]
    ) <los_discretization>


    The quantization error introduced by raytracing and binning for a particular grid and view geometry selection can be isolated by disabling instrument effects and comparing the output against an analytically computed column density.  Test density datasets like those in @datasets do not have a closed-form analytic column density, so a simpler model from @1d_retrieval is necessary.

      As the framework developed in this thesis allows for differing grids and view geometries in the simulator and retrieval, these pipelines (shown in @discretization_flowchart) must be validated separately against an analytically derived ground truth, shown in @discretization_flowchart.

    #figure(
        image("figures/discretization_err_flowchart.png", width: 30em),
        caption: [Discretization error of the simulator/retrieval pipelines can be analyzed by comparing their noiseless column density outputs against an analytic integration]
    ) <discretization_flowchart>

    @discretization_err_wfi and @discretization_err_nfi show the total error introduced by discretization in the simulator for several radial science pixel binning resolutions and radial grid resolutions, which are the main drivers of error.

    #rt([FIXME: Lara: plots are simulator only.  include results for retrieval stage?])

    #figure(
        image("figures/discretization_err_wfi.png", height: 40%),
        caption: "WFI Column density error induced by discretization of simulator 3D grid and view geometry.  Tests sweep over 3D grid radial bins and view geometry radial bins"
    ) <discretization_err_wfi>

    #figure(
        image("figures/discretization_err_nfi.png", height: 40%),
        caption: "NFI Column density error induced by discretization of simulator 3D grid and view geometry.  Tests sweep over 3D grid radial bins and view geometry radial bins"
    ) <discretization_err_nfi>

    @discretization_table summarizes the choice of 3D density discretization grid and science pixel binning resolution selected based on @discretization_err_wfi and @discretization_err_nfi and similar experiments for the retrieval stage, with experimentally swept parameters in bold.

    #figure(
        table(
            columns: 4,
            table.header([Algorithm Stage], [WFI View Geom.\ (rad.×azi.)], [NFI View Geom.\ (rad.×azi.)], [Grid\ (rad.×elev.×azi.)]),
            table.hline(stroke: 2pt),
            align: horizon,
            [Simulator], [*100*\×50], [*100*\×50], [*500*\×45×60],
            [Retrieval], [*100*\×50], [*100*\×50], [*200*\×45×60],
        ),
        caption: [Forward operator discretization parameters chosen to minimize errors while keeping in mind computational burden]
    ) <discretization_table>


        // image("figures/discretization_err_nfi.png"),

    In summary, the above analysis suggests that a science pixel resolution of 100×50 pixels and hydrogen density resolution of 200 radial bins during retrieval, while 500 radial bins are necessary during simulation ..... #rt([FIXME])

    === Reconstruction Model Parameters

      #let bdiag_args = ("width": 7em, "left_sep": 1em)
      #let zoennchen_param_table = table(
          columns: (bdiag_args.width, auto, auto, auto),
          align: horizon + center,
          bdiagbox(..bdiag_args)[*L*][*cpoints*], [8], [12], [16],
          ..csv("zoennchen_params.csv").flatten(),
      );
      #let vmb_param_table = table(
          columns: (bdiag_args.width, auto, auto, auto),
          align: horizon + center,
          bdiagbox(..bdiag_args)[*L*][*cpoints*], [8], [12], [16],
          ..csv("vmb_params.csv").flatten(),
      );
      #figure(
          grid(columns: 2, gutter: 1em,
              subfigure(zoennchen_param_table, "dfit", [Zoennchen Dataset (2015)]),
              subfigure(vmb_param_table, "dfit", [VMB Dataset (quiet)]),
          ),
          caption: "Direct Fit Maximum % Error over all voxels. Sweeping spherical harmonic truncation order and number of control points"
      ) <parameter_sweep>





= Proposed Future Work <future_work>

  == Uncertainty Quantification
  #rt([

  - analysts rely on knowledge of uncertainty to make accurate decisions when reasoning about model predictions from real data
  - example: mass predicted in a patient based on MRI data.  knowledge of uncertainty in prediction will be a factor in whether to operate or take further tests to improve confidence in presence of the mass
  - in context of carruthers: reporting error margins can assist scientists in determining if features in retrievals or artifacts of the model prediction
  - mission requires derivation of uncertainty with reported reconstructions
  - should take into consideration measurement noise and forward operator (determined by view geometry) to estimate posterior
  - montecarlo simulation is an approach - can estimate distribution of measurements, run ensemble of reconstructions
    - impractical when retrievals are expensive

  - epistemic uncertainty - sometimes called "model uncertainty". generally can be improved with access to more training data
  - aleatoric uncertainty -

  - important note: uncertainty quantification does not attempt to quantify error due to model mismatch ()
      - https://tristanvanleeuwen.github.io/IP_and_Im_Lectures/statistical_perspective.html
  - foo
  ])

  Analysts rely on knowledge of uncertainty to make accurate decisions when reasoning about predictions from models.  For example, a doctor looking at a a mass detected in a patient will want to the the level of certainty in the prediction before deciding to operate or conduct additional tests.

  In the context of Carruthers, knowledge of model uncertainty will assist physicists studying the exosphere in determining if certain features in exospheric retrievals are real or artifacts of model prediction.  Carruthers reporting requirements include a map of expected errors to be delivered with retrievals.

  Literature often divides uncertainty in model predictions into two types: _epistemic_ and _aleatoric_.  To avoid an overly broad or philosphical description of uncertainty, the following definition focuses on the case of a model $m_theta$ whose parameters are learned in a supervised manner from a dataset $D$:

  - _epistemic uncertainty_ - Uncertainty due to training dataset, either due to noise in the dataset or lack of coverage to all inputs.  Generally can be reduced by introducing more data into training dataset.
  - _aleatoric uncertainty_ - Uncertainty in measurement constraints which are propagated through the retrieval algorithm.

  #figure(
      grid(columns: 2, column-gutter: 2em,
          subfigure(image("figures/epistemic_scratch.png"), "epiale", [Epistemic uncertainty due to limited training data]),
          subfigure(image("figures/aleatoric_scratch.png", height: 12em), "epiale", [Aleatoric uncertainty due to noise in measurement constraints]),
      ),
      caption: [Simple example demonstrating different types of uncertainty in the solution a 1D inverse problem]
  ) <uncertainty_types>

  == Dynamic Retrieval <dynamic_retrieval>
  - dynamic system identification
      - greybox/semi-physical modelling - incorporate justifiable knowledge about system without going into too much detail
      - whitebox modelling - complete system identification

// = Dynamic Retrieval Validation
//   - reconstruction requirements
//       - spatial resolution requirements and reporting interval - same as static
//       - change detection requirement

= Conclusion

  #set heading(numbering: "A.1.1", supplement: [Appendix])
  #counter(heading).update(0)

= Appendix - Instrument Simulation and Background Subtraction <appendix_sim>

  // #rt([FIXME: delete this paragraph]) This section describes a statistical model for the instrument noise and background signals present in the NFI and WFI cameras during measurement of exospheric Lyman-α.  Modelling these processes is important for converting raw sensor measurements in digital numbers (DN) as telemetered by the spacecraft to corresponding radiances that can be used for tomographic reconstruction, summarized in @calibration.   A statistical model is also important for generating synthetic noisy measurements to validate the performance of retrieval algorithms.  As a result, the Carruthers cameras have undergone extensive laboratory characterization to determine instrument model parameters and periodic on-orbit calibration is planned to account for parameter drift due to exposure to the space environment.

= Instrument Model and Central Limit Theorem <instrument_clt>


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

  #math.equation(
      $Z tilde.op cal(N)(M N^2mu_Y, M N^2sigma_Y^2)$
  )

  where

  #math.equation(
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
  ) <clt_mean>

  To validate the accuracy of the approximation, we compare a Monte Carlo simulation of @stackingequationdirect to a Normal distribution with parameters given above by the #gls("CLT"), shown in @clt. (#rt("FIXME: give specific parameters used in simulation below?  at least mention if photon rate is worst case"))

  #figure(
      grid(columns: 2, column-gutter: 1pt,
          subfigure(image("figures/montecarlo_1_NFI.png"), "mcclt", [NFI, 2x2 binning, 1 frame]),
          subfigure(image("figures/montecarlo_1_WFI.png"), "mcclt", [WFI, 4x4 binning, 1 frame]),
          subfigure(image("figures/montecarlo_14400_NFI.png"), "mcclt", [NFI, 2x2 binning, 14400 frames]),
          subfigure(image("figures/montecarlo_28800_WFI.png"), "mcclt", [WFI, 4x4 binning, 28800 frames]),
      ),
      caption: [Monte Carlo simulation of actual distribution of $Z$ vs #gls("CLT") approximation]
  ) <clt>

  //   == Berry-Esseen Bound

  //   #rt([Section incomplete.  Derivation of $ex(|Y^3|)$ is very complicated but necessary for Berry-Essen bound which gives theoretically guarantees about the convergence of the blue→red lines in the plots above.  May remove this section])

  //   #math.equation(numbering:none,
  //       $ex(Y^2) = var(Y) + ex(Y)^2$
  //   )

  // - #rt("WIP: CLT verification experiments")
  // - #rt("WIP: Berry-Esseen Theorem")

  //   #image("figures/clt_wfi.png", height: 15em)
  //   #image("figures/clt_nfi.png", height: 15em)

    == Calibration and Subtraction <calibration>

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

  #rt([FIXME: remove this if not using Berry-Esseen bound])

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

= Appendix - Discretization Considerations <discretization_considerations>

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

= Appendix - Extra Physics <appendix_extra_physics>

  - #rt([come up with a better name - "second-order correction terms"])


= Glossary



  // - *#acr("LOS")* - #acrfull("LOS")

  #print-glossary(glossary)

  // - *MCP* - microchannel plate 23, 24
  // - *SNR* - signal-to-noise ratio 26, 27
  // - *UV* - ultraviolet: description goes here 22, 23



#bibliography("refs.bib")