Cleo vs PySDM 1-D KiD: setup comparison status
===============================================

This page doesn't document how to run a model with the KiD. Instead it lists
the similarities and differences between the Cleo and PySDM KiD Setups.

TL;DR Most of the dynamics and initial superdroplets now match between Cleo and PySDM;
the collision kernel and fall speeds are the largest remaining differences.


The same in Cleo and PySDM
---------------------------

**Dynamics**

- **Updraught:** the same ρ\_d·w(t) sinusoid, giving the same Courant number
  rhod_w(t + dt/2)·dt/dz, and the same ρ\_d g-factor.
- **Advection scheme:** ``cmn.get_options(is_accurate=True)`` in both.
- **Initial profiles:** q\_v and θ come from ``common.py``; ρ\_d comes from the
  same ``cmn.hydrostatic_equilibrium`` by default. With ``is_exner_novapour``
  we use θ_d = θ, yet calculate the pressure with accounting for vapour pressure.
  This means also with ``is_exner_novapour``, Cleo's ρ\_d, T, and P match those
  calculated by PySDM, however accounting for vapour in the pressure differs from
  the original Shipway and Hill (2012) KiD (see ``exner_function_thermodynamics``).
- **θ\_d is fixed in time in both.** Latent heating is felt inside a coupling step
  and then thrown away: Cleo's is reset by ``set_thermo``, PySDM's by ``sync()``
  with ``update_thd=False``.
- **Vapour:** changes from condensation are kept from one step to the next, and
  neither model feeds an advective tendency into its condensation substeps.

**Initial superdroplets**

- **Aerosol spectrum:** default is a lognormal with 40 nm geometric mean radius
  and geometric standard deviation 1.4, sampled over the same 9.52–168 nm
  range (``cmn.get_initial_spectrum_and_size_range``).
- **Number concentration:** a uniform 50 cm⁻³ when ``xi_by_pressure`` and
  ``xi_by_dry_air_density`` are both ``False``.
- **Radius sampling:** α = 1.0 alpha sampling with quasirandom radii, and
  independent seeds for radii, positions and collisions.

**Diagnostics:** cloud water is 1–40 µm and rain is ≥ 40 µm.

**Collisions:** coalescence efficiency is 1 in both.


Differences current settings cannot remove
--------------------------------------------

- **Collision kernel**
    - Cleo (``LongHydroProb``): geometric kernel × Long's collision efficiency (Simmel et al. 2002), with Simmel fall speeds
    - PySDM (``Hydrodynamic`` or ``Long1974``): Berry 1967 efficiency or Long 1974 polynomial fit; neither matches Cleo
    - Probably this is the largest difference, especially on when rain starts and how much
- **Fall speeds**
    - Cleo (Simmel et al. (2002) in the collision kernel; Rogers / Gunn & Kinzer 1949 for sedimentation)
    - PySDM (Gunn & Kinzer 1949; ``Formulae`` default)
    - Affects collisions and to a probably negligible extent sedimentation
- **if ``xi_by_pressure != False and xi_by_dry_air_density != False``**
    - (###WIP) the two models don't scale the number concentration the same way yet.
- **Particle reservoir below the surface**
    - Cleo uses circa. 25m grid-box which refreshes each time-step
    - PySDM uses gridboxes down to circa. -950m which supplies aerosol from below
    - may cause minor difference in surface ``qvap`` and initial (haze) particle sizes
- **Superdroplet positions**
    - Cleo initialises positions per gridbox
    - PySDM initialises random positions throughout the whole domain
    - PySDM may have extra noise in droplet positions (about ±6% scatter in number per gridbox)

For a fair comparison, the kernel and fall-speed mismatch matters most. To match
Cleo exactly, PySDM would need a custom kernel with Long's efficiency and Simmel
fall speeds.


Differences in how output is recorded
------------------------------------------

These need care when comparing the two models' output.

- **Order within a step.** Cleo advects and then runs the microphysics. PySDM does
  the reverse, and its vapour, RH, T and p products are recorded *before* that
  step's advection, so PySDM's vapour output lags by one advection step.
- **Temperature output.** Cleo records ``temp`` after its microphysics, so it
  includes up to 1.25 s of latent heating. PySDM's ``T`` comes from the fixed θ\_d.
