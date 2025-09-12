# Project Proposal

## Project Goal

Recently, a low signal gravitational wave event was reported by LVK ([S250818k](https://gcn.nasa.gov/circulars/events/ligovirgokagra-s250818k)). While it is likely that this detection is terrestrial in nature, if the source is indeed astrophysical, then preliminary analysis places a $>99%$ likelihood on a primary mass $<1 \,M_\odot$ with a sub-solar mass secondary. If confirmed, this would be the first sub-solar mass object detected via gravitational waves. In addition, a transient EM signal was detected within the sky localization band of the source, which would make S250818k the second multi-messenger event after GW170817.  

If a compact object of $<1\,M_\odot$ were confirmed, it would require an overhaul of current understanding of the formation of such objects. As such, it is imperative to check all possible sources for an apparent sub-solar mass detection. We seek to explore the possibility that a Neutron Star-White Dwarf (NS-WD) binary may produce a S250818k-like signature. 

We seek to develop a semi-analytic NS-WD binary model to test this possibility. In particular, we aim to implement 2.5 Post-Newtonian effects alongside a point mass treatment of the stars and a fluid mechanical treatment of non-conservative unstable mass transfer. Our ultimate goal is to produce gravitational waveforms which can be analyzed following the LIGO pipeline to determine measured parameters from such an event, in particular whether such the WD companion may be labeled as sub-solar mass either correctly or incorrectly. If time allows, we may investigate EM counterparts due to super-Eddington accretion on the NS primary. 

## Description of the inputs (0.5 point)

- Initial orbital parameters
  + NS and WD mass
  + Semi-major axis
  + Eccentricity

## Description of the outputs (0.5 point)

- Gravitational waveform
- Derived parameters from waveform
  + Chirp mass up to inclination
  + Optimistically, primary and secondary mass
  + Luminosity distance
- EM counterpart signal (time permitting) 

## How will you test that your code is giving accurate results? (1 point)

INSERT TEXT

## What will you vary to explore performance as a function the problem size?  Describe the relevant range of problem sizes. (1 point)

INSERT TEXT

## Complete list of what computer architectures you will use, as well as the justification for your choices (1 point)

Computer architectures to be used:
- Serial (i.e., one core of one CPU)
- Multi-core, Shared-memory (i.e., one compute node with multiple CPU cores which have access to a shared memory system)
- CPU + GPU: (avaliable via Lynx cluster)

Rationale for above choices:

Serial and multi-core, shared memory architectures are required components of the project. More broadly speaking, a serial architecture is an out-of-box architecture that is sufficient for normal day-to-day operations and is meant for code that doesn't need to be parallelized. Multi-core, shared memory is an architecture used for parallelized code and will be our second version of code for the project, but it has issues with different threads trying to write to the same memory location simultaneously. Out of the options for the third architecture, it is between CPU + GPU and multi-core, distributed memory (the others aren't readily available in Lynx and require extra steps to set up).

For our improved parallelized version of the multi-core, shared memory architecture, we choose the CPU + GPU architecture due to ease of implementation in Julia and its fit to the equations of fluid mechanics. Specifically, using a CPU + GPU architecture removes the need to explicitly split arrays into subarrays and assign each to different cores, which a multi-core, distributed memory system would require. CPU + GPU architectures also don't need explicit communication between cores to be implemented by the programmer, whereas a multi-core, distributed memory system does. Therefore, debugging is easier with a CPU + GPU architecture. Furthermore, a multi-core, distributed memory system is best when large amounts of memory are required for each core, but with the fluid mechanics implementation we envision, the size of our 2D/3D grids probably won't necessitate such high memory costs.

## Complete list of what programming languages, libraries and packages you will use, as well as the rationale of your choices.  (1 point)

Programming language(s) to be used:

Julia 

List of any libraries/packages that you know will be needed for your project:  

- Base.Threads for parallelization (in Julia)
- CUDA.jl (for NVIDIA GPUs) or AMDGPU.jl (for AMD GPUs) or oneAPI.jl (for Intel GPUs) (in Julia)
- DifferentialEquations.jl for ODE/PDE solvers
- Plots.jl for visualization

Optionally, list of any functionality that you need and suspect that you can find in an existing library/package, but you have not yet identified a specific library/package to meet that need:

NA

Rationale for above choices:

- Base.Threads is needed for parallelization
- CUDA.jl / AMDGPU.jl / oneAPI.jl is needed for utilizing GPUs (our chosen architecture is CPU + GPU). The specific package used will depend on the GPU we have access to.
- DifferentialEquations.jl is necessary for the differential equations in fluid mechanics.
- Plots.jl is necessary for visualization of results, tests, etc. 

## Collaborators 
Will you be collaborating with any other students on the project?  If so, who?  

The team is composed of Tristan Weaver and Victor Liu.
