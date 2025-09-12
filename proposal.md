# Project Proposal

## Project Goal

Recently, a low signal gravitational wave event was reported by LVK ([S250818k](https://gcn.nasa.gov/circulars/events/ligovirgokagra-s250818k)). While it is likely that this detection is terrestrial in nature, if the source is indeed astrophysical, then preliminary analysis places a $>99%$ likelihood on a primary mass $<1 \,M_\odot$ with a sub-solar mass secondary. If confirmed, this would be the first sub-solar mass object detected via gravitational waves. In addition, a transient EM signal was detected within the sky localization band of the source, which would make S250818k the second multi-messenger event after GW170817.  

If a compact object of $<1\,M_\odot$ were confirmed, it would require an overhaul of current understanding of the formation of such objects. As such, it is imperative to check all  possible sources for an apparent sub-solar mass detection. We seek to explore the possibility that a Neutron Star-White Dwarf (NS-WD) binary may produce a S250818k-like signature. 

We seek to develop a semi-analytic NS-WD binary model to test this possibility. In particular, we aim to implement 2.5 Post-Newtonian effects alongside a point mass treatment of the stars and a fluid mechanical treatment of non-conservative unstable mass transfer. Our ultimate goal is to produce gravitational waveforms which can be analyzed following the LIGO pipeline to determine measured parameters from such an event, in particular whether such the WD companion may be labeled as sub-solar mass either correctly or incorrectly. If time allows, we may investigate EM counterparts due to super-eddington accretion on the NS primary. 

## Description of the inputs (0.5 point)

- Initial orbital parameters
 + NS and WD mass
 + Semi-major axis
 + Eccentricity

## Description of the outputs (0.5 point)

- Gravitational waveform
- Derived parameters from waveform
 + Chirp mass (primary and secondary mass if degeneracies can be broken) 
 + Luminosity distance
- EM counterpart signal (time permitting) 

## How will you test that your code is giving accurate results? (1 point)

INSERT TEXT

## What will you vary to explore performance as a function the problem size?  Describe the relevant range of problem sizes. (1 point)

INSERT TEXT

## Complete list of what computer architectures you will use, as well as the justification for your choices (1 point)


Edit the list below (most deleting lines that you won't use) as appropriate for your project.

Computer architectures to be used:
- Serial (i.e., one core of one CPU)
- Multi-core, Shared-memory (i.e., one compute node with multiple CPU cores which have access to a shared memory system)
- For a third architecutre, I anticipate that most students will pick one of the following:
  + Multi-core, Distributed-memory  (i.e., multiple CPU cores, each of which has access to its own memory system avaliable via Lynx cluster)
  + CPU + GPU: (avaliable via Lynx cluster)
  + CPU + Multiple GPUs: (potentially avaliable via Roar Collab, but may require a special request)
  + CPU + TPU: (see [this post](https://gweb-cloudblog-publish.appspot.com/products/ai-machine-learning/now-you-can-train-ml-models-faster-and-lower-cost-cloud-tpu-pods/amp/) for info, but would be on own for setting up access to tensor processing units)
  + Multi-core, Distributed-memory: Cloud (e.g., JuilaHub, AWS, Azure, Domino Data Lab; but you'd be on your own for setup and any cost)

Rationale for above choices:

I anticipate that most students will be able to have their computational needs met by the ICDS Lynx cluster.
If you anticipate needing additional computing resources, describe your needs in more detail. 
Do you need help getting access to those resources?

## Complete list of what programming languages, libraries and packages you will use, as well as the rationale of your choices.  (1 point)

Programming language(s) to be used:

Julia 

List of any libraries/packages that you know will be needed for your project:  

INSERT TEXT

Optionally, list of any functionality that you need and suspect that you can find in an existing library/package, but you have not yet identified a specific library/package to meet that need:

INSERT TEXT

Rationale for above choices:

INSERT TEXT

## Collaborators 
Will you be collaborating with any other students on the project?  If so, who?  

INSERT TEXT

