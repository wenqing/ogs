---
title: OpenGeoSys
authors:
  - Lars Bilke

date: 2017-01-13T14:24:23+01:00

hero:
  headline: Open-source multi-physics
  textline: OpenGeoSys (OGS) is a scientific [open source project](https://gitlab.opengeosys.org/ogs/ogs) for the development of numerical methods for the simulation of thermo-hydro-mechanical-chemical (THMC) processes in porous and fractured media. Current version is OpenGeoSys-6 which is documented on this page. For information about OpenGeoSys-5, see [its dedicated section](/ogs-5). OGS has been successfully applied in the fields of regional, contaminant and coastal hydrology, fundamental and geothermal energy systems, geotechnical engineering, energy storage, CO2 sequestration/storage and nuclear waste management and disposal.
  quicklinks:
    - name: Tools
      anchor: tools
    - name: HPC
      anchor: hpc
    - name: Vis
      anchor: vis

feature_intro: OpenGeoSys' adaptable and modular architecture enables a wide variety of use cases and flexible workflows. In the following we highlight some of its most important features.

features:
- headline: Comprehensive Pre-Processing Tools
  textline: |
    A wide range of helper tools exist to get your model up and running with OpenGeoSys.

    Convert your existing data sets into appropriate OGS data formats and structures.

    Create meshes approximating geometrically the domain of interest. Analyze mesh quality, cleanup the mesh or adding layers to it.

    Parametrize the model with material parameters, boundary conditions and source terms.
  visual:
    permalink: "/docs/tools/meshing-submeshes/extract-surface/TopBottomSideSurface.png"
    alt: Extracted surfaces
  # links:
  #   - text: See Docs
  #     url: "/docs/tools"
  layout: left
  anchor: tools

- headline: Process Coupling
  textline: |
    A coupled system of equations can be either solved in a fully coupled way of the monolithic method, or in the sequential manner of the staggered scheme. The monolithic scheme is applied for all coupled processes, while the staggered scheme are available for the coupled processes of thermo-hydraulic, hydro-mechanical, and phase field mechanical problems.
  visual:
    permalink: /features/StaggeredCouplingScheme.png
    alt: Staggered coupling scheme
  layout: right

- headline: Data integration
  textline: |
    Integrate and visualize data sets for OpenGeoSys by using the OpenGeoSys Data Explorer. It provides functionality to visually assess the data and see possible artefacts, inconsistencies between data sets or missing information.
  layout: left
  carousel:
    slides:
      - permalink: /features/vis/chaohu_paper_mesh.png
        caption: Foo
      - permalink: /features/vis/ogsgui.png
        caption: Bar
      - permalink: /features/vis/DEMSelke3D.png
    play: false
#  links:
#    - text: "Learn More: Data Explorer"
#      url: "/features/data-explorer"

- headline: Visualize results
  textline: |
    By using VTK data formats visualizing simulation result data sets becomes an easy task. The de-facto standard software for scientific visualtions [ParaView](https://www.paraview.org) can be used to explore and analyze complex data in a visual way.

    [Virtual reality enabled visualization](https://www.ufz.de/vislab) brings your data onto the large screen for intuitive exploration and assessment.
  layout: right
  carousel:
    slides:
      - permalink: /features/vis/layeredview2.png
        caption: Bar
      - permalink: /features/vis/vislab.png
        caption: Interactive exploration of river flow phenomenae at the TESSIN VISLab of the Helmholtz Centre for Environmental Research – UFZ
      - permalink: /features/vis/contours2-bw.png
    play: true
  anchor: vis

- headline: High performance computing
  textline: |
    High performance computing (HPC) has became a necessity in the modelling of environmental and geotechnical problems for better characterization of the complexity of geo-systems as well as predicting their evolution in time. Parallel computing is the most efficient method in the high performance computing. In OGS, the parallalization of the finite element (FE) computation is based on the domain decomposition method (DDC).

    Decomposed global matricies and vectors are handled by PETSc and the system of linear equations are solved by the performant PETSc solver. PETSc builds upon the Message Passing Interface (MPI) suitable for a wide variety of parallel computing architectures.

    Parallelization is implemented for single processes as well problems with coupled processes which are using the same order of element for each process.
  visual:
    permalink: "/features/HPC-DDC.png"
    alt: Domain decomposition for parallel processing
    rounded: true
# links:
#    - text: Learn more
#      url: "/features/hpc"
  anchor: hpc

- headline: Transparent development workflows
  textline: |
    OpenGeoSys is an open-source project developed by a community of researchers. We try to be
    open-minded and and make team decisions. We try to help users and developers as best as we can.

    We invite you to take part in this journey, shape the future of OpenGeoSys together and happily welcome any contribution.
  subfeatures:
    - headline: Setup a development environment
      textline: |
        [Learn how to](/docs/devguide) obtain the source code, how to install required other software (e.g. compilers and tools), how to configure the software and how to generate the application binary.
    - headline: Contribute code
      textline: |
        Implement your new feature and [let the CI system](/docs/devguide/development-workflows/continuous-integration) run sophisticated tests automatically for you incorporating multiple computing platforms, a magnitude of software configurations and a whole array of CPU intensive complex test simulation runs.
    - headline: Get help from core developers
      textline: |
        Once your feature is ready the [code review process](/docs/devguide/development-workflows/code-reviews/) starts. A helpful [core developer](https://gitlab.opengeosys.org/ogs/ogs/-/graphs/master) checks the proposed change for general acceptance and may give hints for improvement (of e.g. the computational performance or the code structure). Once the iterative feedback loop between you, code reviewer(s) and the automated test system satisfies all aspects the proposed change is merged into the main development line.
  visual:
    permalink: "/features/OGS-Software-Engineering-Small.png"
    alt: Dev workflow
  layout: right
#   class: inverse

- headline: Ready to dive in?
  textline: |
  subfeatures:
    - headline: Using OpenGeoSys
      textline: |
        Start using OpenGeoSys by [downloading](/releases) a prebuilt package.
#       links:
#         - text: <i class="far fa-download"></i> Download OpenGeoSys
#           url: /releases
    - headline: Developing OpenGeoSys
      textline: |
        Getting started developing OpenGeoSys at the [Developer Guide](/docs/devguide).
#      links:
#        - text: <i class="far fa-book"></i> Read the Developer Guide
#          url: /docs/devguide
    - headline: Become part of the Community
      textline: |
        Get in touch with the OpenGeoSys Community via our [Discussion forum](https://discourse.opengeosys.org), [GitLab](https://gitlab.opengeosys.org/ogs/ogs) or by [email](mailto:info@opengeosys.org).
#      links:
#        - text: <i class="fab fa-gitlab"></i> GitLab
#          url: https://gitlab.opengeosys.org/ogs/ogs
  layout: vertical
#   class: inverse
---
