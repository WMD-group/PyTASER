---
title: 'PyTASER: A Python package to simulate transient absorption spectroscopy (TAS) for bulk crystals from first principles'
tags:
  - Python
  - materials science
  - first-principles
  - optics
  - spectroscopy
  - density functional theory
authors:
  - name: Savyasanchi Aggarwal
    equal-contrib: false
    affiliation: "1, 4" # (Multiple affiliations must be quoted)
  - name: Young Won Woo
    equal-contrib: false
    affiliation: "2,4"
  - name: Seán R. Kavanagh
    orcid: 0000-0003-4577-9647
    equal-contrib: false
    affiliation: "1, 4" # (Multiple affiliations must be quoted)
  - name: Alex M. Ganose
    orcid: 0000-0002-4486-3321
    equal-contrib: false
    affiliation: "3" # (Multiple affiliations must be quoted)
  - name: Aron Walsh
    orcid: 0000-0001-5460-7033
    corresponding: true
    affiliation: "4"
affiliations:
 - name: Department of Chemistry, University College London, London, United Kingdom
   index: 1
 - name: Department of Materials Science and Engineering, Yonsei University, Seoul, Korea
   index: 2
 - name: Department of Chemistry, Imperial College London, London, United Kingdom
   index: 3
 - name: Department of Materials, Imperial College London, London, United Kingdom
   index: 4
date: 22 July 2023
bibliography: paper.bib
---

# Summary


PyTASER is a Python package for simulating transient absorption spectra for pristine crystalline materials from first-principles calculations. This facilitates mapping between the electronic structure of the material to its optical properties, using the easily accessible bandstructure and electronic density-of-states of the material, both of which can be obtained through density functional theory. It preferentially allows input of locally-obtained data, but is also interfaced with the Materials Project database to allow use by non-computationally familiar experimentalists, or high-throughput method users. 

Transient Absorption Spectroscopy (TAS) is a powerful tool used to characterise the properties and time-evolution of electronically excited states of materials. This is achieved by comparing the optical absorption spectra between the ground 'dark' state, and a monochromatic-laser excited 'light' state, with a variable time delay between the two measurements to incorporate time evolution of the excited state. 
Thus, TAS is widely used to understand microscopic processes in photochemical and electrochemical transformations, including phenomena such as electron trapping and carrier recombination, which are important properties when considering technological application of a material. 
The drawback of modern TAS is that the experimental setup for TAS is expensive and convoluted, and spectra are difficult to interpret. This is especially the case for crystals, where specific valence and conduction band structures can give rise to complex features. As a result, it is difficult for experimentalists to compare spectra between different materials, and test optical properties for new materials.
With `PyTASER`, we provide a simple yet powerful method to simulate TAS for pristine materials, using data obtained from first principles calculations. This will not only assist experimentalists in comparing their data with theoretical estimates, but also encourage the mapping of a material's electronic structure with its optical properties, providing a wider understanding for newly discovered materials. 

PyTASER uses the principle of allowed vertical transitions between material bands to identify the possible electronic excitations that can occur in the respective ground 'dark' and excited 'light' stages. It does this by finding the Effective Absorption for each state; a product of the material's Joint Density of States (JDOS) and the transition probability for each band transition. These are determined from pre-calculated density functional theory calculations - once calculated, `PyTASER` then compares the resulting possible electronic excitations between the dark and light states. 

JDOS is defined as the density of allowed vertical band-band transitions based on the energetics and occupancy of each band - formalised by the equation defined in `Eqn. X`. This is fairly trivial to resolve for dark states, as the fully occupied valence bands and unoccupied conduction bands of such states lead to well defined fermi energies. It has seen common use in packages such as `AbiPy` and `OptaDOS`.

[Eqn X - JDOS dark equation]

However, determining the JDOS for the light state is much more difficult, as the initial monochromatic excitation leads to partial occupancies in both the valence and conduction bands. The resulting intra-band transitions lead to numerous fermi levels, which complicates the equation shown in `Eqn. X `.
`PyTASER` overcomes this by using quasi-fermi levels within the bands to treat the intra-band and inter-band transitions separately (`Eqn. XI`), with the partial occupancies estimated using the Fermi-Dirac approximation. This further introduces an element of control regarding both the temperature and the concentration of free carriers present in the material. The latter can be considered inversely analagous to the pump-probe time delay of experimental TAS. 

[Eqn. XI - JDOS equation used in PyTASER]

Alongside this, `PyTASER` also computes the transition probability from the `WAVEDER` file (when run with Vienna Ab-initio Simulation Package (VASP) electronic structure calculation code). This defines the probability of an optical transition actually taking place, and is found by computing the difference between the orbital derivatives found in `WAVEDER` to calculate the transition dipole moment. The transition probability can then be computed using `Eqn. XII`, and is multiplied into the JDOS to plot the effective absorption for each band transition. 

[Eqn XII - equation relating the transition probability with the orbital derivatives and transition dipole moment]

By directly comparing between the 'light' and 'dark' effective absorption values as shown in  `Fig. Z`, `PyTASER` can offer a full TAS profile with time-evolution for a range of many pristine materials, showing excellent qualitative accuracy compared to literature-reported spectra. 

[Fig. Z - Aron's schematic showing how light and dark transitions result in the overall TAS profile.]


# Statement of need

`PyTASER` is a Python package for simulating the TAS profiles of pristine materials in different temperature and time-delay conditions, using pre-calculated electronic structure data calculated from density functional theory. This makes it a powerful tool to not only scope out optical properties of new materials, but also to easily identify the cause of specific optical properties in difficult-to-interpret. While other software is being developed to also simulate TAS spectra (ref OPTADOS), they rely on performing additional complex density-functional-perturbation-theory (DFPT) calculations, and are difficult to interpret to unfamiliar users. To the best of our knowledge, `PyTASER` is unique in providing a full, detailed and qualitatively accurate TAS spectrum using pre-calculated DFT data, with the ability to account for complex optical processes such as stimulated emission, ground state bleaching and photo-induced absorption. 

`PyTASER` is designed to be usable by users with little experience in computational modelling methods and optical techniques, enabled by the following features:
- Use of Python as the programming language due to its low entry-barrier, flexibility and popularity in the materials modelling field.
- Full documentation and easy-to-follow workflows, with good unit-test coverage.
- Ability to produce publication ready figures, with a high level of manoeuvrability in plotting.
- Interfaced with the popular materials analysis package `pymatgen`.
- Compatible with VASP, the most popular electronic structure calculation code. Other codes can also be added with small amount of coding effort. 

A notable feature in `PyTASER` is the ability to plot a 'decomposed' TAS spectrum, plotting the individual band transitions alongside the overall spectrum (Fig. ZI). As a result, users can determine the exact band contributions towards different spectral features, which greatly simplifies analysis when identifying the cause of different optical processes. In experiment, this is a difficult process, requiring numerous stages of Fourier analysis and a deep understanding of quantum mechanics - `PyTASER` greatly delimits this technical barrier, providing a wider understanding of the material being investigated. In addition to this, `PyTASER` can also plot a comparison between the light and dark states to highlight the differences between them. 

[Fig. ZI - The decomposed TAS spectrum for CdTe ]

Furthermore,`PyTASER` is interfaced with the Materials Project (MP) database, containing over 100,000 relevant materials. This allows users to produce spectra for materials that they have not locally calculated and reduces the barrier for beginners to utilise the package for an immense range of different materials. It is important to note that the spectra calculated using this will be of a lower accuracy compared to locally calculated data, as the MP database does not include orbital derivative data - the resulting TAS spectrum will be independent of transition probabilities and rely solely on the JDOS. Despite this, the output will still provide a deep insight into the properties of the material. 


This code has been used by a number of users for a range of materials with differing electronic structures. It is currently being used in a study performed by the authors to identify common causes of recurring spectral features in various semiconducting materials.


# Acknowledgements

The authors are grateful for feature suggestion and testing from Anahita Manchala, as well as feedback from Liam Harnett-Caulfield. This work was supported by ...(financial acknowledgements)

# References
