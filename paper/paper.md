---
title: 'PyTASER: A Python package to simulate transient absorption spectroscopy (TAS) for crystals from first principles'
tags:
  - Python
  - materials science
  - first-principles
  - optics
  - differential-absorption 
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
date: 23 July 2023
bibliography: paper.bib
---

# Summary

Differential absorption concerns the varying absorption of electromagnetic radiation by different substances, leading to distinctive spectral features used for identification or analysis. Transient absorption spectroscopy (TAS) is one such technique that optically probes short-lived excited states. 

PyTASER is a Python package for simulating TAS spectra for crystalline materials from first-principles electronic structure calculations. This facilitates mapping between the band structure of the material to changes in its optical properties. The code preferentially allows input of data from local density functional theory (DFT) calculations, but is also interfaced with the Materials Project database to allow rapid or high-throughput predictions. 

TAS is a powerful tool used to characterise the properties and time-evolution of electronically excited states of materials. This is achieved by comparing the optical absorption spectra between the ground 'dark' state, and a monochromatic-laser excited 'light' state, with a variable time delay between the two measurements to incorporate time evolution of the excited state [@porter:1997; @berera:2009]. The approach is widely used to understand microscopic processes in photochemical and electrochemical transformations, including phenomena such as electron trapping and carrier recombination [@clarke:2010; @kafizas:2016], which are important properties when considering technological application of a material. The drawback of modern TAS is that the spectra are difficult to interpret. This is especially the case for crystals, where specific valence and conduction band structures can give rise to complex features. As a result, it is difficult for experimentalists to compare spectra between different materials, and test optical properties for new materials.

With `PyTASER`, we provide a simple yet powerful method to simulate TAS for pristine materials, using data obtained from first principles calculations. This will not only assist experimentalists in comparing their data with theoretical estimates, but also encourage the mapping of a material's electronic structure with its optical properties, providing a wider understanding for complex materials. 

PyTASER uses the principle of allowed vertical transitions between material bands to identify the possible electronic excitations that can occur in the respective ground 'dark' and excited 'light' stages. It does this by finding the effective absorption for each state; a product of the material's joint density of states (JDOS) and the transition probability for each band transition. These are determined from pre-calculated density functional theory calculations. Once calculated, `PyTASER` then compares the resulting possible electronic excitations between the dark and light states. 

![Schematics of the ground and excited electronic structures and optical profiles. The ground 'dark' state is at the top, showing full occupancy and unoccupancy(blue and orange, respectively) bands. The excited `light` state shows partial occupancy in a similar plot at the bottom. The overall DA plot is displayed to the right, the difference between the effective absorption plots for light and dark. \label{fig:figure1}](Fig1_rev.pdf){width=100mm}

## JDOS method

JDOS is defined as the density of allowed vertical band-band transitions based on the energetics and occupancy of each band, formalised by the equation defined in `Equation [@jdos_dark]}`. This is straightforward to resolve for dark states, as the fully occupied valence bands and unoccupied conduction bands of such states lead to well defined Fermi energies. It has seen common use in packages such as [AbiPy](https://github.com/abinit/abipy) [@abipy] and [OptaDOS](https://github.com/optados-developers/optados) [@optados].

::: {jdos_dark .math}
  \rho(\varepsilon)=\frac{2}{8 \pi^3} \sum_v \sum_c \int \delta\left[\varepsilon_{c, \boldsymbol{k}}-\varepsilon_{v, \boldsymbol{k}}-\varepsilon\right] d^3 \boldsymbol{k}
:::
Here, $c$ and $v$ refer to the valence and conduction bands. $\varepsilon$ refers to the energy of the respective band at kpoint $\boldsymbol{k}$.

However, determining the JDOS for the light state is more difficult, as the initial monochromatic excitation leads to partial occupancies in both the valence and conduction bands. The resulting intra-band transitions lead to numerous fermi levels, which complicates the equation shown in `Equation [@jdos_dark]`.
`PyTASER` overcomes this by using quasi-fermi levels [@katahara:2014; @reddy:2016] within the bands to treat the intra-band and inter-band transitions separately (`Equation [@jdos_pytaser]`), with the partial occupancies estimated using the Fermi-Dirac approximation [@zitter:1987]. This further introduces an element of control regarding both the temperature and the concentration of free carriers present in the material. The latter can be considered inversely analogous to the pump-probe time delay of experimental TAS. 

::: {jdos_pytaser .math}
  \begin{gathered}
  \rho\left(\varepsilon, \varepsilon_{F, h}, \varepsilon_{F, e}, T\right) \\
  =\frac{2}{8 \pi^3} \sum_i \sum_{f>i} \int \delta\left[\varepsilon_{f, \boldsymbol{k}}-\varepsilon_{i, \boldsymbol{k}}-\varepsilon\right] f_{i, \boldsymbol{k}}\left(1-f_{f, \boldsymbol{k}}\right) d^3 \boldsymbol{k}
  \end{gathered}
:::
In this equation, the subscripts $\varepsilon_{F, h}$ and $\varepsilon_{F, e}$  refer to the quasi-hole and quasi-electron fermi levels, respectively. The subscripts $i$ and $f$ refer to the initial and final band states. The $f$ variable is the occupancy at the respective band, at kpoint $\boldsymbol{k}$.

## Optics method

Alongside the JDOS method, `PyTASER` also computes the optical transition probability from the frequency dependent dielectric tensors. The band-to-band transition dipole moment can be computed using `Equation [@transition_probability] `, and is multiplied into the JDOS to estimate the effective absorption for each band transition.  By directly comparing between the 'light' and 'dark' optical absorption values as shown in \autoref{fig:figure1}, `PyTASER` can offer a more realistic, albeit more computationally expensive, TAS profile with good agreement compared to literature-reported spectra. 

::: {transition_probability .math}
  \lambda_{i,f}= \left[\left\langle \phi_{i} | \mu_{T} | \phi_{f} \right\rangle \right]^{2}
:::

[Eqn 3 - equation relating the transition probability with the orbital derivatives and transition dipole moment]

# Statement of need

`PyTASER` is a Python package for simulating TAS profiles of crystals in different temperature and time-delay conditions. This makes it a powerful tool to not only scope out optical response of new materials, but also to assign the cause of specific optical features. `PyTASER` is unique in providing a detailed TAS spectrum using pre-calculated DFT data, with the ability to account for complex optical processes such as stimulated emission, ground state bleaching, and photo-induced absorption [@berera:2009; @kafizas:2016]. 

`PyTASER` is designed to be usable by users with moderate experience in computational methods and optical techniques, enabled by the following features:

- Use of Python as the programming language due to its low entry-barrier, flexibility and popularity in the materials modelling field.
- Full documentation and easy-to-follow workflows, with good unit-test coverage.
- Ability to produce publication ready figures, with a high level of manoeuvrability in plotting.
- Interfaced with the popular materials analysis package [`pymatgen`](https://pymatgen.org/) [@pymatgen]. 
- Currently compatible with [VASP](https://www.vasp.at/wiki/index.php/The_VASP_Manual), the most popular electronic structure calculation code [@vasp]. Other codes can also be added with a small amount of coding effort. 

A notable feature in `PyTASER` is the ability to plot a 'decomposed' TAS spectrum, plotting the individual band transitions alongside the overall spectrum (\autoref{fig:figure2}). As a result, users can determine the exact band contributions towards different spectral features, which greatly simplifies analysis when identifying the cause of different optical processes. In experiment, this is a difficult process, requiring numerous stages of Fourier analysis and a deep understanding of quantum mechanics [@wang:2015; @kafizas:2016] - `PyTASER` greatly delimits this technical barrier, providing a wider understanding of the material being investigated. In addition to the overall spectrum, `PyTASER` can also plot the individual contributions for the separate light and dark states, showcasing their respectively contributing bands. 

![A 'decomposed' TAS spectrum for GaAs, calculated using data from the Materials Project. This plot indicates that the '-1,0' (orange) and '-2,-1' (pink) band transitions contribute most towards the TAS spectrum in the bandgap (blue, dashed) region, in these conditions. \label{fig:figure2}](decomposed_tas_gaas.png){width=100mm}

Furthermore,`PyTASER` is interfaced with the [Materials Project](https://next-gen.materialsproject.org/) database [@materials_project]. This allows users to produce spectra for systems that they have not locally calculated and reduces the barrier for beginners to utilise the package for an immense range of different materials (within the JDOS approximation). The code has been tested for a range of materials with differing electronic structures, using both locally calculated and Materials Project data. 

# Acknowledgements

The authors are grateful for feature suggestion and testing from Anahita Manchala, feedback from Liam Harnett-Caulfield, and advice from Artem Bakulin and James Durrant. This work was supported by the Royal Society and the Leverhulme Trust.

# References
