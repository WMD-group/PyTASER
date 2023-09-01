---
title: 'PyTASER: Simulating transient absorption spectroscopy (TAS) for crystals from first principles'

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
    orcid: 0009-0007-7128-3465
    equal-contrib: false
    affiliation: "1, 2" # (Multiple affiliations must be quoted)	
  - name: Seán R. Kavanagh
    orcid: 0000-0003-4577-9647
    equal-contrib: false
    affiliation: "1, 2" # (Multiple affiliations must be quoted)	
  - name: Young Won Woo
    equal-contrib: false
    affiliation: "1, 3"	
  - name: Lucas G. Verga
    equal-contrib: false
    affiliation: "1" # (Multiple affiliations must be quoted)		
  - name: Alex M. Ganose
    orcid: 0000-0002-4486-3321
    equal-contrib: false
    affiliation: "4" # (Multiple affiliations must be quoted)	
  - name: Aron Walsh
    orcid: 0000-0001-5460-7033
    corresponding: true
    affiliation: "1"
	
affiliations:
 - name: Thomas Young Centre and Department of Materials, Imperial College London, London, United Kingdom
   index: 1
 - name: Thomas Young Centre and Department of Chemistry, University College London, London, United Kingdom
   index: 2
 - name: Department of Materials Science and Engineering, Yonsei University, Seoul, Korea
   index: 3
 - name: Department of Chemistry, Imperial College London, London, United Kingdom
   index: 4

date: 1 September 2023

bibliography: paper.bib
---

# Summary

Differential absorption (DA) concerns the varying absorption of electromagnetic radiation by different substances, leading to distinctive spectral features used for identification or analysis. Transient absorption spectroscopy (TAS) is one such technique that optically probes short-lived excited states. 

PyTASER is a Python package for simulating TAS spectra for crystalline materials from first-principles electronic structure calculations. This facilitates mapping from the electronic band structure of the material to changes in its optical properties upon excitation. The code preferentially allows input of data from local density functional theory (DFT) calculations, but is also interfaced with the Materials Project database [@materials_project] to allow rapid or high-throughput predictions based on the pre-computed electronic structure. 

TAS is a powerful tool used to characterise the optical properties and time-evolution of electronically excited states in materials. This is achieved by comparing the differences in optical absorption spectra between the ground 'dark' state and a photo-excited 'light' state, with a variable time delay between the two measurements to incorporate the time evolution of the excited state [@porter:1997; @berera:2009]. \autoref{eq:tas_spectrum} shows how this can be calculated.

\begin{equation}\label{eq:tas_spectrum}
  \\Delta\alpha(\Delta t) =  \alpha_{light}(\Delta t) - \alpha_{dark}
\end{equation}
Here, $\alpha$ refers to the effective absorption of the material system at the specified light conditions. $\Delta t$ refers to the time-delay between initial photo-excitation and absorption measurement. 

The approach is widely used to understand microscopic processes in photochemical and electrochemical transformations, including complex phenomena such as electron trapping and carrier recombination [@clarke:2010; @kafizas:2016], which dictate performance across many optoelectronic applications, such as solar cells [@kavanagh2021rapid], photocatalysts [@pastor2022electronic] and LEDs [@huang2021perovskite]. 
The drawback of modern TAS is that the spectra are often difficult to interpret - especially for crystals, where specific valence and conduction band structures can give rise to complex features. As a result, it is difficult for experimentalists to test optical properties for materials, especially between different experimental setups.

With `PyTASER`, we provide a simple yet powerful method to simulate TAS for crystalline materials, using data obtained from first principles calculations. This will not only assist experimentalists in comparing their data with theoretical estimates, but also encourage the mapping of a material's electronic structure with its optical properties under excitation, providing a wider understanding of complex materials. 

PyTASER uses the principle of allowed vertical optical transitions between material bands to identify the possible excitations that can occur in the respective ground 'dark' and excited 'light' stages. It does this by calculating the effective absorption for each state; a product of the material's joint density of states (JDOS) and the transition probability for each band transition. These are based on post-processing of density functional theory calculations. Once calculated, `PyTASER` then compares the change in electronic transitions between the dark and light states. 

![Schematics of the ground and excited state electronic structures and optical profiles. The ground 'dark' state is at the top, showing full occupancy and unoccupancy (blue, orange) for the conduction and valence bands respectively. The excited `light` state shows partial occupancy in a similar plot at the bottom. The overall DA plot is displayed to the right, the difference between the dark and light effective absorption plots. \label{fig:figure1}](Fig1.pdf){width=100mm}

## JDOS method

JDOS is defined as the density of allowed vertical band-to-band transitions based on the energy separation and occupancy of each band, formalised by the equation defined in \autoref{eq:jdos_dark}. This is straightforward to resolve for dark states in insulating crystals, as the fully occupied valence bands and unoccupied conduction bands of such states lead to well-defined Fermi levels. It has seen common use in packages such as [AbiPy](https://github.com/abinit/abipy) [@abipy] and [OptaDOS](https://github.com/optados-developers/optados) [@optados].

\begin{equation}\label{eq:jdos_dark}
  \rho(\varepsilon)=\frac{2}{8 \pi^3} \sum_v \sum_c \int \delta\left[\varepsilon_{c, \boldsymbol{k}}-\varepsilon_{v, \boldsymbol{k}}-\varepsilon\right] d^3 \boldsymbol{k}
\end{equation}
Here, $c$ and $v$ refer to the conduction and valence bands respectively. $\varepsilon$ refers to the energy of the respective band at kpoint $\boldsymbol{k}$.

Determining the JDOS for the light state is more difficult, as the initial 'pump' excitation leads to partial occupancies in both the valence and conduction bands, which can contribute to additional *intra*-band optical transitions.
`PyTASER` uses quasi-Fermi levels [@nelson2003physics; @dresselhaus2001solid] within the bands to address both intra-band and inter-band transitions (\autoref{eq:jdos_pytaser}). The partial occupancies ($f_{i,k}$ and $f_{f,k}$ in \autoref{eq:jdos_pytaser}) within the bands can be estimated by using the Fermi-Dirac distribution [@zannoni1999quantization; dirac1926theory] centred at these quasi-Fermi levels, as the light excitation causes (rapid) thermalisation of the material, resulting in excess charge carriers (holes and electrons).
The use of Fermi-Dirac statistics introduces two variables; the effective temperature and concentration of free carriers in the material. The latter is related to the strength of the initial pump, as well as the pump-probe time delay. These can be used to understand the time-evolution of the material's excited state.
This approach introduces two variables; the effective temperature and concentration of free carriers in the material. The latter is related to the strength of the initial pump, as well as the pump-probe time delay.
%
\begin{equation}
\label{eq:jdos_pytaser}
  \rho\left(\varepsilon, \varepsilon_{F, h}, \varepsilon_{F, e}, T\right) = \frac{2}{8 \pi^3} \sum_i \sum_{f>i} \int \delta\left[\varepsilon_{f, \boldsymbol{k}}-\varepsilon_{i, \boldsymbol{k}}-\varepsilon\right] f_{i, \boldsymbol{k}}\left(1-f_{f, \boldsymbol{k}}\right) d^3 \boldsymbol{k}
\end{equation}
%
Here, the subscripts $\varepsilon_{F, h}$ and $\varepsilon_{F, e}$  refer to the quasi-hole and quasi-electron Fermi levels, respectively. The subscripts $i$ and $f$ refer to the initial and final band states. The $f$ variable is the occupancy at the respective band at k-point $\boldsymbol{k}$.

## Optics method

`PyTASER` can also compute the optical transition probability from the frequency-dependent dielectric functions. The transition dipole moment, \autoref{eq:transition_probability}, is used to estimate the effective absorption strength for each band-to-band transition. 
%
\begin{equation}
\label{eq:transition_probability}
  \lambda_{i,f}= \left[\left\langle \phi_{i} | \mu_{T} | \phi_{f} \right\rangle \right]^{2}
\end{equation}
%
This naturally takes into account dipole selection rules, for example, forbidden $g$ → $g$ transitions at the $\Gamma$ point in a centrosymmetric crystal. Spin selection rules are also enforced for open-shell materials. By directly comparing the 'no pump' and 'pump' optical absorption values, shown in \autoref{fig:figure1}, this approach can offer a more realistic, albeit more computationally expensive, TAS profile.

## Differential absorption 

Beyond TAS, we have also included a function to calculate a direct differential absorption spectrum, i.e. 
%
\begin{equation}
\label{eq:da}
  \Delta A(\lambda) = A_{final}(\lambda) - A_{initial}(\lambda)
\end{equation}
%
Here, A_{final} could represent any change in the environment of the system such as de-lithiation in a battery or a transition state in a catalytic cycle. One caveat we emphasise is that the reliability of the predicted spectra depends on the quality of the underlying optical absorption spectra. There will certainly be cases where the inclusion of excitonic, thermal, and/or relativistic effects is necessary.

# Statement of need

`PyTASER` is a Python package for simulating spectral profiles of crystals in different temperature and time-delay conditions. This makes it a powerful tool to not only predict the transient optical response of new materials, but also to assign the origin of specific optical features. `PyTASER` is unique in providing detailed transient absorption spectra using pre-calculated DFT data, with the ability to account for complex optical processes such as stimulated emission, ground state bleaching, and photo-induced absorption [@berera:2009; @kafizas:2016]. 

`PyTASER` is designed for users with moderate experience in computational methods and optical techniques, enabled by the following features:

- Use of Python as the programming language due to its low entry barrier, flexibility and popularity in the materials modelling field.
- Documentation and easy-to-follow workflows with unit-test coverage.
- Ability to produce publication-ready figures, with flexibility in plotting.
- Interfaced with the popular materials analysis package [`pymatgen`](https://pymatgen.org/) [@pymatgen]. 
- Currently compatible with [VASP](https://www.vasp.at/wiki/index.php/The_VASP_Manual), the most popular electronic structure calculation code [@vasp]. Other codes can also be added with a small amount of coding effort. 

A notable feature in `PyTASER` is the ability to decompose individual band transitions alongside the overall spectrum (\autoref{fig:figure2}). As a result, users can determine the band contributions to different spectral features. In experiments, this is a difficult process, requiring numerous stages of analysis and a deep understanding of the material's electronic structure [@wang:2015; @kafizas:2016]. By identifying the key bands, this feature greatly simplifies the process of identifying atomistic origins of electronic states, which is important for identifying the decay kinetics of such states.
In addition to the overall spectrum, `PyTASER` can plot individual contributions for the separate light and dark states, showcasing their respective contributing transitions. 

![A 'decomposed' TAS spectrum for GaAs, calculated using data from the Materials Project. This plot indicates that the '-1,0' (orange) and '-2,-1' (pink) band transitions contribute most towards the TAS spectrum in the bandgap (blue, dashed) region, in these conditions. \label{fig:figure2}](Fig2.png){width=100mm}

Furthermore,`PyTASER` is interfaced with the [Materials Project](https://next-gen.materialsproject.org/) database. This allows users to produce spectra for systems that they have not locally calculated and reduces the barrier for beginners to utilise the package (within the JDOS approximation). The code has been tested for a range of materials with differing electronic structures, using both locally calculated and Materials Project data. 

# Acknowledgements

The authors are grateful for feature suggestions and testing from Anahita Manchala, feedback from Liam Harnett-Caulfield, and advice from Artem Bakulin and James Durrant. This work receieved support from the Royal Society and the Leverhulme Trust. S.R.K. acknowledges the EPSRC Centre for Doctoral Training in the Advanced Characterisation of Materials (CDT-ACM) (EP/S023259/1) for funding a Ph.D. studentship. Via membership of the UK’s HEC Materials Chemistry Consortium, which is funded by the EPSRC (EP/L000202, EP/R029431, and EP/T022213), this work used the ARCHER2 UK National Supercomputing Service (www.archer2.ac.uk) and the UK Materials and Molecular Modelling (MMM) Hub (Young EP/T022213).

# References
