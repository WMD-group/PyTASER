Change log
==========

v2.1.1
~~~~~~
- Remove broken link from docs

v2.1.0
~~~~~~
- Add `bg` option to `TASGenerator.from_vasp_outputs()` to allow scissor shifting of the bandgap to match experiment.
- Implement vectorised band filtering functions, massive speedup in processing times.
- Implement memory sharing for large arrays within multiprocessing.

v2.0.1
~~~~~~
- Minor documentation and website updates.

v2.0.0
~~~~~~
- Reduce the sum over band transitions to the energy mesh (min(-6, -energy_max), max(6, energy_max)) to
  make parsing of VASP optics more efficient.
- Switch to semantic versioning.

v23.6.3
~~~~~~~
- Remake PyTASER MP example workflow
- Update PyPI releasing workflow
- Update minor MPRester API issue

v23.6.2
~~~~~~~
- Add multiprocessing functionality to dramatically speed up parsing of large kpoint sets / large numbers of
  electronic bands, and add example notebook and HPC script for this.

v23.6.1
~~~~~~~
- Add functionality to parse orbital derivatives from `VASP` optics calculations, to incorporate oscillator strengths
  in the prediction of TAS spectra.
- Add functionality to save and load `Tas` objects to `json` file, to save time on re-analysing later on.
- Update plotting to avoid overlapping lines, ensure colour-consistency when the number of transitions changes,
  remove redundant transitions from the legend, pass kwargs to `plt.legend()` etc.
- Add the `from_vasp_outputs()` class method for `TASGenerator()` to directly generate from `VASP` outputs (streamlines
  this process).
- Add `yaxis` = "jdos_diff", "tas_absorption_only" and "alpha" plotting options.
- Add tests and examples for all the above changes.
- `temp` and `conc` are passed through from the initial user specification, to avoid confusion if changing later in
  plotting functions and expecting TAS to be re-calculated.
- General clean-up, avoidance of unnecessary warnings, addition of error catches, docstrings expansion...
