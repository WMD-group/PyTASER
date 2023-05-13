Change log
==========

v0.2.0
------
- Add functionality to parse orbital derivatives from `VASP` optics calculations, to incorporate oscillator strengths
  in the prediction of TAS spectra.
- Add functionality to save and load `Tas` objects to `json` file, to save time on re-analysing later on.
- Update plotting to avoid overlapping lines, ensure colour-consistency when the number of transitions changes,
  remove redundant transitions from the legend, pass kwargs to `plt.legend()` etc.
- Add the `from_vasp_outputs()` class method for `TASGenerator()` to directly generate from `VASP` outputs.
- Add `yaxis` = "jdos_diff", "tas_absorption_only" and "alpha" plotting options.
- Add tests and examples for all the above changes.
- General clean-up, avoidance of unnecessary warnings, addition of error catches, docstrings expansion...
