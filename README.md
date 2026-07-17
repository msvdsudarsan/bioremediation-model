# Bioremediation Model

This repository accompanies our study of a singularly perturbed nonlinear
dynamical model for the microbial bioremediation of pharmaceutical pollutants
in soil. The model couples the decay of a pharmaceutical pollutant to the
growth of a soil microbial community through a Monod uptake term together with
a toxic feedback, and the analysis characterises its stability, bistability,
and bifurcation behaviour.

The code and figures collected here reproduce every numerical result reported
in the manuscript, so that the findings can be inspected and re-run
independently.

## Authors

- M. S. V. D. Sudarsan (corresponding author), Department of Mathematics,
  School of Sciences, Humanities and Management, Dr. RVR & NRI Institute of
  Technology (Deemed to be University), Vijayawada, Andhra Pradesh, India.
  ORCID: 0009-0001-2126-6428
- Prabhakara Reddy Deevi Reddy, Mathematics and Computing Skills Unit,
  Preparatory Studies Centre, University of Technology and Applied Sciences,
  Nizwa, Sultanate of Oman. ORCID: 0000-0001-9726-2944
- Anitha Deevi Reddy, Mathematics and Computing Skills Unit, Preparatory
  Studies Centre, University of Technology and Applied Sciences, Nizwa,
  Sultanate of Oman. ORCID: 0000-0002-4283-4789
- Y. Rajesh Yadav, Department of Mathematics, Sri Venkateswara University,
  Tirupati, Andhra Pradesh, India. ORCID: 0009-0000-1815-9608

## Repository structure

- `code/` — the MATLAB and Python scripts that generate the results.
  - `bioremediation.m` — the main MATLAB script that integrates the model and
    produces the figures.
  - `bioremediation_tests.m` — the accompanying consistency checks (for
    example, the closed-form threshold against the numerical solution).
  - `robustness_tests.py` — a Python check of the model's qualitative
    behaviour under parameter variation.
- `figures/` — the figures reported in the paper, in PDF form.
- `results/` — the numerical output and comparison tables referenced in the
  text.

## Requirements

- MATLAB (R2021a or later is sufficient) for the `.m` scripts.
- Python 3.9 or later for `robustness_tests.py`; see `requirements.txt` for the
  Python dependencies.

## Running the code

From the `code/` directory, run `bioremediation.m` in MATLAB to reproduce the
figures, and `bioremediation_tests.m` to reproduce the numerical checks. The
Python script can be run with `python robustness_tests.py`.

## Citation

If you use this code, please cite the accompanying paper. Citation metadata is
provided in `CITATION.cff`.

## License

This project is released under the MIT License; see `LICENSE` for the full
text.
