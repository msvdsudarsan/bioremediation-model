# Bioremediation Model

Code and numerical output for the paper:

> **Two–Timescale Dynamics of Microbial Degradation of Pharmaceutical Pollutants in Soil: Bistability, Tipping Thresholds, and Bifurcation Analysis**
> Sri Venkata Durga Sudarsan Madhyannapu, Prabhakara Reddy Deevi Reddy, Anitha Deevi Reddy, Y. Rajesh Yadav.

The model couples a pharmaceutical pollutant `P` with the microbial biomass `M` that degrades it. Microbial uptake follows a saturating (Monod) response and the pollutant exerts a toxic feedback on the degraders. Because pollutant kinetics are fast relative to population turnover, the system is analysed as a singularly perturbed problem. The analysis and simulations reveal a bistable regime (recovery vs. washout) bounded by a saddle-node bifurcation in the toxicity coefficient.

## Repository layout

```
code/       MATLAB and Python scripts
figures/    figures reproduced by the scripts (PDF)
results/    saved console output and comparison tables
```

## Code

| File | Purpose |
|------|---------|
| `code/bioremediation.m` | Equilibria, stability, phase portrait, time series and the toxicity bifurcation (Figures 1-3). |
| `code/bioremediation_tests.m` | Supplementary Tests A-C: singular-perturbation convergence, closed-form fold boundary, and basin of attraction. |
| `code/bioremediation_hill_test.m` | Sensitivity to the toxic functional form: bilinear `gamma*P*M` vs. saturating (Hill) `gamma*P*M/(K_T+P)`. Reports endpoint states, recovery-basin fractions, and the washout-state eigenvalue under each law. |
| `code/robustness_tests.py` | Independent Python cross-check of the equilibria and thresholds. |

The MATLAB scripts run without modification in MATLAB Online and in GNU Octave.

## Key numerical findings

- Nominal parameters: `sigma=1.0, delta=0.1, alpha=0.8, K=1.0, theta=0.6, mu=0.2, gamma=0.05`.
- Three equilibria: washout `E0=(10,0)` (stable), saddle `E2=(6.37,0.52)`, recovery `E1=(0.63,3.04)` (stable).
- Closed-form tipping toxicity `gamma_c = (sqrt(theta)-sqrt(mu))^2 / K = 0.10718` at `K=1`; agreement with the numerically located fold is within `5e-6`.
- Recovery basin covers `89.1%` of the sampled window under the bilinear toxic term.
- **Toxic-form sensitivity:** with a saturating (Hill) toxic term (`K_T=1`) the washout state loses stability (transverse eigenvalue `-0.1545` -> `+0.30`), every sampled trajectory recovers, and the recovery basin grows to `97.1%`. The bilinear model is therefore a conservative, worst-case description of remediation failure.

## Reproducing the results

MATLAB / Octave:

```matlab
run('code/bioremediation.m')
run('code/bioremediation_tests.m')
run('code/bioremediation_hill_test.m')
```

Python cross-check:

```bash
pip install -r code/requirements.txt
python code/robustness_tests.py
```

Saved console output for each run is in `results/`.

## Citation

If you use this code, please cite the paper. Metadata is provided in `CITATION.cff`.

## License

Released under the MIT License. See `LICENSE`.
