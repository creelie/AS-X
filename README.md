# Black Holes Carry Positive Discrete Information

Paper (`paper/main.tex`, `paper/main.pdf`) and the programs that check it.

- `code/twin_exact.py 21`: the twining genera from the fixed points, exactly in Q(ζ_N).
- `code/dmz_check.py 100`: the identity (B² − E₄A²)/(144A) = A₂,₁ − 2𝓗 up to q⁹⁹.
- `code/fitW.py`, `code/sturm.py`: the closed form W_g = f_g φ₋₂,₁ + h_g φ₀,₁ up to the Sturm bound (data in `code/Wfit.pkl`).
- `code/kser2.py 48200`: the coefficient condition for the fourteen comparison series, m ≤ 12045 (output in `code/kser2_48200.log`, about 12 minutes).
- `code/bound2.py`, `code/bound_check.py`: the analytic bound for m ≥ 11000, with 40-digit arithmetic.
- `code/scan.py`, `code/figdata.py`: S_a for D ≤ 600 (data in `code/C_600.pkl`).
- `paper/figures/make_figs.py`: the figures (run inside `paper/figures`).

Run the programs inside `code/`. Requires Python 3 with python-flint, sympy, mpmath, numpy and matplotlib.

Author: Deep Bhattacharjee. Licence: see `LICENSE`.

DOI: https://doi.org/10.5281/zenodo.24123192 (all versions)
