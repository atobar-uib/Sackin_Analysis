# Sackin Index on Full k-ary Trees

This repository computes the distribution of the Sackin index on unlabelled, unordered (non-plane), full k-ary rooted trees. It provides exact integer multiplicities, frequency-weighted Beta approximations, diagnostic tables, and figures. This results are used in the paper: **Sackin Index Distributions on Full \(k\)-ary Trees: Spectrum, Multiplicities and Beta-Type Behaviour**

## Repository contents

| File or directory | Purpose |
|---|---|
| `sackin_analysis.py` | Exact-integer implementation and command-line interface documented here |
| `check_andres.py` | Independent total-count for unlabelled, unordered (non-plane), full k-ary rooted trees.   |
| `compare_count.py` | check that the sum of obtained values coincides with the total number of trees |
| `requirements.txt` | required packages |
| `results_2/` | outputs of sackin_analysis.py used in the paper |


## Installation

From the repository directory, install python if necessary, create an environment and install the dependencies:

```bash
uv python install 3.14
uv venv --python 3.14
uv pip install -r requirements.txt
```

Check that the script runs

```bash
uv run sackin_analysis.py --help
```


## Quick start

```bash
# Print the exact distribution for binary trees with 10 leaves.
uv run sackin_analysis.py frequencies --n 10

# Fit a Beta distribution and print the diagnostics.
uv run sackin_analysis.py fit --n 20

# Save a density comparison figure.
uv run sackin_analysis.py plot --n 20 --output results_2/sackin20.pdf
```

## Command-line help

Relative output paths are resolved from the current working directory.

```bash
uv run sackin_analysis.py --help
uv run sackin_analysis.py fit --help
```

Every command supports `-h` or `--help`.

## What n, k and p mean

| Parameter | Meaning |
|---|---|
| `n` | Number of leaves |
| `k` | Number of children of each internal node |
| `p` | Number of internal nodes |


## Commands at a glance

| Command | Purpose |
|---|---|
| `frequencies` | Compute counts and probabilities for one `(n, k)` |
| `fit` | Fit a Beta distribution and print diagnostics |
| `plot` | Save a distribution plot with the fitted Beta density |
| `binary-table` | Fit binary trees for `n = 10, 15, ..., 100` |
| `higher-table` | Compare branching degrees at one internal-node count |
| `grid` | Compare multiple internal-node counts and branching degrees |
| `paper` | Generate the configured manuscript tables, regressions and figures |

## 1. frequencies — exact counts and probabilities

```bash
uv run sackin_analysis.py frequencies --n 10
uv run sackin_analysis.py frequencies --n 21 --k 3 --output results/ternary_frequencies.csv
```

| Option | Required? | Default | Meaning |
|---|---|---|---|
| `--n` | Yes | — | Number of leaves |
| `--k` | No | `2` | Branching degree |
| `--output` | No | No file | CSV destination |

Without `--output`, the table is printed. With it, the table is saved as CSV, and missing parent directories are created.

Columns:

- `Sackin`: Sackin-index value.
- `absolute_frequency`: number of tree shapes with that value.
- `relative_frequency`: absolute frequency divided by the total number of shapes.

## 2. fit — Beta parameters and fit diagnostics

```bash
uv run sackin_analysis.py fit --n 20 --k 2
uv run sackin_analysis.py fit --n 21 --k 3 --cdf ternary_cdf.pdf
```

| Option | Required? | Default | Meaning |
|---|---|---|---|
| `--n` | Yes | — | Number of leaves |
| `--k` | No | `2` | Branching degree |
| `--cdf` | No | No plot | PDF/PNG destination for a cumulative-distribution comparison |

Results are printed as JSON. There is no `--output` option for this command. To save the JSON, use shell redirection:

```bash
uv run sackin_analysis.py fit --n 20 > fit_n20.json
```

If `--cdf` includes a directory, create that directory first: this command does not create it automatically.

### Interpreting the results

The Sackin index is normalized as `X = (S - S_min) / (S_max - S_min)`. The Beta fit uses the conditional distribution given `0 < X < 1`, excluding both endpoints.

| Field | Meaning |
|---|---|
| `n`, `k`, `internal_nodes` | Tree-size parameters |
| `alpha`, `beta` | Fitted Beta parameters |
| `mu_exact`, `var_exact` | Mean and variance of the conditional interior distribution, excluding endpoints |
| `mu_beta`, `var_beta` | Mean and variance of the fitted Beta distribution |
| `error_mu`, `error_var` | Absolute differences between those moments |
| `D` | Kolmogorov distance between the conditional interior CDF and the Beta CDF |
| `D_full` | Corresponding distance for the full distribution |
| `x_max`, `side_max` | Location of the largest conditional CDF difference and side of the discrete jump |
| `extreme_mass` | Combined probability at the two endpoints |
| `optimizer_success`, `optimizer_message` | Optimization status |

Smaller CDF distances indicate closer agreement. The moment errors compare conditional interior moments with the fitted Beta moments. Counts use exact integers; probabilities, moments and fitting use floating-point arithmetic. The code returns diagnostics even if the optimizer reports failure, so inspect its status.

## 3. plot — distribution and Beta density

```bash
uv run sackin_analysis.py plot --n 20 --k 2 --output results/sackin20.pdf
```

| Option | Required? | Default | Meaning |
|---|---|---|---|
| `--n` | Yes | — | Number of leaves |
| `--k` | No | `2` | Branching degree |
| `--output` | Yes | — | Figure destination, such as `.pdf` or `.png` |

Saves a figure and prints the fit diagnostics as JSON. Missing parent directories are created.

The black curve shows the full exact probabilities divided by the normalized spacing, so its heights can be compared with a density. The red curve is the Beta density fitted to the interior. For a cumulative-distribution plot instead, use `fit --cdf`.

## 4. binary-table — fixed binary size range

```bash
uv run sackin_analysis.py binary-table --output results/binary.csv
```

| Option | Required? | Default |
|---|---|---|
| `--output` | No | `binary_beta_table.csv` |

Fits binary trees with `n = 10, 15, ..., 100`, saves a CSV and prints the table. Columns include `n`, fitted parameters, exact and Beta moments, moment errors, and `D`.

The parser provides no option to change this leaf-count range.

## 5. higher-table — fixed internal-node count, varying k

```bash
uv run sackin_analysis.py higher-table --p 10 --k-min 3 --k-max 5 --output results/higher.csv
```

| Option | Required? | Default | Meaning |
|---|---|---|---|
| `--p` | No | `20` | Internal nodes |
| `--k-min` | No | `3` | First branching degree |
| `--k-max` | No | `10` | Last branching degree, included |
| `--output` | No | `higher_k_table.csv` | CSV destination |

It saves and prints a table containing `p`, `k`, `n`, fitted parameters, moments, errors and `D`.

## 6. grid — several p and k values

```bash
uv run sackin_analysis.py grid --ps 10 20 --k-min 3 --k-max 5 --output results/grid.csv --max-output results/max_errors.csv
```

| Option | Required? | Default | Meaning |
|---|---|---|---|
| `--ps` | No | `10 20 30 40` | One or more internal-node counts, separated by spaces |
| `--k-min` | No | `3` | First branching degree |
| `--k-max` | No | `10` | Last branching degree, included |
| `--output` | No | `higher_k_grid.csv` | All fitted cases |
| `--max-output` | No | `higher_k_max_errors.csv` | Maximum-error summary |

Computes every combination of the supplied `p` values and the inclusive `k` range. The example produces six cases.

For each `p`, the summary reports `max_D`, `max_error_mu` and `max_error_var` over the selected branching degrees. These maxima may come from different `k` values. Both CSV files are saved, and the summary is printed.

## 7. paper — complete configured manuscript output

```bash
uv run sackin_analysis.py paper --output-dir results
```

| Option | Required? | Default |
|---|---|---|
| `--output-dir` | No | `results` |

Creates the destination directory and generates:

- `binary_beta_table.csv`: binary fits for `n = 10, 15, ..., 100`.
- `binary_power_law_regressions.json`: log-log regressions for the fitted binary parameters.
- `Sackin10.pdf`, `Sackin15.pdf`, `Sackin20.pdf`, `Sackin25.pdf`: binary distribution figures.
- `higher_k_p20_table.csv`: fits for `p = 20`, `k = 3, ..., 10`.
- `Sackin41_3.pdf`, `Sackin61_4.pdf`, `Sackin81_5.pdf`, `Sackin101_6.pdf`: higher-degree distribution figures.
- `higher_k_grid.csv`: fits for `p = 10, 20, 30, 40` and `k = 3, ..., 10`.
- `higher_k_max_errors.csv`: maximum errors for each `p`.
- `table_binary.tex`, `table_higher_k_p20.tex`, `table_max_errors.tex`: LaTeX table fragments.

These ranges are fixed in the code. Exact enumeration can be expensive, particularly for the full grid and `paper` command.


## Checking total tree counts

`compare_counts.py` compares the sum of Sackin multiplicities with a separately implemented Andres-based counting recurrence in `check_andres.py`. Comparisons use exact integers and raise an error if a mismatch occurs.

```bash
# Default: k=2,...,10 and p=0,...,15.
uv run compare_counts.py

# Select branching degrees and the maximum internal-node count.
uv run compare_counts.py --ks 3 4 5 --p-max 20

# Display both counts for every case.
uv run compare_counts.py --ks 4 --p-max 10 --verbose
```

- `--ks`: branching degrees to check.
- `--p-max`: maximum number of internal nodes, inclusive. The corresponding leaf count is `n = 1 + (k - 1) * p`.
- `--verbose`: print each comparison.

