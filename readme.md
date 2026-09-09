# Sackin Index on Full k-ary Trees

This repository computes the distribution of the Sackin index on unlabelled, unordered (non-plane), full k-ary rooted trees. It provides exact integer multiplicities, frequency-weighted Beta approximations, diagnostic tables, and figures. This results are used in the paper: **Sackin Index Distributions on Full \(k\)-ary Trees: Spectrum, Multiplicities and Beta-Type Behaviour**

## Repository contents

| File or directory | Purpose |
|---|---|
| `sackin_analysis_2.py` | Exact-integer implementation and command-line interface documented here |
| `sackin_analysis.py` | Same implementation as sackin_analysis_2.py but with float arithmetic to accelerate computation |
| `requirements.txt` | required packages |
| `results_2/` | outputs of sackin_analysis_2.py used in the paper |


## Installation

From the repository directory, install python if necessary, create an environment and install the dependencies:

```bash
uv python install 3.14
uv venv --python 3.14
uv pip install -r requirements.txt
```

Check that the script runs

```bash
uv run sackin_analysis_2.py --help
```


## Quick start

```bash
# Print the exact distribution for binary trees with 10 leaves.
uv run sackin_analysis_2.py frequencies --n 10

# Fit a Beta distribution and print the diagnostics.
uv run sackin_analysis_2.py fit --n 20

# Save a density comparison figure.
uv run sackin_analysis_2.py plot --n 20 --output results_2/sackin20.pdf
```

## Computing $A_{n,k}(s)$

The multiplicity $A_{n,k}(s)$ is the number of unlabelled, unordered full $k$-ary rooted tree shapes with $n$ leaves and Sackin index $s$. Every internal node has exactly $k$ children. 


### 1. Base case and recursive decomposition

The tree consisting of a single leaf has Sackin index zero. Therefore,

$$
A_{1,k}(0)=1.
$$

For a larger tree, let its root subtrees be $T_1,\ldots,T_k$, with leaf counts $n_1,\ldots,n_k$. They satisfy

$$
n_1+\cdots+n_k=n.
$$

Attaching the subtrees to a new root increases the depth of every leaf by one. Consequently,

$$
S(T)=n+\sum_{j=1}^{k}S(T_j).
$$

Thus, combining subtree distributions requires adding their Sackin indices and then shifting the result by $n$.

### 2. Enumerating possible subtree sizes

A full $k$-ary tree with $n$ leaves has

$$
p=\frac{n-1}{k-1}
$$

internal nodes. After removing the root, its subtrees contain a total of $p-1$ internal nodes.

The code enumerates all tuples

$$
0\le q_1\le\cdots\le q_k,
\qquad
q_1+\cdots+q_k=p-1.
$$

Here, $q_j$ is the number of internal nodes in the $j$-th subtree, whose leaf count is

$$
n_j=1+(k-1)q_j.
$$

Requiring the tuple to be nondecreasing ensures that each unordered collection of subtree sizes is considered once. However, repeated sizes still require special treatment: several different subtree shapes may have the same number of leaves.

### 3. Representing distributions as polynomials

For each admissible leaf count $j$, the code stores the polynomial

$$
C_{j,k}(z)=\sum_{q\ge0}A_{j,k}\!\left(S_{\min}(j,k)+q\right)z^q.
$$

The coefficient at position $q$ counts shapes whose Sackin index exceeds the minimum by $q$. Storing relative increments avoids keeping the initial zero coefficients of an absolute-index polynomial.

Polynomial multiplication combines independent choices of subtrees. If

$$
P(z)=\sum_a p_a z^a,
\qquad
Q(z)=\sum_b q_b z^b,
$$

then

$$
[z^t]P(z)Q(z)=\sum_{a+b=t}p_aq_b.
$$

This is precisely the number of pairs of choices with total increment \(t\). The code performs this multiplication using integer convolution.

### 4. Handling repeated subtree sizes

Suppose the root has $r$ children with the same leaf count $j$. 

Simply raising $C_{j,k}(z)$ to the power $r$ would count ordered selections and would therefore overcount.

Instead, define $H_{j,r}(z)$ to count unordered multisets of $r$ shapes of size $j$, grouped by their total relative Sackin increment. The code computes these polynomials using

$$
H_{j,0}(z)=1,
$$

$$
mH_{j,m}(z)
=
\sum_{i=1}^{m}
C_{j,k}(z^i)\,H_{j,m-i}(z).
$$

This is Newton’s identity for complete homogeneous symmetric functions. The substitution $z\mapsto z^i$ multiplies every increment by $i$, while leaving the number of available shapes at each increment unchanged.


### 5. Combining the groups and restoring the Sackin offset

For a fixed collection of root-subtree sizes, let $r_j$ be the number of children with $j$ leaves. The code multiplies the corresponding multiset polynomials:

$$
G(z)=\prod_{j:\,r_j>0}H_{j,r_j}(z).
$$

If the coefficient of $z^q$ in $G(z)$ is $c_q$, it contributes $c_q$ tree shapes at the absolute Sackin index

$$
s=n+\sum_j r_jS_{\min}(j,k)+q.
$$

The contribution is added to the parent’s coefficient array at position

$$
s-S_{\min}(n,k).
$$

Summing these contributions over all admissible collections of subtree sizes gives $A_{n,k}(s)$.

There is no additional ordering factor: each rooted unordered tree is determined by the multiset of its root subtrees, and the construction already counts that multiset once.

### 6. Specialized binary recurrence

For $k=2$, the code uses a simpler implementation. To state its recurrence, write the absolute-index generating polynomial as

$$
F_n(z)=\sum_s A_{n,2}(s)z^s.
$$

Starting from $F_1(z)=1$, the recurrence for $n\ge2$ is

$$
F_n(z)=z^n
\left(
\sum_{i=1}^{\lfloor(n-1)/2\rfloor}F_i(z)F_{n-i}(z)
+
\begin{cases}
\displaystyle
\frac{F_{n/2}(z)^2+F_{n/2}(z^2)}{2},
& n\text{ even},\\[6pt]
0,
& n\text{ odd}.
\end{cases}
\right).
$$

For unequal subtree sizes, ordinary multiplication counts each pair once because the sizes distinguish the two subtrees. For equal sizes, the multiset formula removes the overcounting caused by exchanging them. The factor $z^n$ accounts for the increase in leaf depths under the new root.



### 7. Exact arithmetic and reuse of intermediate results

In `sackin_analysis_2.py`, multiplicities are stored as Python arbitrary-precision integers. Polynomial convolution uses exact integer arithmetic, and every division in the multiset recurrence is checked for a zero remainder.

Previously computed subtree distributions and multiset polynomials are cached and reused. Floating-point arithmetic is introduced when converting counts to probabilities and when calculating moments or fitting Beta distributions.


### Reference

Macdonald, I. G. (1995). *Symmetric Functions and Hall Polynomials*, 2nd ed. Clarendon Press, Oxford. 

## Command-line help

Relative output paths are resolved from the current working directory.

```bash
uv run sackin_analysis_2.py --help
uv run sackin_analysis_2.py fit --help
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
uv run sackin_analysis_2.py frequencies --n 10
uv run sackin_analysis_2.py frequencies --n 21 --k 3 --output results/ternary_frequencies.csv
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
uv run sackin_analysis_2.py fit --n 20 --k 2
uv run sackin_analysis_2.py fit --n 21 --k 3 --cdf ternary_cdf.pdf
```

| Option | Required? | Default | Meaning |
|---|---|---|---|
| `--n` | Yes | — | Number of leaves |
| `--k` | No | `2` | Branching degree |
| `--cdf` | No | No plot | PDF/PNG destination for a cumulative-distribution comparison |

Results are printed as JSON. There is no `--output` option for this command. To save the JSON, use shell redirection:

```bash
uv run sackin_analysis_2.py fit --n 20 > fit_n20.json
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
uv run sackin_analysis_2.py plot --n 20 --k 2 --output results/sackin20.pdf
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
uv run sackin_analysis_2.py binary-table --output results/binary.csv
```

| Option | Required? | Default |
|---|---|---|
| `--output` | No | `binary_beta_table.csv` |

Fits binary trees with `n = 10, 15, ..., 100`, saves a CSV and prints the table. Columns include `n`, fitted parameters, exact and Beta moments, moment errors, and `D`.

The parser provides no option to change this leaf-count range.

## 5. higher-table — fixed internal-node count, varying k

```bash
uv run sackin_analysis_2.py higher-table --p 10 --k-min 3 --k-max 5 --output results/higher.csv
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
uv run sackin_analysis_2.py grid --ps 10 20 --k-min 3 --k-max 5 --output results/grid.csv --max-output results/max_errors.csv
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
uv run sackin_analysis_2.py paper --output-dir results
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

## Global compatibility option

`--fft-threshold INTEGER` defaults to `1000000000`. It is accepted for compatibility but ignored: only used in sackin_analysis.py

If supplied, it belongs before the command:

```bash
uv run sackin_analysis.py --fft-threshold 1000000 frequencies --n 10
```
