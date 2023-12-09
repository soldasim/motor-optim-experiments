
# How to run the scripts:

## First-time:
- Evaluate the 4 random initial data points in ANSYS and fill them into '/data/mle/mle_data.jl' and '/data/bi/bi_data.jl'. (Both are the same.)

## Initialize:
- Open terminal at the root of this repo.
- Start julia: `julia`
- Activate the environment: `] activate .`
- Load the scripts: `include("run.jl")`

## Workflow:
- Run `run_mle()` or `run_bi()` in julia to get new `x`.
- Evalute ANSYS at `x` to get new `y`.
- Add the concatenated data point `xy` as a new row in '/data/mle/mle_data.jl' or '/data/bi/bi_data.jl'.
- Repeat ...
