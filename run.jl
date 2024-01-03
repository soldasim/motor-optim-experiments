include("boss.jl")
include("data.jl")

using Revise
includet("data/mle/mle_data.jl")
includet("data/bi/bi_data.jl")

# - - - RUN ONE OF THESE - - - - -
function run_mle()
    run_(; model_fitter_mode=:MLE, surrogate_mode=:GP)
end
function run_bi()
    run_(; model_fitter_mode=:BI, surrogate_mode=:GP)
end
# - - - - - - - -

function run_(;
    model_fitter_mode,
    surrogate_mode,
    parallel=true,
    info=true,
)
    X, Y = load_experiment_data(model_fitter_mode; info)
    problem = get_problem(X, Y; surrogate_mode)
    x, acq = test_script(problem; model_fitter_mode, surrogate_mode, parallel, info)
    
    info && @info "Saving data ..."
    save_iter_data(x, acq, problem.data; model_fitter_mode)
    save_new_x(model_fitter_mode, x)
    return x
end
