include("boss.jl")
include("data.jl")

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
)
    X, Y = load_experiment_data(model_fitter_mode)
    problem = get_problem(X, Y; surrogate_mode)
    x, acq = test_script(problem; model_fitter_mode, surrogate_mode, parallel)
    
    save_iter_data(x, acq, problem.data; model_fitter_mode)
    return x
end
