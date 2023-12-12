using BOSS
using Distributions
using OptimizationOptimJL, PRIMA

include("Surrogate_Q_volne_parametry.jl")

function get_acquisition()
    BOSS.ExpectedImprovement(;
        cons_safe=true,
        ϵ_samples=1,  # Does nothing as the `LinFitness` is used.
    )
end

function get_domain()
    return BOSS.Domain(;
        bounds = ModelParam.domain(),
        discrete = ModelParam.discrete_dims(),
        cons = (x) -> [ModelParam.check_feas(x...)...],
    )
end

function get_priors(x_dim, y_dim)
    # θ: A1, A2, A3
    param_range = ModelParam.param_range()
    θ_priors = [truncated(Normal(mid, dif/3); lower=0.) for (mid, dif) in zip(mean(param_range), param_range[2]-param_range[1])]
    
    # x: nk, dk, Ds, Q
    domain = ModelParam.domain()
    length_scale_priors = fill(Product([truncated(Normal(0., dif/3); lower=0.) for dif in (domain[2][i]-domain[1][i] for i in 1:x_dim)]), y_dim)
    
    # y: dP, Tav
    # noise_var_priors = fill(Dirac(1e-8), y_dim)
    noise_var_priors = [truncated(Normal(0., 1.); lower=0.), truncated(Normal(0., 1.); lower=0.)]

    return θ_priors, length_scale_priors, noise_var_priors
end

function get_problem(X, Y;
    surrogate_mode,  # :Semipar, :GP, :Param
)
    x_dim, y_dim = size(X)[1], size(Y)[1]

    objective = missing  # On cluster, the ansys model was used as objective.
    domain = get_domain()
    θ_priors, length_scale_priors, noise_var_priors = get_priors(x_dim, y_dim)
    model = get_surrogate(Val(surrogate_mode), θ_priors, length_scale_priors)

    BOSS.OptimizationProblem(;
        fitness = BOSS.LinFitness([0., -1.]),
        f = objective,
        domain,
        y_max = ModelParam.y_max(),
        model,
        noise_var_priors,
        data = BOSS.ExperimentDataPrior(X, Y),
    )
end

function get_surrogate(::Val{:Param}, θ_priors, length_scale_priors)
    BOSS.NonlinModel(;
        predict = (x, θ) -> ModelParam.calc(x..., θ...),
        param_priors = θ_priors,
    )
end
function get_surrogate(::Val{:GP}, θ_priors, length_scale_priors)
    BOSS.Nonparametric(;
        length_scale_priors,
    )
end
function get_surrogate(::Val{:Semipar}, θ_priors, length_scale_priors)
    BOSS.Semiparametric(
        BOSS.NonlinModel(;
            predict = (x, θ) -> ModelParam.calc(x..., θ...),
            param_priors = θ_priors,
        ),
        BOSS.Nonparametric(;
            length_scale_priors,
        ),
    )
end

function get_model_fitter(::Val{:MLE}, surrogate_mode; parallel=true)
    BOSS.OptimizationMLE(;
        algorithm=NelderMead(),
        multistart=200,
        parallel,
        apply_softplus=true,
        softplus_params=get_softplus_params(surrogate_mode),
        x_tol=1e-3,
    )
end
function get_model_fitter(::Val{:BI}, surrogate_mode; parallel=true)
    BOSS.TuringBI(;
        sampler=BOSS.PG(20),
        warmup=400,
        samples_in_chain=10,
        chain_count=8,
        leap_size=5,
        parallel,
    )
end
function get_model_fitter(::Val{:Random}, surrogate_mode; parallel=true)
    BOSS.RandomMLE()
end

get_softplus_params(::Val{:Param}) = fill(true, 3)
get_softplus_params(::Val{:Semipar}) = fill(true, 3)
get_softplus_params(::Val{:GP}) = nothing

function get_acq_maximizer(::Val{:Random}; parallel=true)
    BOSS.RandomSelectAM()
end
function get_acq_maximizer(::Val{:optim}; parallel=true)
    # The old implementation of PRIMA from `NLopt.jl` was used on the cluster
    # due to the new one from `PRIMA.jl` causing `StackOverflowError`.
    BOSS.CobylaAM(PRIMA;
        multistart=200,
        parallel,
        rhoend=1e-3,
    )
end

function check_surrogate_mode(surrogate_mode, model)
    if (surrogate_mode == :Semipar)
        @assert model isa BOSS.Semiparametric
    elseif (surrogate_mode == :GP)
        @assert model isa BOSS.Nonparametric
    elseif (surrogate_mode == :Param)
        @assert model isa BOSS.NonlinModel
    else
        @assert false
    end
end

"""
Run BOSS on the motor problem.

# Keywords
- `model_fitter_mode`: Defines the model parameter estimation technique. Choose from `:MLE` and `:BI`.
- `acq_maximizer_mode`: Defines the acquisition function maximization technique. Choose from `:optim` and `:Random`.
- `surrogate_mode`: Defines the surrogate model used. Choose from `:Semipar` and `:GP`.
"""
function test_script(problem;
    model_fitter_mode,  # :MLE, :BI
    acq_maximizer_mode=:optim,  # :optim, :Random
    surrogate_mode,  # :Semipar, :GP, :Param
    parallel=true,
)
    if acq_maximizer_mode == :Random
        model_fitter_mode = :Random
    end
    check_surrogate_mode(surrogate_mode, problem.model)

    model_fitter = get_model_fitter(Val(model_fitter_mode), Val(surrogate_mode); parallel)
    acq_maximizer = get_acq_maximizer(Val(acq_maximizer_mode); parallel)
    acquisition = get_acquisition()

    options = BOSS.BossOptions(;
        info=true,
        debug=false,
    )

    x, acq = boss!(problem; model_fitter, acq_maximizer, acquisition, options)
    return x, acq
end
