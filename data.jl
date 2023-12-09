using JLD2

include("./data/mle/mle_data.jl")
include("./data/bi/bi_data.jl")

save_path(::Val{:MLE}) = "./data/mle"
save_path(::Val{:BI}) = "./data/bi"

function load_experiment_data(model_fitter_mode)
    data = experiment_data(Val(model_fitter_mode))
    X = data[:, 1:4]' |> collect
    Y = data[:, 5:6]' |> collect
    return X, Y
end

function save_iter_data(x, acq, data; model_fitter_mode)
    path = save_path(Val(model_fitter_mode)) * "/iters"
    idx = length(readdir(path))

    data_ = data_dict(x, acq, data)
    save(path * "/iter_$(idx).jld2", data_)
end

"""
Convert `BOSS.ExperimentDataPost` to a disctionary better suited for serialization into files.
"""
function data_dict(x, acq, data)
    return Dict(
        "x" => x,
        "acq" => acq,
        "data" => Dict(
            "X" => data.X,
            "Y" => data.Y,
            "θ" => data.θ,
            "length_scales" => data.length_scales,
            "noise_vars" => data.noise_vars,
        ),
    )
end
