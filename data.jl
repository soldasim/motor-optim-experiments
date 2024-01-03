using JLD2

const ROUND_DIGITS = 4

save_path(::Val{:MLE}) = "./data/mle"
save_path(::Val{:BI}) = "./data/bi"

data_path(::Val{:MLE}) = "./data/mle/mle_data.jl"
data_path(::Val{:BI}) = "./data/bi/bi_data.jl"

function load_experiment_data(model_fitter_mode; info=true)
    data = experiment_data(Val(model_fitter_mode))
    info && @info "Loaded $(size(data)[1]) data points."
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

function save_new_x(model_fitter_mode, x)
    path = data_path(Val(model_fitter_mode))
    
    data = read(path, String)
    idx = findlast(']', data)
    data = data[1:idx-1]
    x_str = join(string.(round.(x; digits=ROUND_DIGITS)), ' ')
    data *=  x_str * " \n]\n"

    open(path, "w") do file
        write(file, data)
    end
end
