using Plots

include("../data.jl")
include("../Surrogate_Q_volne_parametry.jl")

@warn "Check which data files are loaded!"
run = "data_1"
include("../"*run*"/mle/mle_data.jl")
include("../"*run*"/bi/bi_data.jl")

function plot_runs(; init_data=4)
    @warn "Check that `init_data = $init_data` is correct!"
    X_mle, Y_mle = load_experiment_data(:MLE)
    X_bi, Y_bi = load_experiment_data(:BI)

    bfs_mle, fs_mle = bsf_series(X_mle, Y_mle; init_data)
    bfs_bi, fs_bi = bsf_series(X_bi, Y_bi; init_data)

    p = plot(; xlabel="iter", ylabel="Tav")
    iters = 0:length(bfs_mle)-1
    mle_color = :blue
    bi_color = :orange
    plot!(p, iters, bfs_mle; label="mle", color=mle_color)
    plot!(p, iters, bfs_bi; label="bi", color=bi_color)
    scatter!(p, iters, fs_mle; label="", color=mle_color)
    scatter!(p, iters, fs_bi; label="", color=bi_color)
    return p
end

function bsf_series(X, Y; init_data)
    @assert size(X)[2] == size(Y)[2]

    bsf = [max(fitness.(eachcol(X[:,1:init_data]), eachcol(Y[:,1:init_data]))...)]
    fs = [max(fitness.(eachcol(X[:,1:init_data]), eachcol(Y[:,1:init_data]))...)]
    for i in (init_data+1):size(X)[2]
        f = fitness(X[:,i], Y[:,i])
        if f > last(bsf)
            push!(bsf, f)
        else
            push!(bsf, last(bsf))
        end
        push!(fs, f)
    end

    @assert length(bsf) == 1 + (size(X)[2] - init_data)
    return -bsf, -fs
end

function fitness(x, y)
    @assert in_domain(x) && feasible_x(x)
    feasible_y(y) || return -Inf
    return -y[2]
end

in_domain(x) = all(ModelParam.domain()[1] .<= x .<= ModelParam.domain()[2])
feasible_x(x) = all(ModelParam.check_feas(x...) .>= 0.)
feasible_y(y) = all(y .<= ModelParam.y_max())
