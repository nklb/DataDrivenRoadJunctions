using Flux

function neuralnormalize(data, n_shift::Vector{<:Real}, n_factor::Vector{<:Real})
    N = length(n_shift)
    tomatrixrows( [(data[k, :] .- n_shift[k]) ./ n_factor[k] for k=1:N] )
end


function neuraldenormalize(data, n_shift::Vector{<:Real}, n_factor::Vector{<:Real})
    N = length(n_shift)
    tomatrixrows( [ data[k, :] .* n_factor[k] .+ n_shift[k] for k=1:N] )
end


struct NeuralRS
    M::Vector{GreenshieldModel}
    NNmod::Chain
    n_shift::Vector{Real}
    n_factor::Vector{Real}
end


function flux(RS::NeuralRS, ρ_in, nrzdinput::Bool=false)
    if !nrzdinput
        ρ = ρ_in
    else
        ρ = neuraldenormalize(ρ_in, RS.n_shift[1:3], RS.n_factor[1:3])
    end

    f_in = [flux(RS.M[k], ρ[k, :]) for k=1:3]

    if !nrzdinput
        NNin = neuralnormalize([ρ_in; tomatrixrows(f_in,)], RS.n_shift, RS.n_factor)
    else
        NNin = [ρ_in; neuralnormalize(tomatrixrows(f_in),
                                      RS.n_shift[4:6], RS.n_factor[4:6])]
    end

    NNout = RS.NNmod(NNin)
    [NNout[1,:]'; (NNout[2,:]-NNout[1,:])'; NNout[2, :]']
end


function flux(RS::NeuralRS, D::MacroData, nrzdinput::Bool=false)
    flux(RS, shiftedrho(D), nrzdinput)
end


function flux(RS::NeuralRS, D::Vector{MacroData}, nrzdinput::Bool=false)
    flux(RS, shiftedrho(D), nrzdinput)
end


function dspenalty(RS::NeuralRS, ρ_in, nrzdinput::Bool, reg::Bool)
    if !nrzdinput
        ρ = ρ_in
    else
        ρ = neuraldenormalize(ρ_in, RS.n_shift[1:3], RS.n_factor[1:3])
    end
    dspenalty(RS, ρ, reg)
end

