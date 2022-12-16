module DataDrivenRoadJunctions

include("vehicles.jl")
include("data.jl")
include("dataimport.jl")
include("models.jl")
include("neural_networks.jl")
include("CentralNetworkScheme.jl")

export GreenshieldModel, MacroData, FlowMaxRS, LinearRS, NeuralRS
export loadmacrodata, flux, velocity, shiftedrho, shiftedv, shiftedflux
export fluxerror, Chain, Dense, rawshifteddata, rmsd, AMSGrad, matimport
export gitroot, gitshorthead, tolatexfloat, NNlib, σ, softplus, Flux
export neuraldenormalize, neuralnormalize, dspenalty, rms, rmssp, Vehicle
export serialize, deserialize, getdsetids, getvehicledset, observationtime
export getvolumes, numberdensities, defaulttrange, inoutfluxes, xinout
export ZeroFlux, HomNeumann, Outgoing, Periodic, Problem21, centralcoupling
export constantinitialdata21, FirstOrder, SecondOrder, stats, FixedFluxOutgoing
export varyFD, negpenalty, fdpenalty

function tolatexfloat(x)
    s = @sprintf("%1.3e", x)
    m = s[1:5]
    sg = s[7]
    ep = parse(Int, s[8:end])
    out = string("\$", m)
    if ep != 0
        out = string(out, " \\times 10^{")
        if sg == '-'
            out = string(out, "-")
        end
        out = string(out, ep, "}")
    end
    out = string(out, "\$")
end

end
