# this needs to be defined in the main name space, otherwise the NN model won't be loaded
function softplus(x::Real)
    log(exp(x) + 1)
end

using Plots
using BSON
using BlackBoxOptim
using Distributions
using DelimitedFiles
using DataDrivenRoadJunctions

gr()
series = "BB-global-incl-LaneIDs"
xio, dist = xinout(series)

dshash = "53e0af7"
M = loadmacrodata(dshash, :FD)

C_models = [string(gitroot(), "/models/FlowMax_model1_8a8c781.bson"),
            string(gitroot(), "/models/FlowMax_model2_8a8c781.bson")]

D_models = [string(gitroot(), "/models/LinReg_model1_8a8c781.bson"),
            string(gitroot(), "/models/LinReg_model2_8a8c781.bson"),
            string(gitroot(), "/models/neural_model_8a8c781.bson")]

model_files = vcat(C_models, D_models)

# setup of flux functions on incoming edges and outgoing edge
f(u::Real) = flux(M[1], u)
g(u::Real) = flux(M[2], u)
h(u::Real) = flux(M[3], u)

# relaxation speed and time instances for plotting
λ = 15
tsteps= LinRange(0, 10, 4)
bc = HomNeumann

for km = 1:length(model_files)
    mfile = model_files[km]
    cpl = BSON.load(mfile)["model"]

    # setup of the problem and the initial condition
    fcut = km in [3, 4]
    p = Problem21(f, g, h, λ, cpl, bc, false, false, false, fcut)
    global s = constantinitialdata21(0.4 * M[1].rhomax, 0.5 * M[2].rhomax,
                                     0.8 * M[3].rhomax, dist/2)

    # initialize subplots array and gr plotting framework
    plots = []
    outfile = string(gitroot(), "/out/prediction_M", km)
    # simulation and plotting loop
    for k=1:length(tsteps)
        t = tsteps[k]
        global s = centralcoupling(p, s, t, FirstOrder, 0.24)
        incoming = plot(collect(s.xc₁), 1000 * s.u₁, color = "red", label = "",
                        title = string("t=", t, " incoming roads"),
                        ylims = (0, 1000 * max(M[1].rhomax, M[2].rhomax)))
        plot!(incoming, s.xc₂, 1000 * s.u₂, color = "blue", linestyle = :dash, label = "")
        push!(plots, incoming)
        push!(plots, plot(s.xc₃, 1000 * s.u₃, color = "black", label = "",
                          title = string("outgoing road"), ylims = (0, 1000 * M[3].rhomax)))
        open(string(outfile, "_",  k, ".dat"), "w") do io
            writedlm(io, [s.xc₁ s.xc₂ s.xc₃ 1000*s.u₁ 1000*s.u₂ 1000*s.u₃])
        end
    end

    savefig(plot(plots..., layout = (4, 2), size = (1800, 1500)), string(outfile, ".pdf"))

end
