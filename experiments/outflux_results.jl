using DataDrivenRoadJunctions
using DelimitedFiles
using Printf
using Plots
using BSON
using Distributions
using Interpolations

dshash = "53e0af7"
dset_selection = string(gitroot(), "/out/dset_selection_", dshash, ".jls")
apids = deserialize(dset_selection)["applic_ids"]
gr()

function vis_outflux_data(km, id)

    outfile = string(gitroot(), "/out/outflux_D", id, "_M", km)

    bfile = string(outfile, ".bson")

    if !isfile(bfile)
        println("Data file ", bfile, " not found.")
        return
    end

    prog = BSON.load(bfile)[:progression]
    t = BSON.load(bfile)[:t]
    outflux = BSON.load(bfile)[:outflux]

    N = length(prog)
    plots = []
    for k = 1:N
        s = prog[k]
        incoming = plot(collect(s.xc₁), 1000 * s.u₁, color = "red", label = "",
                        title = string("t=", s.t, " incoming roads"))
        plot!(incoming, s.xc₂, 1000 * s.u₂, color = "blue", linestyle = :dash, label = "")
        push!(plots, incoming)
        push!(plots, plot(s.xc₃, 1000 * s.u₃, color = "black", label = "",
                          title = string("outgoing road"))) 
        open(string(outfile, "_",  k, ".dat"), "w") do io
            writedlm(io, [s.xc₁ s.xc₂ s.xc₃ 1000*s.u₁ 1000*s.u₂ 1000*s.u₃])
        end
    end

    savefig(plot(plots..., layout = (N, 2), size = (1800, 1500)), string(outfile, "_progression.pdf"))
    
    savefig(plot(t, outflux, label="outflux"), string(outfile, ".pdf"))
end

    
for id in apids
    for km = 1:7
        ofile = string(gitroot(), "/out/outflux_D", id, "_M", km, ".bson")
        if isfile(ofile)
            println("Processing ", ofile, "...")
            vis_outflux_data(km, id)
        end
    end
end
