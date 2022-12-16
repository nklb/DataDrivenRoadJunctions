using BSON
using Optim
using BlackBoxOptim
using DataDrivenRoadJunctions

dshash = "53e0af7"
D_train = loadmacrodata(dshash, :training)
M = convert(Vector{GreenshieldModel}, BSON.load(string(gitroot(), "/out/data_", dshash, ".bson"))[:FD])

ofile1 = string(gitroot(), "/models/FlowMax_model1_", gitshorthead(), ".bson")
if !isfile(ofile1)
    
    function objectiveerror1(x)
        E_f = fluxerror(FlowMaxRS(M, x[1]), D_train)
        E_f[4]
    end

    gloptimum = bboptimize(objectiveerror1; SearchRange = [(0, 0.8)],
                           MaxFuncEvals = 10_000)
    x = best_candidate(gloptimum)


    bson(ofile1, Dict([("model", FlowMaxRS(M, x[1])), ("optimum", gloptimum),
                       ("argmin", x), ("fitness", best_fitness(gloptimum)),
                       ("ds", dshash)]))

end

ofile2 = string(gitroot(), "/models/FlowMax_model2_", gitshorthead(), ".bson")

if !isfile(ofile2)
    function objectiveerror2(x)
        E_f = fluxerror(FlowMaxRS(varyFD(M, x[2:4]), x[1]), D_train)
        E_f[4]
    end

    println("Optimizing simplified second order model (varying vmax)...")
    gloptimum = bboptimize(objectiveerror2; SearchRange = [(0, 1), (0.2, 5),
                                                      (0.2, 5), (0.2, 5)],
                           MaxFuncEvals = 20_000)

    x = best_candidate(gloptimum)

    bson(ofile2, Dict([("model", FlowMaxRS(varyFD(M, x[2:4]), x[1])), ("optimum", gloptimum),
                       ("argmin", x), ("fitness", best_fitness(gloptimum)),
                       ("ds", dshash)]))
end
