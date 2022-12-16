using BSON
using Optim
using BlackBoxOptim
using DataDrivenRoadJunctions


dshash = "53e0af7"
D_train = loadmacrodata(dshash, :training)
M = convert(Vector{GreenshieldModel}, BSON.load(string(gitroot(), "/out/data_", dshash, ".bson"))[:FD])

ofile1 = string(gitroot(), "/models/LinReg_model1_", gitshorthead(), ".bson")
if !isfile(ofile1)

    function objectiveerror1(x)
        RS = LinearRS(M, vec(x))
        E = fluxerror(RS, D_train)
        
        E[4] + dspenalty(RS, shiftedrho(D_train), true) +
            negpenalty(RS, shiftedrho(D_train), true)
    end

    println("Optimizing linear first order model...")
    srange = [(-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5),
              (-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5)]
    gloptimum = bboptimize(objectiveerror1; SearchRange = srange,
                           MaxFuncEvals = 150_000)
    x = best_candidate(gloptimum)


    bson(ofile1, Dict([("model", LinearRS(M, vec(x))), ("optimum", gloptimum),
                       ("argmin", x), ("fitness", best_fitness(gloptimum)),
                       ("ds", dshash)]))
end

ofile2 = string(gitroot(), "/models/LinReg_model2_", gitshorthead(), ".bson")

if !isfile(ofile2)
    function objectiveerror2(x)
        RS = LinearRS(varyFD(M, x[13:15]), vec(x[1:12]))
        E = fluxerror(RS, D_train)
        
        E[4] + dspenalty(RS, shiftedrho(D_train), true) +
            negpenalty(RS, shiftedrho(D_train), true) +
            fdpenalty(RS, D_train)
    end

    println("Optimizing linear second order model (varying vmax)...")
    srange = [(-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5),
              (-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5), (-5, 5),
              (0.2, 5), (0.2, 5), (0.2, 5)]
    gloptimum = bboptimize(objectiveerror2; SearchRange = srange,
                           MaxFuncEvals = 150_000)
    x = best_candidate(gloptimum)

    bson(ofile2, Dict([("model",  LinearRS(varyFD(M, x[13:15]), vec(x[1:12]))), ("optimum", gloptimum),
                       ("argmin", x), ("fitness", best_fitness(gloptimum)),
                       ("ds", dshash)]))
end
