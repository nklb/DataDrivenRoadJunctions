# this needs to be defined in the main name space, otherwise the NN model won't be loaded
function softplus(x::Real)
    log(exp(x) + 1)
end

using DataDrivenRoadJunctions
using Plots
using BSON
using BlackBoxOptim
using Distributions
using DelimitedFiles
using Printf

dshash = "53e0af7"
D_train = loadmacrodata(dshash, :training)
D_test = loadmacrodata(dshash, :test)
M = loadmacrodata(dshash, :FD)

C_models = [string(gitroot(), "/models/FlowMax_model1_8a8c781.bson"),
            string(gitroot(), "/models/FlowMax_model2_8a8c781.bson")]

D_models = [string(gitroot(), "/models/LinReg_model1_8a8c781.bson"),
            string(gitroot(), "/models/LinReg_model2_8a8c781.bson"),
            string(gitroot(), "/models/neural_model_8a8c781.bson")]

model_files = vcat(C_models, D_models)

merror(mfile, D) = fluxerror(BSON.load(mfile)["model"], D)

D = D_test

for k=1:4
    print(" | ")
    for j =1:length(model_files)
        @printf("%s | ", tolatexfloat(merror(model_files[j], D)[k]))
    end
    print(" \n")
end

L1 = BSON.load(D_models[1])["model"]
L2 = BSON.load(D_models[2])["model"]

for A in [L1.A, L2.A]
    for k =1:3
        print(" | ")
        for j =1:6
            @printf("%1.3f | ", A[k, j])
        end
        print(" \n")
    end
    
end

ofile = string(gitroot(), "/out/Example_flx", gitshorthead(), "_")

D = MacroData("BB-global-incl-LaneIDs", "10314")
FxD = shiftedflux(D)
Nt = round(Int, size(FxD, 2)/2)
FM = plot(D.t[1:Nt], FxD[1, 1:Nt], label="f_E")
plot!(FM, D.t[1:Nt], FxD[2, 1:Nt], label="f_T")
plot!(FM, D.t[1:Nt], FxD[3, 1:Nt], label="f_O")

open(string(ofile, "data.dat"), "w") do io
    writedlm(io, [D.t[1:Nt] FxD[1, 1:Nt] FxD[2, 1:Nt] FxD[3, 1:Nt]])
end

for k=1:length(model_files)
    mfile = model_files[k]
    FxM = flux(BSON.load(mfile)["model"], D)

    open(string(ofile, "M", k, ".dat"), "w") do io
        writedlm(io, [D.t[1:Nt] FxM[1, 1:Nt] FxM[2, 1:Nt] FxM[3, 1:Nt]])
    end
end
