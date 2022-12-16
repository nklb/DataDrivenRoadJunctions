using BSON
using Random
using DataDrivenRoadJunctions


dshash = "53e0af7"
D_train = loadmacrodata(dshash, :training)
D_test = loadmacrodata(dshash, :test)
M = convert(Vector{GreenshieldModel}, BSON.load(string(gitroot(), "/out/data_", dshash, ".bson"))[:FD])

n_shift = vec(mean(rawshifteddata(D_train), dims=2))
n_factor = vec(std(rawshifteddata(D_train), dims=2))

x_train = neuralnormalize(shiftedrho(D_train), n_shift[1:3], n_factor[1:3])
y_train = shiftedflux(D_train)

x_test = neuralnormalize(shiftedrho(D_test), n_shift[1:3], n_factor[1:3])
y_test = shiftedflux(D_test)

model =  Chain(Dense(6 => 12, σ), Dense(12 => 2, softplus)) 

function loss(x, y)
    RS = NeuralRS(M, model, n_shift, n_factor)
    f_out = flux(RS, x, true)

    # penalty term for negative fluxes
    npenalty = rmssp(-vec(f_out))
    penalty = dspenalty(RS, x, true, true) +  npenalty

    rmsd(vec(f_out), vec(y)) +  penalty
end

function pureloss(x, y, model)
    RS = NeuralRS(M, model, n_shift, n_factor)
    f_out = flux(RS, x, true)
    rmsd(vec(f_out), vec(y))
end

opt = AMSGrad(1e-3)
train_data = [(x_train, y_train)]
parameters = Flux.params(model)

N_train = size(x_train, 2)
idx = collect(1:N_train)

println("initial training loss: ", loss(x_train, y_train), ", initial test loss: ",
        loss(x_test, y_test))


for epoch in 1:60
    shuffle!(idx)
    for i in idx[1:round(Integer, N_train/2)]
        (x, y) = (x_train[:, i], y_train[:, i])
        gs = Flux.gradient(parameters) do
            loss(x,y)
        end
        Flux.Optimise.update!(opt, parameters, gs)
    end
    tloss = loss(x_train, y_train)
    println("training epoch: ", epoch, ", training loss: ", loss(x_train, y_train),
            ", test loss: ", loss(x_test, y_test), ", training flux error:", pureloss(x_train, y_train, model),
            ", test flux loss: ", pureloss(x_test, y_test, model))
end


ofile = string(gitroot(), "/models/neural_model_", gitshorthead(), ".bson")

bson(ofile, Dict([("model", NeuralRS(M, model, n_shift, n_factor)), ("x_train", x_train),
                  ("x_test", x_test), ("etrain", pureloss(x_train, y_train, model)),
                  ("etest", pureloss(x_test, y_test, model),  ("ds", dshash))]))

