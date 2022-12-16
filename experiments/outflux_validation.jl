# to run this script navigate into the experiments folder and run the following command:
# > julia --project=.. -p n outflux_validation.jl
# n is to be replaced with the number of parallel processes

using Distributed

@everywhere using Pkg;
@everywhere Pkg.activate("..") 
@everywhere using DataDrivenRoadJunctions

@everywhere using BSON
@everywhere using Distributions
@everywhere using Interpolations
@everywhere using BlackBoxOptim

@everywhere begin
        
    dshash = "53e0af7"
    dset_selection = string(gitroot(), "/out/dset_selection_", dshash, ".jls")
    apids = deserialize(dset_selection)["applic_ids"]

    function model_outflux(km, id, dshash)
        series = "BB-global-incl-LaneIDs"
        xio, dist = xinout(series)
        
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

        vehicles = getvehicledset(series, id) 
        trange, e_in, m_in, out = inoutfluxes(series, vehicles)
        Fe = linear_interpolation(trange, e_in, extrapolation_bc=Line()) 
        Fm = linear_interpolation(trange, m_in, extrapolation_bc=Line()) 

        # relaxation speed and time instances for plotting
        λ = 15
        T = trange[end] + 5
        tsteps= LinRange(0, T, round(Integer, 4*T))
        oF = zeros(length(tsteps))
        # boundary condition, choose between ZeroFlux HomNeumann Outgoing Periodic
        bc = FixedFluxOutgoing

        # coupling model, choose between CentralRelaxationLimit and TrafficFlowMaximization
        mfile = model_files[km]
        cpl = BSON.load(mfile)["model"]

        # setup of the problem and the initial condition
        fcut = km in [5, 6]
        p = Problem21(f, g, h, λ, cpl, bc, false, Fe, Fm, fcut)
        s = constantinitialdata21(0, 0, 0, dist/2)

        # initialize subplots array and gr plotting framework
        plots = []
        instances = []

        ns = 10
        t_plot = 0
        
        outfile = string(gitroot(), "/out/outflux_D", id, "_M", km)
        # simulation and plotting loop
        for k = 1:length(tsteps)
            t = tsteps[k]
            s = centralcoupling(p, s, t, FirstOrder, 0.24)
            oF[k] = s.oF

            if k in 1:round(Integer, length(tsteps)/ns):length(tsteps)
                push!(instances, deepcopy(s))
            end
        end
        bson(string(outfile, ".bson"), progression = instances, t = tsteps, outflux = oF)
    end
end


@sync @distributed for id in apids
    for km = 1:7
        ofile = string(gitroot(), "/out/outflux_D", id, "_M", km, ".bson")
        if !isfile(ofile)
            model_outflux(km, id, dshash)
        end
    end
end
