using MAT
using MATLAB

function matimport(series::String, importtemplates::Bool=false)
    dsets = readdir(string(gitroot(), "/FCD/", series))
    templates = Dict{Integer, VehicleTemplate}()
    ids = Set{Integer}()

    for dset in dsets
        if dset[1] != '1' || dset[end-2:end] != "mat"
            continue
        end
        id = dset[1:end-4]
        println("Dataset: ", id)
        dpath = string(gitroot(), "/FCD/", series, "/", dset)
        eval_string(string("mx_veh=load('", dpath,"');"))
        L = matopen(dpath)
        mat_veh = read(L, "vehicles")
        N = length(mat_veh)
        vehicles = Vector{Vehicle}()
        for i = 1:N
            println("vehicle ", i)
            v = copy(mat_veh[i])
            setindex!(v, convert(Int64, v["vehicle_id"]), "vehicle_id")
            setindex!(v, convert(Int64, v["template_id"]), "template_id")
            # transform discrete_trajectory and contour table
            for table in ["discrete_trajectory", "contour"]
                eval_string(string("tdata=table2cell(mx_veh.vehicles{", string(i),
                                   "}.", table,");"))
                eval_string(string("theader=mx_veh.vehicles{", string(i),
                                   "}.", table,".Properties.VariableNames;"))
                @mget tdata
                @mget theader
                setindex!(v, DataFrame([theader[c]=>tdata[:,c]
                                        for c in 1:length(theader)]),
                          table)
            end
            # transform fitted functions
            if v["fitTypeE"] != "poly3" || v["fitTypeN"] != "poly3"
                error("FitType is no poly3 but ", v["fitTypeE"])
            end
            for f in ["traj_EtGlob", "traj_tEGlob", "traj_NtGlob", "traj_tNGlob"]
                eval_string(string("fx=mx_veh.vehicles{", string(i),
                                   "}.", f, ";"))
                eval_string(string("coef=[fx.p4; fx.p3; fx.p2; fx.p1];"))
                @mget coef
                setindex!(v, Polynomial(coef), f)
            end
            push!(vehicles, Vehicle(v["vehicle_id"], v["discrete_trajectory"],
                                    v["fitTypeE"], v["traj_EtGlob"], v["traj_tEGlob"], v["fitTypeN"],
                                    v["traj_NtGlob"], v["traj_tNGlob"], v["fitted_trajectoryGlob"],
                                    v["fitted_speedGlob"], v["fitted_accGlob"], v["template_id"]))
            if !in(ids, v["template_id"]) && importtemplates
                setindex!(templates, VehicleTemplate(v["Width"], v["Length"], v["contour"]), v["template_id"])
                push!(ids, v["template_id"])
            end
        end
        io = open(string(dpath[1:end-3], "jls.xz"), "w")
        io = TranscodingStream(XzCompressor(), io)
        serialize(io, vehicles)
        close(io)
    end
    if importtemplates
        io = open(string(gitroot(), "/FCD/", series, "/VehicleTemplates.jls.xz"), "w")
        io = TranscodingStream(XzCompressor(), io)
        serialize(io, sort(templates))
        close(io)
    end
end
