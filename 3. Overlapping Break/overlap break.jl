#############################################################################
# VARIABLE INDEX
# verbose = if true, solution is saved as .txt
# rule = precedence matrix rule: m, mp, ms, mps
# instance = from 0 to 10
# warmstart = grka, single, or none
# sortmethod = kids or top
# timelimit = maximum allowed solver time in seconds
# tugs = k = number of tugs
# rins = Relaxation Induced Neighborhood Search (RINS) heuristic frequency (default = -1, off = 0)
# cuts = Global cut control (-1 = default/automatic, 0 = off, 1-3 = least to most aggressive)
#############################################################################


######################################## Packages #############################################

using DataFrames, XLSX, Gurobi, JuMP, Test, CSV, DataStructures, ProgressMeter, Dates, LightGraphs, PrettyTables
import MathOptInterface

######################################## Dual Cycling #############################################

function dualcycling(; verbose = true, verbose_ws = true, rule = "m", instance = "3", warmstart = "grka", sortmethod = "kids", 
                    timelimit = 20, rins = -1, cuts = 2,  tugs = 4, breakstart = 10, breaklength = 6)

    path_to_input = joinpath(@__DIR__, "instances/input/"*rule*"_layout"*instance*".xlsx")

    if warmstart == "grka"
        filename = string("overlapBreak_", rule, instance, "_", warmstart, "_", sortmethod)
    else
        filename = string("overlapBreak_", rule, instance, "_", warmstart)
    end
    

######################################## Parameters #############################################

    # Set of trailers to be discharged (S)
    unloadlist = DataFrame(XLSX.readtable(path_to_input, "position_sequence"))
    S = unloadlist[!, :position_onboard] 
        # Number of trailers to be unloaded:
	numu = size(unloadlist, 1)
    
    # Set of trailers to be loaded (Q)
    loadlist = DataFrame(XLSX.readtable(path_to_input, "position_sequence"))
    Q = loadlist[!, :position_onboard]
        # Number of trailers to be loaded:
	numl= size(loadlist, 1)  

    # Break start and break length
    bs = breakstart
    bl = breaklength

    # Set of time steps (T)
    upper_bound = numl+numu+bl+1
    T = [1:upper_bound;]
    
    # Set of trailer handling precedence pairs (P)
    # (i,i') \in P
    P_s = DataFrame(XLSX.readtable(path_to_input, "precedence_unloading"))
    pred_u = P_s[!, :pred] 
    suc_u = P_s[!, :suc]
    numpu = size(P_s,1)
    
    # (j,j') \in P
    P_q = DataFrame(XLSX.readtable(path_to_input, "precedence_loading")) 
    pred_l = P_q[!, :pred] 
    suc_l = P_q[!, :suc] 
    numpl = size(P_q,1)

    # Number of tugs k
    k = tugs
    
######################################### Model #################################

    cd("output")
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit"=>timelimit, "RINS"=>rins, "LogFile"=>filename, "Cuts"=>cuts))

#################################### Decision Variables ########################################## 

    # x[i,t] == 1 if trailer is unloaded from position i at timestamp t
    @variable(model, x[i in S, t in T], Bin, start = 0)
    
    # y[j,t] == 1 if a cargo is loaded onto position j at timestamp t
	@variable(model, y[j in Q, t in T], Bin, start = 0)

    # wsq_t = number of tugs travelling from ship to quey at t
    @variable(model, 0 <= wvq[t in T] <= 4, Int, start = 0)

    # wqs_t = number of tugs travelling from quey to ship at t
    @variable(model, 0 <= wqv[t in T] <= 4, Int, start = 0)

    # wqq_t = number of tugs idling on quey at t
    @variable(model, 0 <= wqq[t in T] <= 4, Int, start = 0)

    # (4) Lower bound for makespan
    u_lb = (numu+numl)/k+1

    # (5) Upper bound for makespan
    u_ub = upper_bound

    # Makespan
    @variable(model, u_lb <= u <= u_ub, Int)


#################################### Objective & Constraints ##########################################    

    # (1) Objective: minimize the makespan
    @objective(model, Min, u) 

    # (2) Latest finishing time for loading for each trailer position
    @constraint(model, c2[j in Q, t in T], t*y[j,t] <= u)

    # (3) Latest finishing time for discharging for each trailer position
    @constraint(model, c3[i in S, t in T], t*x[i,t] <= u)

    # (4) and (5) are defined with the decision variables

    # (6) Each trailer position is discharged exactly once
    @constraint(model, c6[i in S], sum(x[i,t] for t in T) == 1) 

    # (7) Each trailer position is loaded exactly once
    @constraint(model, c7[j in Q], sum(y[j,t] for t in T) == 1) 

    # (8) Enforces discharging-discharging operations
    @constraint(model, c8[p in 1:numpu], sum((t-1)*x[suc_u[p],t] for t in T) >= sum(t*x[pred_u[p],t] for t in T)) 

    # (9) Enforces loading-loading operations
    @constraint(model, c9[p in 1:numpl], sum((t-1)*y[suc_l[p],t] for t in T) >= sum(t*y[pred_l[p],t] for t in T))

    # (10) Enforces loading-discharging operations
    @constraint(model, c10[p in 1:numu], sum(t*y[S[p], t] for t in T) >= sum(t*x[S[p],t] for t in T))

    # (11) Ensures that tugs do not idle on ship
    @constraint(model, c11[t in T[2:end]], sum(x[i,t] for i in S) + wvq[t] == sum(y[j,t-1] for j in Q) + wqv[t-1])

    # (12) Ensures that there are always k tugs
    @constraint(model, c12[t in T], sum(x[i,t] for i in S) + wvq[t] + sum(y[j,t] for j in Q) + wqv[t] + wqq[t] == k) #13 actually 12

    # (13) Enforces maximum of k/2 tugs are discharging a trailer at once
    @constraint(model, c13[t in T], sum(x[i,t] for i in S) <= 2)

    # (14) Enforces maximum of k/2 tugs are loading a trailer at once
    @constraint(model, c14[t in T], sum(y[j,t] for j in Q) <= 2)

    # (15) Ensures all tugs are on the quey at the first time step
    @constraint(model, sum(y[j,1] for j in Q) + wqv[1] + wqq[1] == k)

    # (16) Ensures no tugs are discharging at the first time step
    for i in S
        JuMP.fix(x[i, 1], 0; force=true)
    end

    # (17) Ensures no tugs are travelling from ship to quey at first time step
    JuMP.fix(wvq[1], 0; force=true)

    # (18) Ensures that x[i,t] and y[j,t] are binary (not necessary since we defined them as such)

    # (19) Ensures that wvq, wqv, wqq, and u are integers (not necessary since we defined them as such)

    # (20) Ensures that half of the tugs remain idle at bs and bs+bl
    @constraint(model, c20[t in T[bs]], wqq[t] >= k/2)
    @constraint(model, c21[t in T[bs+bl]], wqq[t] >= k/2)

    # (20) Ensures that all tugs remain idle starting at bs until bs+bl
    @constraint(model, c22[t in T[bs+1:bs+bl-1]], wqq[t] == k)

#################################### Warm Start ##########################################    
   
    if warmstart == "single"
        include("3. Overlapping Break/overlap break single.jl")
        ws_objective = single_warmstart(path_to_input, filename, verbose_ws, upper_bound, tugs, T, S, P_s, P_q, bs, bl, x, y, wqq, wqv, wvq, u)
    end

    if warmstart == "grka"
        include("3. Overlapping Break/overlap break grka.jl")
        ws_objective = grka_warmstart(path_to_input, filename, verbose_ws, sortmethod, upper_bound, tugs, T, S, Q, P_s, P_q, bs, bl, x, y, wqq, wqv, wvq, u)
    end
    
#################################### Optimize ##########################################    
   
    JuMP.optimize!(model)
    obj = Int(JuMP.objective_value(model))

    #################################### Save Results ##########################################    
   
    save_results(model, rule , instance, warmstart, sortmethod, ws_objective, "overlap break", "results.csv")
    if verbose == true
        save_loading_plan(tugs, upper_bound, x, y, wvq, wqv, wqq, T, S, filename)
    end

    cd("..")
end

function save_results(model, rule, instance, warmstart_type, warmstart_sort, warmstart_obj, model_type, save_file)
    best_objective = round(JuMP.objective_value(model))
    gap = round(JuMP.relative_gap(model) * 100, digits=2)
    bound = round(JuMP.objective_bound(model))
    solve_time = round(JuMP.solve_time(model), digits=0)
    save_array = [model_type, instance, rule, warmstart_type, warmstart_sort, warmstart_obj, solve_time, best_objective, bound, gap]
    
    # Create a DataFrame row
    new_row = DataFrame(Model = [model_type], Instance = [instance], Rule = [rule], Warmstart = [warmstart_type], Warmstart_sort = [warmstart_sort], Warmstart_objective = [warmstart_obj], Solve_time = [solve_time], Best_objective = [best_objective], Best_bound = [bound], Gap = [gap])
    
    # Check if the file exists and is not empty
    if isfile(save_file) && filesize(save_file) > 0
        # Append the new row to the existing file
        CSV.write(save_file, new_row, append = true, header = false)
    else
        # Write the new row with headers
        CSV.write(save_file, new_row)
    end
end

function save_loading_plan(tugs, upper_bound, x, y, wvq, wqv, wqq, T, S, filename)
    header = ["Time", "V to Q", "Q to V", "IDLE", "", "", "", ""]
    data = Matrix(undef, 0, 2*tugs)
    for t in 1:(upper_bound)
        row = ["$(T[t])"]
        push!(row, "V to Q = $(JuMP.value(wvq[t]))")
        push!(row, "Q to V = $(JuMP.value(wqv[t]))")
        push!(row, "IDLE = $(JuMP.value(wqq[t]))")
        ul_actions = tugs
        for i in S
            if ((JuMP.value(x[i,t]))==1)
                push!(row, "unload $(i)")
                ul_actions = ul_actions - 1
            end 
            if ((JuMP.value(y[i,t]))==1)
                push!(row, "load $(i)")
                ul_actions = ul_actions - 1
            end                 
        end
        for j in 1:ul_actions
            push!(row, "")
        end
        data = vcat(data, permutedims(row))
    end
        open(filename*"_loading_plan.txt", "w") do f
        pretty_table(f, data; header=header, alignment=:l)
    end
end