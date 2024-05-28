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
# presolve = Default is -1, 2 is maybe better?
# mipfocus = default is 0, maybe 2 is better
#############################################################################


######################################## Packages #############################################

using DataFrames, XLSX, Gurobi, JuMP, Test, CSV
using DataStructures
using ProgressMeter, Dates, LightGraphs
using PrettyTables
import MathOptInterface

######################################## Dual Cycling #############################################

function dualcycling(; verbose = true, verbose_ws = true, rule = "mp", instance = "0", warmstart = "grka", sortmethod = "kids", 
                    timelimit = 20, rins = -1, cuts = 2, tugs = 4, break_length = 4)

    path_to_input = joinpath(@__DIR__, "instances/input/"*rule*"_layout"*instance*".xlsx")

    if warmstart == "grka"
        filename = string("longBreak_", rule, instance, "_", warmstart, "_", sortmethod)
    else
        filename = string("longBreak_", rule, instance, "_", warmstart)
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
	numl = size(loadlist, 1)  

    # Break length
    b = break_length

    # Calculate upper bound on number of time steps
    upper_bound = numl+numu+b+1

    # Set of time steps (T)
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

    # Set of tugs (K)
    K = [1:tugs;]
    
######################################### Model #################################

    cd("output")
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit"=>timelimit, "RINS"=>rins, "LogFile"=>filename, "Cuts"=>cuts))

#################################### Decision Variables ##########################################    
    
    # x[i,t] == 1 if trailer is unloaded from position i at timestamp t
    @variable(model, x[i in S, t in T, k in K], Bin, start = 0)

    # y[j,t] == 1 if a cargo is loaded onto position j at timestamp t
    @variable(model, y[j in Q, t in T, k in K], Bin, start = 0)

    # wvq_tk = number of tugs travelling from ship to quey at t
    @variable(model, wvq[t in T, k in K], Bin, start = 0)

    # wqv_tk = number of tugs travelling from quey to ship at t
    @variable(model, wqv[t in T, k in K], Bin, start = 0)

    # wqq_tk = number of tugs idling on quey at t
    @variable(model, wqq[t in T, k in K], Bin, start = 0)

    # wqb_tk = 1 if tug k is on breeak at t
    @variable(model, wqb[t in T, k in K], Bin, start = 0)

    # bs_tk = 1 if break for tug k starts at t
    @variable(model, sb[t in T, k in K], Bin, start = 0)

    # (4) Lower bound for makespan
    u_lb = ((numu+numl)/tugs+1) + b + 1

    # (5) Upper bound for makespan
    u_ub = upper_bound

    # Makespan
    @variable(model, u_lb <= u <= u_ub, Int)

#################################### Objective & Constraints ##########################################    

    # (1) Objective: minimize the makespan
    @objective(model, Min, u) 

    # (2) Latest finishing time for loading for each trailer position
    @constraint(model, c2[k in K, j in Q, t in T], t*y[j,t,k] <= u)

    # (3) Latest finishing time for discharging for each trailer position
    @constraint(model, c3[k in K, i in S, t in T], t*x[i,t,k] <= u)

    # (4) and (5) are defined with the decision variables

    # (6) Each trailer position is discharged exactly once
    @constraint(model, c6[i in S], sum(sum(x[i,t,k] for t in T) for k in K) == 1) 

    # (7) Each trailer position is loaded exactly once
    @constraint(model, c7[j in Q], sum(sum(y[j,t,k] for t in T) for k in K) == 1) 

    # (8) Enforces discharging-discharging operations
    @constraint(model, c8[p in 1:numpu], sum(sum((t-1)*x[suc_u[p],t,k] for t in T) for k in K) >= sum(sum(t*x[pred_u[p],t,k] for t in T) for k in K))

    # (9) Enforces loading-loading operations
    @constraint(model, c9[p in 1:numpl], sum(sum((t-1)*y[suc_l[p],t,k] for t in T) for k in K) >= sum(sum(t*y[pred_l[p],t,k] for t in T) for k in K))

    # (10) Enforces loading-discharging operations
    @constraint(model, c10[p in 1:numu], sum(sum(t*y[S[p],t,k] for t in T) for k in K) >= sum(sum(t*x[S[p],t,k] for t in T) for k in K))

    # (11) Ensures that tugs do not idle on ship
    @constraint(model, c12[k in K, t in T[2:end]], sum(x[i,t,k] for i in S) + wvq[t,k] == sum(y[j,t-1,k] for j in Q) + wqv[t-1,k])

    # (12) Ensures that there are always k tugs
    @constraint(model, c13[t in T, k in K], (sum(x[i,t,k] for i in S) + wvq[t,k] + sum(y[j,t,k] for j in Q) + wqv[t,k] + wqq[t,k] + wqb[t,k]) == 1)

    # (13) Enforces maximum of k/2 tugs are discharging a trailer at once
    @constraint(model, c14[t in T], sum(sum(x[i,t,k] for i in S) for k in K) <= tugs/2)

    # (14) Enforces maximum of k/2 tugs are loading a trailer at once
    @constraint(model, c15[t in T], sum(sum(y[j,t,k] for j in Q) for k in K) <= tugs/2)

    # (15) Ensures all tugs are on the quey at the first time step
    @constraint(model, c16[k in K], (sum(y[j,1,k] for j in Q) + wqv[1,k] + wqq[1,k]) == 1) 

    # (16) Ensures no tugs are discharging at the first time step
    for i in S
        for k in K
            JuMP.fix(x[i,1,k], 0; force=true)
        end
    end

    # (17) Ensures no tugs are travelling from ship to quey at first time step
    for k in K
        JuMP.fix(wvq[1,k], 0; force=true)
    end

    # (18) Ensures that x[i,t] and y[j,t] are binary (not necessary since we defined them as such) (also ensure that sb is binary, and also all the rest (wqq, wqv, etc.)))

    # (19) Ensures that wvq, wqv, wqq, and u are integers (not necessary since we defined them as such)

    # (20) Ensures that no tugs are on break at t = 1
    for k in K
        JuMP.fix(wqb[1,k], 0; force=true)
    end

    # (21) Ensures that breaks are consecutive
    @constraint(model, c21[t in T[2:end], k in K], sb[t,k] >= wqb[t,k]-wqb[(t-1),k])

    # (22) Ensures that total time on break is equal to break length
    @constraint(model, c22[k in K], sum(sb[t,k] for t in T) == 1)

    # (23) Ensures that total time on break is equal to break length
    @constraint(model, c23[k in K], sum(wqb[t,k] for t in [1:u_lb;]) == b) 

#################################### Warm Start ##########################################    
    
if warmstart == "single"
    include("5. Long Break/long break single.jl")
    ws_objective = single_warmstart(path_to_input, filename, verbose_ws, upper_bound, tugs, K, T, S, P_s, P_q, b, x, y, wvq, wqv, wqq, wqb, sb, u)
end

if warmstart == "grka"
    include("5. Long Break/long break grka.jl")
    ws_objective = grka_warmstart(path_to_input, filename, verbose_ws, sortmethod, upper_bound, u_lb, tugs, K, T, S, Q, P_s, P_q, break_length, x, y, wqq, wqv, wvq, wqb, sb, u)
end

#################################### Optimize ##########################################    
   
    JuMP.optimize!(model)

#################################### Save Results ##########################################    
   
    save_results(model, rule , instance, warmstart, sortmethod, ws_objective, "long break", "results.csv")
    if verbose == true
        save_loading_plan(tugs, upper_bound, x, y, wvq, wqv, wqq, wqb, T, S, K, filename)
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

function save_loading_plan(tugs, upper_bound, x, y, wvq, wqv, wqq, wqb, T, S, K, filename)
    header = ["Time"]
    for k in K
        append!(header, ["Tug $(K[k])"])
    end
    data = Matrix(undef, 0, tugs+1)
    for t in 1:(upper_bound)
        row = ["$(T[t])"]
        for k in K
            if ((JuMP.value(wvq[t,k]))==1)
                push!(row, "$k : V to Q")
            elseif ((JuMP.value(wqv[t,k]))==1)
                push!(row, "$k : Q to V")
            elseif ((JuMP.value(wqq[t,k]))==1)
                push!(row, "$k : IDLE")
            elseif ((JuMP.value(wqb[t,k]))==1)
                push!(row, "$k : BREAK")    
            else
                for i in S
                    if ((JuMP.value(x[i,t,k]))==1)
                        push!(row, "$k : unload $(i)")
                    end 
                
                    if ((JuMP.value(y[i,t,k]))==1)
                        push!(row, "$k : load $(i)")
                    end                 
                end
            end
        end
        data = vcat(data, permutedims(row))
    end
        open(filename*"_loading_plan.txt", "w") do f
        pretty_table(f, data; header=header, alignment=:l)
    end
end