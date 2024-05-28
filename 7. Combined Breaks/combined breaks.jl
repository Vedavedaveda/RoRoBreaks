#############################################################################
# VARIABLE INDEX
# verbose = if true, solution is saved as .txt
# rule = precedence matrix rule: m, mp, ms, mps
# instance = from 0 to 10
# warmstart = Type of warm start, either "greedy" or "single"
# tugs = k = number of tugs
# timelimit = maximum allowed solver time in seconds
# rins = Relaxation Induced Neighborhood Search (RINS) heuristic frequency (default = -1, off = 0)
# cuts = Global cut control (-1 = default/automatic, 0 = off, 1-3 = least to most aggressive)
#############################################################################


######################################## Packages #############################################

using DataFrames, XLSX, Gurobi, JuMP, Test, CSV
using DataStructures
using ProgressMeter, Dates, LightGraphs
using PrettyTables
import MathOptInterface

######################################## Dual Cycling #############################################

function dualcycling(; verbose = true, verbose_ws = true, rule = "m", instance="1", warmstart="greedy", sortmethod = "top",
                    timelimit=20, rins=-1, cuts = 2, tugs = 4, short_breaklength = 2, rl = 4, ru= 8,
                    long_breaklength = 6, long_break_period_start = 25, long_break_period_end = 35)

    path_to_input = joinpath(@__DIR__, "instances/input/"*rule*"_layout"*instance*".xlsx")

    if warmstart == "grka"
        filename = string("combinedBreaks_", rule, instance, "_", warmstart, "_", sortmethod)
    else
        filename = string("combinedBreaks_", rule, instance, "_", warmstart)
    end

######################################## Parameters #############################################
    
    ### TRAILERS ###

    # Set of trailers to be discharged (S)
    unloadlist = DataFrame(XLSX.readtable(path_to_input, "position_sequence"))
    S = unloadlist[!, :position_onboard] 
	numu = size(unloadlist, 1) # size of unload
    
    # Set of trailers to be loaded (Q)
    loadlist = DataFrame(XLSX.readtable(path_to_input, "position_sequence"))
    Q = loadlist[!, :position_onboard]
	numl= size(loadlist, 1)  # size of load

    ### BREAKS ###

    ### LONG BREAKS ###
    lbl = long_breaklength
    lbps = long_break_period_start
    lbpe = long_break_period_end

    ### SHORT BREAKS ###
    sbl = short_breaklength 
    rl = rl # lower bound for the working time before a break
    ru = ru # upper bound for the working time before a break
    N = [1:floor((numl+numu+1)/rl)-1;]
    
    f = floor(lbps/ru) + 1  # The break index that is sure to occur before lunch break

    ### TIME STEPS ###
    #upper_bound = 2*(numl+numu+1)
    #N = [1:(floor(upper_bound/ru) - 1);]  # Set of breaks
    #upper_bound = upper_bound + long_breaklength + short_breaklength*(length(N))

    upper_bound = numl+numu+1+((length(N)+1)*sbl)+lbl
    T = [1:upper_bound;]

    ### PRECEDENCE ###
    # (i,i') in P
    P_s = DataFrame(XLSX.readtable(path_to_input, "precedence_unloading"))
    pred_u = P_s[!, :pred] 
    suc_u = P_s[!, :suc]
    numpu = size(P_s,1)
    
    # (j,j') in P
    P_q = DataFrame(XLSX.readtable(path_to_input, "precedence_loading")) 
    pred_l = P_q[!, :pred] 
    suc_l = P_q[!, :suc] 
    numpl = size(P_q,1)

    ### TUGS ###
    # Set of tugs (K)
    K = [1:tugs;]
    
######################################### Model #################################

    cd("output")
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit"=>timelimit, "RINS"=>rins, "LogFile"=>filename, "Cuts"=>cuts))

#################################### Decision Variables ##########################################   

    ### TUG ACTIONS ###
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

    ### LONG BREAK ###
    # wqb = 1 if tug k is on the long break at time t
    @variable(model, wqb[t in T, k in K], Bin, start = 0)

    # bs_tk = 1 if long break for tug k starts at t
    @variable(model, lbs[t in T, k in K], Bin, start = 0)

    ### SHORT BREAKS ###
    # wqbs_tkn = 1 if tug k is on break n at t
    @variable(model, wqbs[t in T, k in K, n in N], Bin, start = 0)

    # sbs_tkn = 1 if break n for tug k starts at t
    @variable(model, sbs[t in T, k in K, n in N], Bin, start = 0)

    ### BOUNDS ###
    # (4) Lower bound for makespan
    u_lb = (numu+numl)/tugs+1

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
    @constraint(model, c11[k in K, t in T[2:end]], sum(x[i,t,k] for i in S) + wvq[t,k] == sum(y[j,t-1,k] for j in Q) + wqv[t-1,k])

    # (12) Ensures that tug is doing one thing at a time
    @constraint(model, c12[t in T, k in K], sum(x[i,t,k] for i in S) + wvq[t,k] + sum(y[j,t,k] for j in Q) + wqv[t,k] + wqq[t,k] + wqb[t,k] + sum(wqbs[t,k,n] for n in N) == 1)

    # (13) Enforces maximum of k/2 tugs are discharging a trailer at once
    @constraint(model, c13[t in T], sum(sum(x[i,t,k] for i in S) for k in K) <= tugs/2)

    # (14) Enforces maximum of k/2 tugs are loading a trailer at once
    @constraint(model, c14[t in T], sum(sum(y[j,t,k] for j in Q) for k in K) <= tugs/2)

    # (15) Ensures all tugs are on the quay at the first time step
    @constraint(model, c15[k in K], (sum(y[j,1,k] for j in Q) + wqv[1,k] + wqq[1,k]) == 1) 

    # (16) Ensures no tugs are discharging at the first time step
    for i in S
        for k in K
            JuMP.fix(x[i,1,k], 0; force=true)
        end
    end

    # (17) Ensures no tugs are travelling from ship to quey at first time step
    for k in K
        JuMP.fix(wvq[1,k], 0; force=true)
        JuMP.fix(wqb[1,k], 0; force=true)
    end

    # (18) Ensures that x[i,t] and y[j,t] are binary (not necessary since we defined them as such) (also ensure that sb is binary, and also all the rest (wqq, wqv, etc.)))

    # (19) Ensures that wvq, wqv, wqq, and u are integers (not necessary since we defined them as such)

    # (20) Ensures that no tugs are on break at t = 1
    for k in K
        JuMP.fix(wqb[1,k], 0; force=true)
        JuMP.fix(lbs[1,k], 0; force=true)
        for n in N
            JuMP.fix(wqbs[1,k,n], 0; force=true)
            JuMP.fix(sbs[1,k,n], 0; force=true)
        end
    end

    ### BREAKS ###

    ### LONG BREAK CONSTRAINTS ###
    # (21) Ensures that the long break is consecutive
    @constraint(model, c21[t in T[2:end], k in K], lbs[t,k] >= wqb[t,k]-wqb[(t-1),k] )

    # (22) Ensures that the long break starts exactly once
    @constraint(model, c22[k in K], sum(lbs[t,k] for t in T) == 1)

    # (23) Ensures that the entire long break occurs in the designated period
    @constraint(model, c23[k in K], sum(wqb[t,k] for t in [lbps:lbpe;]) == lbl)

    ### SHORT BREAKS CONSTRAINTS ###
    # (24) Ensures that breaks are consecutive
    @constraint(model, c24[t in T[2:end], k in K, n in N], sbs[t,k,n] >= wqbs[t,k,n]-wqbs[(t-1),k,n])

    # (25) Ensures that exactly one short break exists within each short break index
    @constraint(model, c25[k in K, n in N], sum(sbs[t,k,n] for t in T) == 1)

    # (26) Ensures that total time on break is equal to break length
    @constraint(model, c26[k in K, n in N], sum(wqbs[t,k,n] for t in T) == sbl)
    
    # (27) Ensures that one short break happens within the first interval (THIS CONSTRAINT IS DIFFERENT THAN IN SHORT BREAKS?)
    @constraint(model, c27[k in K], 1 <= sum(t*sbs[t,k,1] for t in 1:ru) <= ru)
    
    # (28) Ensures that the first short break after the long break does not begin too soon (defined by rl) after the long break ends
    @constraint(model, c28[k in K], sum(t*lbs[t,k] for t in T) + lbl + rl <= sum(t*sbs[t,k,f] for t in T))

    # (29) Ensures that the first short break after the long break does not begin too late (defined by ru) after the long break ends
    @constraint(model, c29[k in K], sum(t*sbs[t,k,f] for t in T) <= sum(t*lbs[t,k] for t in T) + lbl + ru)

    # (30) Ensures that all breaks (apart from the first break and the first break after lunch) occur within the appropriate intervals
    @constraint(model, c30[k in K, n in setdiff(N, [1,f])], rl <= sum(t*sbs[t,k,n] for t in T) - sum(t*sbs[t,k,(n-1)] for t in T) <= ru)

#################################### Warm Start ##########################################    
   
    if warmstart == "single"
        include("7. Combined Breaks/combined breaks single.jl")
        ws_objective = single_warmstart(path_to_input, filename, verbose_ws, upper_bound, tugs, K, T, S, P_s, P_q, short_breaklength, rl, N, lbl, lbps, f, x, y, wqq, wqv, wvq, wqb, lbs, wqbs, sbs, u)
    end

    if warmstart == "grka"
        include("7. Combined Breaks/combined breaks grka.jl")
        ws_objective = grka_warmstart(path_to_input, filename, verbose_ws, sortmethod, tugs, S, Q, P_s, P_q, x, y, wqq, wqv, wvq, wqb, lbs, wqbs, sbs, u, upper_bound, u_lb, short_breaklength, K, T, ru, rl, N, f, lbl, lbps, lbpe)
    end    

#################################### Optimize ##########################################    
   
    JuMP.optimize!(model)
    
#################################### Save Results ##########################################    
   
    save_results(model, rule , instance, warmstart, sortmethod, ws_objective, "combined breaks", "results.csv")
    if verbose == true
        save_loading_plan(tugs, upper_bound, x, y, wvq, wqv, wqq, wqb, wqbs, N, T, S, K, filename)
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

function save_loading_plan(tugs, upper_bound, x, y, wvq, wqv, wqq, wqb, wqbs, N, T, S, K, filename)
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
                push!(row, "$k : LONG BREAK")    
            elseif (sum(JuMP.value(wqbs[t,k, n]) for n in N)==1)
                push!(row, "$k : SHORT BREAK")
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