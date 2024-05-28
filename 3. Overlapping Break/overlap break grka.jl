function grka_warmstart(path_to_input, filename, verbose_ws, sortmethod, upper_bound, tugs, T, S, Q, P_s, P_q, bs, bl, x, y, wqq, wqv, wvq, u)
    
    ### ASSIGN/GET VARIABLES ###
    t = 1
    job_dict = get_job_dict(S, Q) # make a dictionary, that includes each operation and its unique ID, {operation: ID}
    DAG = get_complete_precedence(job_dict, S, P_s, P_q)  # build DAG - complete precedence (UU, LL, UL)
    job_dict_reversed = Dict(value => key for (key, value) in job_dict) # job dict reversed, {ID:operation}

    ### SORT THE PRECEDENCE ###
    # we want to first unload/load the trailers with the most descendants, that is the tailers that 'block' the most operations
    if sortmethod == "kids"
        job_precedence = sort_kids(DAG)
    elseif sortmethod == "top"
        job_precedence = reverse(topological_sort_by_dfs(DAG))
    end

    ### GET VARIABLES FOR STORING TUGS' ACTIONS ###
    tugs_array = [String[],String[],String[],String[]] # for storing each tug's action 
    for i in (1): (upper_bound) # at the beggining assign each action to being idle ('WQQ')
        for tug in tugs_array
            push!(tug,"WQQ")
        end
    end
    
    ### PERFORM WARMSTART ### 
    # this function includes unloading/loading of the vessel and assigning correct tug actions for each tug
    t = get_warmstart_variables(job_precedence, job_dict_reversed, DAG, tugs_array, t) 
    # correction of the last time step - in case it ends with loading, and the tug never returns to the quay 
    add_timestep = false
    for tug in tugs_array
        if split(tug[t], " ")[1] == "load"
            add_timestep = true
            tug[t+1] = "WVQ"
        end
    end
    if add_timestep t = t + 1 end

    ### ADD BREAKS ###
    break_tugs_array = [String[],String[],String[],String[]]
    add_break(tugs_array, break_tugs_array, bs, bl, tugs)

    ### ASSIGN MODEL'S VARIABLES WITH WARMSTART OUTPUT ###
    #assign_variables(tugs_array, x, y, wqq, wqv, wvq, upper_bound)
    assign_variables(break_tugs_array, x, y, wqq, wqv, wvq, upper_bound)
    set_start_value(u, t+bl)

    ### PRINT THE RESULTS ###
    if verbose_ws == true
        print_warmstart(tugs, upper_bound, x, y, wvq, wqv, wqq, T, S, filename)
        println("ws objective = ", t+bl)
    end

    return t+bl
    
end

function get_job_dict(S, Q) # assign each operation an ID in a dict: {position:ID} (ordered IDs for both S and Q)
    job_dict = Dict{String, Int}()
    for (operation_id_S, position_S) in enumerate(S)
        push!(job_dict, position_S*"U" => operation_id_S)
    end

    len_job_dict = length(job_dict)

    for (operation_id_Q, position_Q) in enumerate(Q)
        push!(job_dict, position_Q*"L" => operation_id_Q+len_job_dict)
    end
    return job_dict
end

function get_complete_precedence(job_dict, S, P_s, P_q) 
    # make a precedence graph DAG (complete precedence for all unloading and loading rules)
    # all the vertices in the graph are operation IDs 
    DAG = DiGraph(length(job_dict))
    @showprogress 1 "Building graph" for (child, parent) in eachrow(P_s) # P_s is the precedence matrix P (i, i')
        add_edge!(DAG, job_dict[parent*"U"], job_dict[child*"U"])
    end
    for (child, parent) in eachrow(P_q)  # P_q is the precedence matrix P (j, j')
        add_edge!(DAG, job_dict[parent*"L"], job_dict[child*"L"])
    end
    for job in S
        add_edge!(DAG, job_dict[job*"L"], job_dict[job*"U"])
    end
    return DAG
end

function sort_kids(DAG) # sort by number of descendants/kids
    job_precedence_kids = Dict{Int64, Int64}()
    for vertex in vertices(DAG)
        shortest_paths = dijkstra_shortest_paths(DAG, vertex).dists # list of shortest paths between 'vertex' and other vertices in DAG
        finite_dist = filter(x -> (x != typemax(x) && x != 0), shortest_paths)
        len_finite_dist = length(finite_dist)
        push!(job_precedence_kids, vertex => len_finite_dist)
    end
    job_precedence_kids = sort(collect(job_precedence_kids), by=x->x[2])
    
    job_precedence = Int64[]
    for (key, value) in job_precedence_kids
        push!(job_precedence, key)
    end
    return job_precedence
end

function assign_variables(tugs_array, x, y, wqq, wqv, wvq, upper_bound)
    for t in 1:upper_bound
        wqq_tugs = 0
        wqv_tugs = 0
        wvq_tugs = 0
        for tug in tugs_array
            if tug[t] == "WQQ"
                wqq_tugs = wqq_tugs + 1
            elseif tug[t] == "WQV"
                wqv_tugs = wqv_tugs + 1
            elseif tug[t] == "WVQ"
                wvq_tugs = wvq_tugs + 1
            elseif split(tug[t], " ")[1] == "load"
                slot_name = split(tug[t], " ")[2]
                set_start_value(y[slot_name, t], 1)
            elseif split(tug[t], " ")[1] == "unload"
                slot_name = split(tug[t], " ")[2]
                set_start_value(x[slot_name, t], 1)
            else
                println("Something wrong.")
            end
        end
        set_start_value(wqq[t], wqq_tugs)
        set_start_value(wqv[t], wqv_tugs)
        set_start_value(wvq[t], wvq_tugs)
    end

end

function get_warmstart_variables(job_precedence, job_dict_reversed, DAG, tugs_array, t) # main loop for the wamstart
    ### ASSIGN VARIABLES ###
    k = 4
    last_load = 0
    last_unload = 0
    load_available = [] # tug saying: hi, I can load something :)
    unload_available = [1, 2, 3, 4] # start with all tugs being available for unloading

    while length(job_precedence) > 0 # loop through the IDs of operations, and deal with them one by one
        number_of_unload = Int64[] # keep track of the IDs of unloaded trailers in each timestamp
        number_of_load = Int64[] # keep track of the IDs of loaded trailers in each timestamp
        t = t + 1
        # loop through the job precedence list, and decide which trailers should be loaded/unloaded in this time step
        check_load_unload(job_precedence, job_dict_reversed, DAG, load_available, unload_available, number_of_unload, number_of_load, tugs_array, k, last_load, last_unload, t)
        @label label2
        w = 0

        assign_missing_load_unload(tugs_array, t) # assign missing values for returning from / coming into the vessel
        assign_tugs_availability(tugs_array, load_available, unload_available, t) # make sure correct tugs are in the load/unload available lists

        # keep track of how many tugs where loading/unloading in this time step before moving to the next one one 
        last_load = length(number_of_load)
        last_unload = length(number_of_unload)
        # remove all loaded/unloaded trailers in this time step from the precedence list
        filter!(e->!(e in number_of_unload), job_precedence)
        filter!(e->!(e in number_of_load), job_precedence)
    end
    return t
end

function check_load_unload(job_precedence, job_dict_reversed, DAG, load_available, unload_available,number_of_unload, number_of_load, tugs_array_in, k, last_load, last_unload, t)
    for (precedence_index, precedence_value) in enumerate(job_precedence)
        # get information about the operation - if it's load or unload
        operation_full = job_dict_reversed[precedence_value]
        operation = operation_full[1:end-1]
        string(operation_full[end]) == "U" ? mode = "U" : mode = "L" # if the operation is unloading then mode = U, mode = L otherwise
        
        # in this loop check whether our current operation is blocked by any other operation
        break_precedence = false
        for job_b in job_precedence
            if job_b in vcat(number_of_unload, number_of_load)
                nothing
            elseif job_b == precedence_value
                nothing
            else
                has_path(DAG, precedence_value, job_b) ? break_precedence=true : nothing # if it is not possible to load/unload take another element
            end
        end
        if break_precedence continue end

        @goto label3 # if exhausion, meaning not blocked by any, then go to label 3 and load/unload
        @label label3

        ########### UNLOADING ########
        if mode == "U" 
            # check if unload-unload constraint is violated
            if length(number_of_unload) > 0
                for pq_done in number_of_unload
                    has_path(DAG, precedence_value, pq_done) ? (@goto label1) : nothing
                end
            end
            
            # if all the tugs have already been assigned an action, move onto next time step
            length(vcat(number_of_unload, number_of_load)) == k || length(number_of_unload)>= k/2 || length(number_of_unload) >= last_load + (k-last_load-last_unload) ? (return) : nothing
            # get an available tug and assign the tug to unloading the current trailer
            isempty(unload_available) ? (continue) : nothing
            unloading_tug = pop!(unload_available)
            # if tug is both in unload and load lists, remove it from here: load available
            if unloading_tug in load_available filter!(e->(e != unloading_tug), load_available) end
            tugs_array_in[unloading_tug][t] = "unload $operation"
            push!(number_of_unload, precedence_value)

        ############ LOADING #########
        else 
            # check if loading-loading or unloading-loading constraint are violated
            if length(number_of_load) > 0
                for pq_done in number_of_load
                    job_dict_reversed[pq_done][1:end-1] == job_dict_reversed[precedence_value][1:end-1] ? continue : nothing
                    has_path(DAG, precedence_value, pq_done) ? (@goto label1) : nothing
                end
            end
            # similarly to unload - if none of the 'tug constraints' are violated, proceed with loading of the current trailer
            length(vcat(number_of_unload,number_of_load)) == k || length(number_of_load)>= k/2 || length(number_of_load) >= last_unload + (k-last_load-last_unload) ? (return) : nothing
            isempty(load_available) ? (continue) : nothing
            loading_tug = pop!(load_available)
            if loading_tug in unload_available
                filter!(e->(e != loading_tug), unload_available)
            end
            tugs_array_in[loading_tug][t] = "load $operation"
            push!(number_of_load, precedence_value)
        end
        @label label1
        unloading_tug = nothing
        loading_tug = nothing
    end
end

function assign_missing_load_unload(tugs_array_in, t) # make sure tugs are returning from / coming to the vessel if necessary
    for tug in tugs_array_in
        if split(tug[t], " ")[1] == "unload" && split(tug[t-1], " ")[1] != "load"
            tug[t-1] = "WQV"
        elseif split(tug[t-1], " ")[1] == "load" && split(tug[t], " ")[1] != "unload"
            tug[t] = "WVQ"
        end
    end
end

function assign_tugs_availability(tugs_array_in, load_available, unload_available, t) # put tugs into the load/unload available lists
    for tug_index in 1:(length(tugs_array_in))
        if split(tugs_array_in[tug_index][t], " ")[1] == "unload"
            tug_index in load_available ? nothing : push!(load_available, tug_index)
        elseif split(tugs_array_in[tug_index][t], " ")[1] == "load"
            tug_index in unload_available ? nothing : push!(unload_available, tug_index)
        elseif tugs_array_in[tug_index][t] == "WVQ"
            tug_index in load_available ? nothing : push!(load_available, tug_index)
        elseif tugs_array_in[tug_index][t] =="WQV"
            nothing
        else
            tugs_array_in[tug_index][t] = "WQQ"
            tug_index in load_available ? nothing : push!(load_available, tug_index)
            tug_index in unload_available ? nothing : push!(unload_available, tug_index)
        end
    end
end

function add_break(tugs_array, break_tugs_array, bs, bl, tugs)
    for (tug_index, tug) in enumerate(tugs_array)
        if tug_index <= tugs/2
            if isodd(bs)
                bs1 = bs - 1
                bl1 = bl
            else
                bs1 = bs
                bl1 = bl - 1
            end
            for (index, action) in enumerate(tug)
                if index == bs1
                    if action == "WQQ" || action == "WQV" || split(action, " ")[1] == "load"
                        push!(break_tugs_array[tug_index], "WQQ")
                        for i in 1:bl1
                            push!(break_tugs_array[tug_index], "WQQ")
                        end
                        push!(break_tugs_array[tug_index], action)
                    else
                        push!(break_tugs_array[tug_index], action)
                        push!(break_tugs_array[tug_index], "WQQ")
                        for i in 1:bl1
                            push!(break_tugs_array[tug_index], "WQQ")
                        end
                    end
                else
                    push!(break_tugs_array[tug_index], action)
                end
            end
        else
            if isodd(bs)
                bs2 = bs
                bl1 = bl
            else
                bs2 = bs + 1
                bl1 = bl - 1
            end
            for (index, action) in enumerate(tug)
                if index == bs2
                    if action == "WQQ" || action == "WQV" || split(action, " ")[1] == "load"
                        push!(break_tugs_array[tug_index], "WQQ")
                        for i in 1:bl1
                            push!(break_tugs_array[tug_index], "WQQ")
                        end
                        push!(break_tugs_array[tug_index], action)
                    else
                        push!(break_tugs_array[tug_index], action)
                        push!(break_tugs_array[tug_index], "WQQ")
                        for i in 1:bl1
                            push!(break_tugs_array[tug_index], "WQQ")
                        end
                    end
                else
                    push!(break_tugs_array[tug_index], action)
                end
            end        
        end
    end
end

function print_warmstart(tugs, upper_bound, x, y, wvq, wqv, wqq, T, S, filename)
    header = ["Time", "V to Q", "Q to V", "IDLE", "", "", "", ""]
    data = Matrix(undef, 0, 2*tugs)
    for t in 1:(upper_bound)
        row = ["$(T[t])"]
        push!(row, "V to Q = $(start_value(wvq[t]))")
        push!(row, "Q to V = $(start_value(wqv[t]))")
        push!(row, "IDLE = $(start_value(wqq[t]))")
        ul_actions = tugs
        for i in S
            if ((start_value(x[i,t]))==1)
                push!(row, "unload $(i)")
                ul_actions = ul_actions - 1
            end 
            if ((start_value(y[i,t]))==1)
                push!(row, "load $(i)")
                ul_actions = ul_actions - 1
            end                 
        end
        for j in 1:ul_actions
            push!(row, "")
        end
        data = vcat(data, permutedims(row))
    end
        open(filename*"_ws_plan.txt", "w") do f
        pretty_table(f, data; header=header, alignment=:l)
    end
end