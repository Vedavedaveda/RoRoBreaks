function single_warmstart(path_to_input, filename, verbose_ws, upper_bound, tugs, T, S, P_s, P_q, bs, bl, x, y, wqq, wqv, wvq, u)
    # Initialise priority queues for loading and discharging sequences
    pq_u, pq_l = get_priority_queues(path_to_input)

    # Make half of the tugs idle at t = 1
    halftugs = floor(Int, tugs/2)
    set_start_value(wqq[1], halftugs)

    # Initialise t
    t = 1

    ### UNLOAD ###
    while length(pq_u) > 0                 # Do until priority queue of tugs to unload is empty
        t = t + 1
        if length(pq_u) <= halftugs                # Checks if: remaining trailers that need to be unlaoded <= half of the tugs
            unload_remaining(P_s, pq_u, halftugs, x, wqq, wqv, t)
        else
            unload_halftugs(P_s, pq_u, halftugs, x, wqv, t)
        end 
        if t == bs
            t = make_break(tugs, halftugs, bs, bl, wqq, t)
        end        
    end

    ### LOAD ###
    while length(pq_l) > 0                 # Do until priority queue of tugs to load is empty
        t = t + 1                                 # Update time step
        if length(pq_l) <= halftugs               # Checks if: remaining trailers that need to be laoded <= half of the tugs
            load_remaining(P_q, pq_l, halftugs, y, wqq, wvq, t)
        else
            load_halftugs(P_q, pq_l, halftugs, y, wqq, wvq, t)
        end
        if t == bs
            t = make_break(tugs, halftugs, bs, bl, wqq, t)
        end
    end
    
    # Make tugs idle for the rest of the time steps in T
    for i in (t+1):upper_bound
        set_start_value(wqq[i], tugs) 
    end

    # Set start value for makespan
    objective = t
    set_start_value(u, objective)

    if verbose_ws == true
        print_single_warmstart(tugs, upper_bound, x, y, wvq, wqv, wqq, T, S, filename)
        println("ws objective = ", start_value(u))
    end
    return objective
    
end

function get_priority_queues(path_to_input)
    # Create priority queues of the unload and load sequences, ordered by precedence constraints
    pos_seq = DataFrame(XLSX.readtable(path_to_input, "position_sequence"))
    ws_u = select(pos_seq, [:position_onboard, :dis_seq])
    ws_l = select(pos_seq, [:position_onboard, :load_seq])

    pq_u = PriorityQueue()
    for (position_onboard, dis_seq) in eachrow(ws_u)
        pq_u[position_onboard] = dis_seq
    end

    pq_l = PriorityQueue()
    for (position_onboard, load_seq) in eachrow(ws_l)
        pq_l[position_onboard] = load_seq
    end

    return pq_u, pq_l
end

function precedence_check(P, now, next)
    prec_violation = 0
    for (suc, pred) in eachrow(P)
        if (now, next) == (suc, pred)
            prec_violation = 1
            break
        end
    end
    return prec_violation
end

function set_idle(idle_tugs, wqq, t)
    for i in 1:idle_tugs  
        set_start_value(wqq[(t-1)], start_value(wqq[(t-1)]) + 1)    
        set_start_value(wqq[t], start_value(wqq[t]) + 1)
    end
end

function make_break(tugs, halftugs, bs, bl, wqq, t)
    set_start_value(wqq[t], start_value(wqq[t]) + halftugs)
    for b in (bs+1):(bs+bl-1)
        set_start_value(wqq[b], tugs)
    end
    set_start_value(wqq[bs+bl], halftugs)
    return bs + bl
end

function unload_halftugs(P_s, pq_u, halftugs, x, wqv, t)
    for i in 1:halftugs                   # For next group of half tugs
        workingtugs = 0
        unload_now = dequeue_pair!(pq_u)[1]          # get next trailer to unload
        set_start_value(wqv[(t-1)], start_value(wqv[(t-1)]) + 1)         # make tug travel to vessel at t-1
        set_start_value(x[unload_now,t], 1)        # make tug unload trailer at t
        workingtugs = workingtugs + 1

        unload_next = peek(pq_u)[1]
        prec_violation = precedence_check(P_s, unload_now, unload_next)
        if prec_violation == 1
            break
        end
    end
end

function unload_remaining(P_s, pq_u, halftugs, x, wqq, wqv, t)
    remaining = length(pq_u)               # Store number of remaining trailers to be unloaded
    for i in 1:remaining                   # For remaining trailers
        workingtugs = 0                             # NOT SURE IF THIS NEEDS TO BE OUTSIDE LOOP
        unload_now = dequeue_pair!(pq_u)[1]          # Get next trailer to be unloaded
        set_start_value(wqv[(t-1)], start_value(wqv[(t-1)]) + 1)         # Tug travels from quay to vessel at t-1
        set_start_value(x[unload_now,t], 1)        # Tug k unloads trailer at t
        workingtugs = workingtugs + 1
        if length(pq_u) == 0 
            break
        end
        unload_next = peek(pq_u)[1]
        prec_violation = precedence_check(P_s, unload_now, unload_next)
        if prec_violation == 1
            break
        end
    end
    set_idle(halftugs-remaining, wqq, t)
end

function load_halftugs(P_q, pq_l, halftugs, y, wqq, wvq, t)
    for i in 1:halftugs                   # For next group of half tugs
        workingtugs = 0
        load_now = dequeue_pair!(pq_l)[1]             # get next trailer to load
        set_start_value(y[load_now,(t-1)], 1)       # make tug load trailer at t-1
        set_start_value(wvq[t], start_value(wvq[t]) + 1)              # make tug travel to quay at t
        workingtugs = workingtugs + 1

        load_next = peek(pq_l)[1]
        prec_violation = precedence_check(P_q, load_now, load_next)
        
        if prec_violation == 1
            set_idle(halftugs - workingtugs, wqq, t)
            break
        end
    end
end

function load_remaining(P_q, pq_l, halftugs, y, wqq, wvq, t)
    remaining = length(pq_l)              # Store number of remaining trailers to be unloaded       
    for i in 1:remaining                  # For remaining trailers    
        workingtugs = 0   
        load_now = dequeue_pair!(pq_l)[1]            # Get next trailer to be loaded
        set_start_value(y[load_now,(t-1)], 1)      # Tug k loads trailer at t-1
        set_start_value(wvq[t], start_value(wvq[t]) + 1)             # Tug travels from vessel to quay at t
        workingtugs = workingtugs + 1
        if length(pq_l) == 0 
            break
        end
        load_next = peek(pq_l)[1]
        prec_violation = precedence_check(P_q, load_now, load_next)

        if prec_violation == 1
            set_idle(halftugs-workingtugs, wqq, t)
        end
    end
    set_idle(halftugs-remaining, wqq, t)
    for i in 1:halftugs                   # For remaining group of half tugs
        set_start_value(wqq[t], start_value(wqq[t]) + 1)              # Make k idle at t
    end
end

function print_single_warmstart(tugs, upper_bound, x, y, wvq, wqv, wqq, T, S, filename)
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