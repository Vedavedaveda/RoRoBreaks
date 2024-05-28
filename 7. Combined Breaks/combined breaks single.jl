function single_warmstart(path_to_input, filename, verbose_ws, upper_bound, tugs, K, T, S, P_s, P_q, short_breaklength, rl, N, lbl, lbps, f, x, y, wqq, wqv, wvq, wqb, lbs, wqbs, sbs, u)

    # Initialise priority queues for loading and discharging sequences
    pq_u, pq_l = get_priority_queues(path_to_input)

    # Create deque of tugs that can be continually recycled
    qk = Deque{Int}()
    for k in K
        pushfirst!(qk, k)
    end

    halftugs = floor(Int, tugs/2)

    # Make half of the tugs idle at t = 1
    for i in 1:halftugs
        k = pop!(qk)                       
        pushfirst!(qk,k)
        set_start_value(wqq[1,k],1)
    end

    # Initialise t
    t = 1
    next_break = rl
    bi = 1              # Break index

    lbo = false            # long break occurance

    ### UNLOAD ###
    while length(pq_u) > 0                  # Do until priority queue of tugs to unload is empty
        t = t + 1                                  # Update time step
        if length(pq_u) <= halftugs                # Checks if: remaining trailers that need to be unlaoded <= half of the tugs
            unload_remaining(P_s, pq_u, qk, halftugs, x, wqq, wqv, t)
        else
            unload_halftugs(P_s, pq_u, qk, halftugs, x, wqv, t)
        end     

        if (t%2 == 0) && (t >= lbps) && (lbo == false)
            lbo = true
            t = make_long_break(qk, halftugs, lbl, wqb, lbs, t)
            next_break = t + rl
        end
        
        if (t%2 == 0) && (t >= next_break) && ((bi < f && lbo == false) || (bi >= f && lbo == true))
            t = make_short_break(qk, halftugs, short_breaklength, wqbs, sbs, bi, t)
            next_break = t + rl
            bi = bi + 1
        end
    end

    ### LOAD ###
    while length(pq_l) > 0                 # Do until priority queue of tugs to load is empty
        t = t + 1                                 # Update time step
        if length(pq_l) <= halftugs               # Checks if: remaining trailers that need to be laoded <= half of the tugs
            load_remaining(P_q, pq_l, qk, halftugs, y, wqq, wvq, t)
        else
            load_halftugs(P_q, pq_l, qk, halftugs, y, wqq, wvq, t)
        end

        if (t%2 == 0) && (t >= lbps) && (lbo == false)
            lbo = true
            t = make_long_break(qk, halftugs, lbl, wqb, lbs, t)
            next_break = t + rl
        end
        if (t%2 == 0) && (t >= next_break) && ((bi < f && lbo == false) || (bi >= f && lbo == true)) && (length(pq_l) > 0)
            t = make_short_break(qk, halftugs, short_breaklength, wqbs, sbs, bi, t)
            next_break = t + rl
            bi = bi + 1
        end
    end

    # Set start value for makespan
    objective = t
    set_start_value(u, objective)

    #=
    for j in 1:halftugs
        k = pop!(qk)                            
        pushfirst!(qk,k)
    end

    # Make tugs idle for the rest of the time steps in T while continuing to include short breaks
    while t <= upper_bound - 1
        t = t + 1
        for j in 1:halftugs
            k = pop!(qk)                            
            pushfirst!(qk,k)
            set_start_value(wqq[t-1,k], 1)
            set_start_value(wqq[t,k], 1)
        end
        if (t >= next_break) && (bi <= length(N))
            for i in 1:halftugs # Give half of the tugs breaks
                k = pop!(qk)                            
                pushfirst!(qk,k)
                set_start_value(sbs[t,k,bi], 1)
                for i in 0:(short_breaklength-1)
                    set_start_value(wqbs[t+i,k,bi], 1)
                end
            end
            t = t + 1
            for i in 1:halftugs # Give the other half of the tugs breaks one time step later
                k = pop!(qk)                            
                pushfirst!(qk,k)
                set_start_value(sbs[t,k,bi], 1)
                for i in 0:(short_breaklength-1)
                    set_start_value(wqbs[t+i,k,bi], 1)
                end
            end
            t = t + short_breaklength - 1
            next_break = next_break + rl
            bi = bi + 1
        end
    end

    for j in 1:halftugs
        k = pop!(qk)                            
        pushfirst!(qk,k)
        set_start_value(wqq[t,k], 1)
    end
    =#
    
    fill_timesteps(upper_bound, tugs, halftugs, qk, short_breaklength, rl, next_break, bi, N, wqq, wqbs, sbs, t)

    if verbose_ws == true
        print_single_warmstart(tugs, upper_bound, x, y, wvq, wqv, wqq, wqb, wqbs, N, T, S, K, filename)
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

function set_idle(idle_tugs, qk, wqq, t)
    for i in 1:idle_tugs            # For leftover tugs that did not unload at t
        k = pop!(qk)                             # Get next tug k
        pushfirst!(qk,k)                         # Send to back of deque
        set_start_value(wqq[(t-1),k], 1)         # Make tug idle at t-1 and t
        set_start_value(wqq[t,k], 1)
    end
end

function make_short_break(qk, halftugs, short_breaklength, wqbs, sbs, bi, t)
    for i in 1:halftugs # Give half of the tugs breaks
        k = pop!(qk)                            
        pushfirst!(qk,k)
        set_start_value(sbs[t,k,bi], 1)
        for i in 0:(short_breaklength-1)
            set_start_value(wqbs[t+i,k,bi], 1)
        end
    end
    t = t + 1
    for i in 1:halftugs # Give the other half of the tugs breaks one time step later
        k = pop!(qk)                            
        pushfirst!(qk,k)
        set_start_value(sbs[t,k,bi], 1)
        for i in 0:(short_breaklength-1)
            set_start_value(wqbs[t+i,k,bi], 1)
        end
    end
    return t + short_breaklength - 1
end

function make_long_break(qk, halftugs, lbl, wqb, lbs, t)
    for i in 1:halftugs # Give half of the tugs breaks
        k = pop!(qk)                            
        pushfirst!(qk,k)
        println("tug $k going on long break at $t")
        set_start_value(lbs[t,k], 1)
        for i in 0:(lbl-1)
            set_start_value(wqb[t+i,k], 1)
        end
    end
    t = t + 1
    for i in 1:halftugs # Give the other half of the tugs breaks one time step later
        k = pop!(qk)                            
        pushfirst!(qk,k)
        println("tug $k going on long break at $t")
        set_start_value(lbs[t,k], 1)
        for i in 0:(lbl-1)
            set_start_value(wqb[t+i,k], 1)
        end
    end
    return t + lbl - 1
end

function unload_halftugs(P_s, pq_u, qk, halftugs, x, wqv, t)
    for i in 1:halftugs                   # For next group of half tugs
        workingtugs = 0
        k = pop!(qk)                             # get next available tug
        pushfirst!(qk,k)                         # add tug to back of queue
        unload_now = dequeue_pair!(pq_u)[1]          # get next trailer to unload
        set_start_value(wqv[(t-1),k], 1)         # make tug travel to vessel at t-1
        set_start_value(x[unload_now,t,k], 1)        # make tug unload trailer at t
        workingtugs = workingtugs + 1

        unload_next = peek(pq_u)[1]
        prec_violation = precedence_check(P_s, unload_now, unload_next)
        if prec_violation == 1
            break
        end
    end 
end

function unload_remaining(P_s, pq_u, qk, halftugs, x, wqq, wqv, t)
    remaining = length(pq_u)               # Store number of remaining trailers to be unloaded
    for i in 1:remaining                   # For remaining trailers
        workingtugs = 0
        k = pop!(qk)                             # Get next available tug k
        pushfirst!(qk,k)                         # Send tug to back of deque
        unload_now = dequeue_pair!(pq_u)[1]          # Get next trailer to be unloaded
        set_start_value(wqv[(t-1),k], 1)         # Tug travels from quay to vessel at t-1
        set_start_value(x[unload_now,t,k], 1)        # Tug k unloads trailer at t
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
    set_idle(halftugs-remaining, qk, wqq, t)
end

function load_halftugs(P_q, pq_l, qk, halftugs, y, wqq, wvq, t)
    for i in 1:halftugs                   # For next group of half tugs
        workingtugs = 0
        k = pop!(qk)                              # get next available tug
        pushfirst!(qk,k)                          # add tug to back of queue
        load_now = dequeue_pair!(pq_l)[1]             # get next trailer to load
        set_start_value(y[load_now,(t-1),k], 1)       # make tug load trailer at t-1
        set_start_value(wvq[t,k], 1)              # make tug travel to quay at t
        workingtugs = workingtugs + 1

        load_next = peek(pq_l)[1]
        prec_violation = precedence_check(P_q, load_now, load_next)

        if prec_violation == 1
            set_idle(halftugs-workingtugs, qk, wqq, t)
            break
        end
    end
end

function load_remaining(P_q, pq_l, qk, halftugs, y, wqq, wvq, t)
    remaining = length(pq_l)              # Store number of remaining trailers to be unloaded       
    for i in 1:remaining                  # For remaining trailers    
        workingtugs = 0   
        k = pop!(qk)                             # Get next available tug k
        pushfirst!(qk,k)                         # Send tug to back of deque
        load_now = dequeue_pair!(pq_l)[1]            # Get next trailer to be loaded
        set_start_value(y[load_now,(t-1),k], 1)      # Tug k loads trailer at t-1
        set_start_value(wvq[t,k], 1)             # Tug travels from vessel to quay at t
        workingtugs = workingtugs + 1
        if length(pq_l) == 0 
            break
        end
        load_next = peek(pq_l)[1]
        prec_violation = precedence_check(P_q, load_now, load_next)

        if prec_violation == 1
            set_idle(halftugs-workingtugs, qk, wqq, t)
        end
    end
    set_idle(halftugs-remaining, qk, wqq, t)
    for i in 1:halftugs                   # For remaining group of half tugs
        k = pop!(qk)                              # Get next tug k
        pushfirst!(qk,k)                          # Send to back of deque
        set_start_value(wqq[t,k], 1)              # Make k idle at t
    end
end

function fill_timesteps(upper_bound, tugs, halftugs, qk, breaklength, rl, next_break, bi, N, wqq, wqb, sb, t)
    while t <= upper_bound - 1
        t = t + 1
        if (t%2 == 0) && (t >= next_break) && (bi <= length(N))
            for i in 1:halftugs # Give half of the tugs breaks
                k = pop!(qk)                            
                pushfirst!(qk,k)
                set_start_value(sb[t,k,bi], 1)
                set_start_value(wqq[t+breaklength,k], 1)
                for j in 0:(breaklength-1)
                    set_start_value(wqb[t+j,k,bi], 1)
                end
            end
            t = t + 1
            for i in 1:halftugs # Give the other half of the tugs breaks one time step later
                k = pop!(qk)                            
                pushfirst!(qk,k)
                set_start_value(wqq[t-1,k], 1)
                set_start_value(sb[t,k,bi], 1)
                for j in 0:(breaklength-1)
                    set_start_value(wqb[t+j,k,bi], 1)
                end
            end
            t = t + breaklength - 1
            next_break = t + rl
            bi = bi + 1
        else
            for j in 1:tugs
                k = pop!(qk)                            
                pushfirst!(qk,k)
                set_start_value(wqq[t,k], 1)
            end
        end
    end

    for j in 1:halftugs
        k = pop!(qk)                            
        pushfirst!(qk,k)
        set_start_value(wqq[t,k], 1)
    end
end

function print_single_warmstart(tugs, upper_bound, x, y, wvq, wqv, wqq, wqb, wqbs, N, T, S, K, filename)
    header = ["Time"]
    for k in K
        append!(header, ["Tug $(K[k])"])
    end
    data = Matrix(undef, 0, tugs+1)
    for t in 1:(upper_bound)
        row = ["$(T[t])"]
        for k in K
            if ((start_value(wvq[t,k]))==1)
                push!(row, "$k : V to Q")
            end
            if ((start_value(wqv[t,k]))==1)
                push!(row, "$k : Q to V")
            end
            if ((start_value(wqq[t,k]))==1)
                push!(row, "$k : IDLE")
            end
            if ((start_value(wqb[t,k]))==1)
                push!(row, "$k : LONG BREAK")    
            end
            if (sum(start_value(wqbs[t,k, n]) for n in N)==1)
                push!(row, "$k : SHORT BREAK")
            else
                for i in S
                    if ((start_value(x[i,t,k]))==1)
                        push!(row, "$k : unload $(i)")
                    end 
                
                    if ((start_value(y[i,t,k]))==1)
                        push!(row, "$k : load $(i)")
                    end                 
                end
            end
        end
        data = vcat(data, permutedims(row))
    end
        open(filename*"_ws_plan.txt", "w") do f
        pretty_table(f, data; header=header, alignment=:l)
    end
end