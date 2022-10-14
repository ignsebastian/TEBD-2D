function trotterGate(s,lattice,Ny,tau)
    #Function to build all trotter gates to perform TEBD
    #s = sites
    #lattice = all nearest neighbor pairs
    #Ny = number of sites on y direction
    #tau = time step
    #-------------Defining all variables-----------------------
    
    diff = [] #find all differences between next nearest neigbor sites (e.g, if there is a pair of site 1 and site 3 as nearest neighbor sites, the difference will be 2)
    pairs = [[],[],[],[],[]] #all possible group of pairs in square and triangular lattices
    
    H = [] #to contain gates for each Hamiltonian terms (H_odd, H_even, etc.). All the terms inside the group must be commute with each other
    
    gate = ITensor[] #To combine all gates
    
    #--------------Find all differences-----------------
    for b in lattice
        p = [b.s1,b.s2]
        d = abs(b.s1 - b.s2)
        if !(d in diff)
            push!(diff,d)
        end

        if d == 1
            push!(pairs[1],p)
        elseif d == Ny-1
            push!(pairs[2],p)
        elseif d == Ny
            push!(pairs[3],p)
        elseif d == Ny+1
            push!(pairs[4],p)
        elseif d == 2*Ny-1
            push!(pairs[5],p)
        end
    end
    
    pairs = [arr for arr in pairs if !all(isempty.(arr))] #to delete empty pairs
    
    #-------------Build trotter gate--------------------
    #Building gate for each possible group of Hamiltonian
    
    if 1 in diff
        h1 = ITensor[]
        h2 = ITensor[]
        for pair in pairs[1]
            s1 = s[pair[1]]
            s2 = s[pair[2]]
            hj = hamiltonianGate(s1,s2)
            Gj = exp(tau/2 * hj)
            if mod(pair[1],2) == 1
                push!(h1, Gj)
            elseif Ny != 2 && mod(pair[1],2) == 0
                push!(h2, Gj)
            end
        end
        
        push!(H,h1)
        push!(H,h2)
    end
    
    if Ny > 1
        if Ny-1 in diff
            if Ny-1 != Ny
                h1 = ITensor[]
                h2 = ITensor[]
                if (Ny-1) in diff
                    for pair in pairs[2]
                        U = swapSequence(s,pair[1],pair[2])
                        s1 = s[pair[1]+(Ny-1-1)]
                        s2 = s[pair[2]]
                        hj = hamiltonianGate(s1,s2)
                        Gj = exp(tau/2 * hj)
                        append!(h,U)
                        push!(h,Gj)
                        append!(h1, reverse(U))
                    end
                end
                push!(H,h1)
            end
        end

        if Ny in diff
            h1 = ITensor[]
            h2 = ITensor[]
            used_pair = []
            for pair in pairs[3]
                U = swapSequence(s,pair[1],pair[2])
                s1 = s[pair[1]+(Ny-1)]
                s2 = s[pair[2]]
                hj = hamiltonianGate(s1,s2)
                Gj = exp(tau/2 * hj)

                if !(pair[1] in used_pair) && !(pair[2] in used_pair)
                    append!(h1,U)
                    push!(h1,Gj)
                    append!(h1,reverse(U))

                    push!(used_pair,pair[1])
                    push!(used_pair,pair[2])
                else
                    append!(h2,U)
                    push!(h2,Gj)
                    append!(h2,reverse(U))
                end
            end
            push!(H,h1)
            push!(H,h2)
        end

        if Ny+1 in diff
            h1 = ITensor[]
            h2 = ITensor[]
            used_pair = []
            for pair in pairs[4]
                U = swapSequence(s,pair[1],pair[2])
                s1 = s[pair[1]+(Ny)]
                s2 = s[pair[2]]
                hj = hamiltonianGate(s1,s2)
                Gj = exp(tau/2 * hj)

                if !(pair[1] in used_pair) && !(pair[2] in used_pair)
                    append!(h1,U)
                    push!(h1,Gj)
                    append!(h1,reverse(U))

                    push!(used_pair,pair[1])
                    push!(used_pair,pair[2])
                else
                    append!(h2,U)
                    push!(h2,Gj)
                    append!(h2,reverse(U))
                end
            end
            push!(H,h1)
            push!(H,h2)
        end
        
        if 2*Ny-1 in diff
            h1 = ITensor[]
            h2 = ITensor[]
            used_pair = []
            for pair in pairs[5]
                U = swapSequence(s,pair[1],pair[2])
                s1 = s[pair[1]+(2*Ny-2)]
                s2 = s[pair[2]]
                hj = hamiltonianGate(s1,s2)
                Gj = exp(tau/2 * hj)

                if !(pair[1] in used_pair) && !(pair[2] in used_pair)
                    append!(h1,U)
                    push!(h1,Gj)
                    append!(h1,reverse(U))

                    push!(used_pair,pair[1])
                    push!(used_pair,pair[2])
                else
                    append!(h2,U)
                    push!(h2,Gj)
                    append!(h2,reverse(U))
                end
            end
            push!(H,h1)
            push!(H,h2)
        end
    end
    
    
    #Combine all group
    
    for H_gate in H
        if !isempty(H_gate)
            append!(gate,H_gate)
        end
    end
    append!(gate,reverse(gate))
    
    return gate
end
