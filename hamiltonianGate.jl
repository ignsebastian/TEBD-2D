function hamiltonianGate(s1,s2)
    #A function to provide the hamiltonian
    hj = op("Sz", s1) * op("Sz", s2) +
              1/2 * op("S+", s1) * op("S-", s2) +
              1/2 * op("S-", s1) * op("S+", s2)
    return hj
end