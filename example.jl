using ITensors
using Random

include("swapFunction.jl")
include("hamiltonianGate.jl")
include("trotterGate.jl")

#The example is referring to https://itensor.github.io/ITensors.jl/dev/tutorials/MPSTimeEvolution.html

Random.seed!(2123)

let
    #-------- Parameters --------#
    Nx = 100
    Ny = 1
    N = Nx*Ny
    cutoff = 1E-8
    tau = 0.1
    ttotal = 5.0

    #-------- Define site
    s = siteinds("S=1/2", N; conserve_qns=true)
    lattice = square_lattice(Nx,Ny,yperiodic=true)

    #------- Make gates (1,2),(2,3),(3,4),...
    gates = trotterGate(s,lattice,Ny,tau)

    #------- Initialize random psi
    psi = randomMPS(s,linkdims=2)

    c = div(N, 2) # center site

    #------- Applying trotter gate and calculating 
    #------- expectation value of Sz on center site
    #
    for t in 0.0:tau:ttotal
        Sz = expect(psi, "Sz"; sites=c)
        println("$t $Sz")
        t≈ttotal && break

        psi = apply(gates, psi; cutoff)
        normalize!(psi)
    end

end
