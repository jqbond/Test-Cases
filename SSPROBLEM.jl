using DifferentialEquations
using Plots

function ODESYSTEM(du, u, p, t)
    #Relabel state variables; thetas are fractional surface coverages of A and B
    thetaA = u[1];
    thetaB = u[2];
    
    #Calculate vacant sites; fractional coverages must sum to 1.
    thetav = 1 - thetaA - thetaB;

    #Fixed bulk concentrations of CA and CB; driving force for adsorption
    CA     = 1;
    CB     = 1;

    k1f, k1r, k2f, k2r, k3f, k3r, ssabstol, ssreltol = p

    #Calculate rates of all 3 reactions; functions of concentration and/or coverage
    r1 = k1f*CA*thetav - k1r*thetaA;
    r2 = k2f*thetaA    - k2r*thetaB;
    r3 = k3f*thetaB    - k3r*CB*thetav;

    RtA =  r1 - r2;
    RtB =  r2 - r3;

    du[1] = RtA;
    du[2] = RtB;
   
    #I added this bit to display when the criteria were met; it's pretty early on in the integration
    #maybe t ~ 1e-2 or so.

    #test1 = (abs(du[1]) > ssabstol) && (abs(du[1]) > abs(u[1])*ssreltol)
    #test2 = (abs(du[2]) > ssabstol) && (abs(du[2]) > abs(u[2])*ssreltol)
    #if !any([test1, test2])
    #   display([t, "terminate integration"])
    #   sleep(20)
    #end 
end

function conditiontemp(u, t, integrator, abstol, reltol)

    if DiffEqBase.isinplace(integrator.sol.prob)
        testval = first(get_tmp_cache(integrator))
        DiffEqBase.get_du!(testval, integrator)
        if typeof(integrator.sol.prob) <: DiffEqBase.DiscreteProblem
            @. testval =  testval - integrator.u
        end
    else
        testval = get_du(integrator)
        if typeof(integrator.sol.prob) <: DiffEqBase.DiscreteProblem
            testval =  testval - integrator.u
        end
    end

    any(abs(d) > abstol && abs(d) > reltol*abs(u) for (d,abstol, reltol, u) =
           zip(testval, Iterators.cycle(abstol), Iterators.cycle(reltol), integrator.u)) && (return false)
    return true
end

function affect!(integrator)
    terminate!(integrator)
end

function main()
    k1f = 1000000
    K1  = 30
    k1r = k1f/K1
    k2f = 1e-20
    K2  = 1e-08
    k2r = k2f/K2
    k3f = 1e-20
    K3  = 1
    k3r = k3f/K3
    ssabstol = 1e-8
    ssreltol = 1e-8

    u0 = [0;0]
    p         = [k1f; k1r; k2f; k2r; k3f; k3r;ssabstol;ssreltol]
    tspan     = (0, 1e15)
    condition = (u, t, integrator) -> conditiontemp(u, t, integrator, ssabstol, ssreltol)
    cb        = DiscreteCallback(condition, affect!)
    prob1     = ODEProblem(ODESYSTEM, u0, tspan, p)

    #Solve for steady state by integrating to t --> inf
    sol1      = solve(prob1, Kvaerno3(), abstol = 1e-8, reltol = 1e-8)

    #Solve for steady state using manual implementation of TerminateSteadyState callback
    sol2      = solve(prob1, Kvaerno3(), abstol = 1e-8, reltol = 1e-8, callback = cb)

    #solve for steady state using DynamicSS solver
    prob3 = SteadyStateProblem(ODESYSTEM, u0, p)
    sol3  = solve(prob3, DynamicSS(Kvaerno3()), abstol = ssabstol, reltol = ssreltol) 

    #had a quirk where I couldn't plot on log10 x scale b/c of 0 initial value. hacky workaround
    tset = exp10.(range(-16, stop=15, length=1000))
    fine = sol1(tset)'
    thetaA = fine[:,1]
    thetaB = fine[:,2]
    #thetaV = 1 .- thetaA .- thetaB

    #Store true steady state solution from integrating to t --> inf
    ssvals = [thetaA[end]; thetaB[end]]

    #create plots of variables to illustrate time scales for reaching steady state
    p1   = plot(tset, thetaA, xscale = :log10, xlabel = "time", ylabel = "thetaA", ylim = (0,1), xlim = (1e-16, 1e15), legend = false)   
    p2   = plot(tset, thetaB, xscale = :log10, xlabel = "time", ylabel = "thetaB", ylim = (0,1.5*maximum(thetaB)), xlim = (1e-16, 1e15), legend = false)
    #p3   = plot(tset, thetaV, xscale = :log10, xlabel = "time", ylabel = "thetaV", ylim = (0,1), xlim = (1e-16, 1e15), legend = false)

    #display plots
    display(plot(p1, p2, layout = (2,1), legend = false))
    
    #display true steady state vs. two obtained using callbacks to terminate integration.
    display([ssvals sol2[end] sol3])
end

main()