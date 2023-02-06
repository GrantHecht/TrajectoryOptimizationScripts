using NLsolve, DifferentialEquations, ForwardDiff
using StaticArrays, AstroEOMs, DiffEqSensitivity
using LinearAlgebra

function solveCRTBPTrajectoryOptimization(λ0,ics,fcs,tspan,ps; acceptanceTol = 10.0)
    # Evaluate shooting function once to test is numerical integration is successful with guess
    Ft      = zeros(7)
    retcode = shootingFunction!(Ft,λ0, ics, fcs, tspan, ps)

    if retcode == :Success && maximum(Ft) < acceptanceTol
        # Solve the boundary value problem with NLsolve
        sol = nlsolve((F,x) -> shootingFunction!(F,x,ics,fcs,tspan,ps),
                λ0, autodiff = :forward, show_trace = true, ftol = 1e-8)

        # Return nlsolve solution
        return (sol.zero, converged(sol))
    else
        print("Numerical integration with guess failed. Not attempting to solve!\n")
        return (zeros(7), false)
    end
end

function shootingFunction!(F,λ0,ics,fcs,tspan,ps)
    # Set state of switching struct
    csc = ps[1].c * ps[2].TU / (1000.0 * ps[2].LU)
    λv  = norm(view(λ0, 4:6))
    S   = -λv*csc / ics[7] - λ0[7] + 1.0
    if S > ps[3].ϵ
        ps[3].state = 0
    elseif S < -ps[3].ϵ
        ps[3].state = 2
    else
        ps[3].state = 1
    end

    # Form initial state vector
    x0  = SVector(ics[1],ics[2],ics[3],ics[4],ics[5],ics[6],ics[7],
            λ0[1],λ0[2],λ0[3],λ0[4],λ0[5],λ0[6],λ0[7])

    # Compute scalled time span
    ts  = (tspan[1]*86400/ps[2].TU, tspan[2]*86400/ps[2].TU)

    # Construct callback
    cb  = ContinuousCallback(AstroEOMs.crtbpMEMFSwitchingCallbackCondition,
                             AstroEOMs.crtbpMEMFSwitchingCallbackAffect!;
                             save_positions = (false,false),
                             interp_points  = 25)

    # Perform numerical integration
    prob = ODEProblem(ODEFunction{false}(AstroEOMs.crtbpMinimumFuelOptimalControlMEMFEOMs),
        x0, ts, ps)
    sol  = solve(prob, Vern9(), reltol = 1e-14, abstol = 1e-14, 
            save_everystep = false, save_start = false, initialize_save = false,
            maxiters = 1e7, callback = cb, 
            sensealg = ForwardDiffSensitivity(convert_tspan=true))
    xf   = sol[end]

    # Compute residuals
    @views F[1:6] .= xf[1:6] .- fcs[1:6]
    F[7] = xf[14] - fcs[7]

    return sol.retcode
end