using DrWatson
@quickactivate "TrajectoryOptimizationScripts"

include(srcdir("solveCRTBPTrajectoryOptimization.jl"))
using AstroEOMs, AstroUtils, DifferentialEquations, StaticArrays, DelimitedFiles

function main()
    # ===== Spacecraft parameters
    m0      = 1500.0    # Initial mass [kg]
    mp      = 1200.0    # Propellant mass [kg]
    tMax    = 10.0      # Max thrust [N]
    Isp     = 3000.0    # Specific impulse [s]

    # ===== Dynamical model parameters (Earth Moon system)
    me      = 5.9742e24     # Mass of Earth [kg]
    mm      = 7.3483e22     # Mass of Moon [kg]
    Re      = 6378.1        # Radius of Earth [km]
    Rm      = 1737.4        # Radius of Moon [km]
    MU      = m0            # Mass unit (setting equal to spacecraft init. mass)
    r12     = 384400.0      # Distance between Earth and Moon [km]

    # ===== Problem parameters
    # Constrained time span
    tspan = (0.0, 8.6404) # [Days]

    # Constrained initial and final states
    ics = [-0.019488511458668, -0.016033479812051, 0.0, # Position 
            8.918881923678198, -4.081793688818725, 0.0, # Velocity
            1.0] # Mass
    fcs = [ 0.823385182067467, 0.0, -0.022277556273235, # Position
            0.0, 0.134184170262437, 0.0, # Velocity
            0.0] # Mass co-state

    # ===== Construct problem parameters
    spaceCraft  = SimpleSpacecraft(m0,mp,tMax,Isp)
    crtbpParams = CRTBPParams(me,mm,Re,Rm,MU,r12)
    sw          = SwitchingStruct(1.0, 0)

    # ===== Guess for unknown initial co-states (for minimum energy problem with Tmax = 10 N)
    λ0g = [18.061015041426327,37.62023638562631,-0.12141606546723561,-0.11720442210452124,0.05108452903629872,-0.00011125433108675963,0.15584636103201133]

    # ==== Solve minimum energy problem first for decreasing values of thrust
    n       = 1000
    tMax_s  = 10.0
    tMax_f  = 0.6
    tMaxs   = range(start = tMax_s, stop = tMax_f, length = n)

    # Curve fit for tMax
    a       = 306.6
    b       = -1.83
    c       = 42.03
    d       = -0.1711
    tfs     = a*exp.(b*tMaxs) .+ c*exp.(d*tMaxs)
    for i in eachindex(tMaxs)
        display(tMaxs[i])
        spaceCraft      = SimpleSpacecraft(m0,mp,tMaxs[i],Isp)
        tspan           = (0.0, tfs[i])
        λme,flag        = solveCRTBPTrajectoryOptimization(λ0g,ics,fcs,tspan,(spaceCraft,crtbpParams,sw))
        λ0g .= λme
    end

    # ==== Begin minimum energy-to-fuel homotopy of previous
    #      solve was successful
    if flag == false
        print("First solve failed to converge. Not performing homotopy continuation.\n")
    else
        n   = 25 # Number of steps in homotopy
        ϵs  = [(j^2 - 1.0)/(n^2 - 1.0) for j in n:-1:1]
        λ0  = copy(λ0g)
        for i in eachindex(ϵs)
            # Set homotopy parameter
            sw.ϵ = ϵs[i]

            # Solve problem with updated homotopy parameter
            λ0n, flag = solveCRTBPTrajectoryOptimization(λ0,ics,fcs,tspan,(spaceCraft,crtbpParams,sw))

            # Check to see if we converged
            if flag == false
                print("Intermediate solve failed during homotopy continuation. Breaking from continuation loop.")

                # Set epsilon to value corresponding to most recent solution
                sw.ϵ = ϵs[i - 1]
                break
            else
                # Update initial costates
                λ0 .= λ0n
            end
        end
        
        # Print message if homotopy continuation was successful
        if sw.ϵ == 0.0
            print("Homotopy continuation was successful at computing minimum fuel solution.")
        end

        # ===== Perform numerical integration and save full trajectory to file
        # Set thrust state
        csc = spaceCraft.c * crtbpParams.TU / (1000.0 * crtbpParams.LU)
        λv  = norm(view(λ0, 4:6))
        S   = -λv*csc / ics[7] - λ0[7] + 1.0
        sw.state = S > 0.0 ? 0 : 2

        # Construct requriements
        ts  = (tspan[1]*86400.0 / crtbpParams.TU, tspan[2]*86400.0 / crtbpParams.TU)
        x0  = SVector(ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[7],
                        λ0[1], λ0[2], λ0[3], λ0[4], λ0[5], λ0[6], λ0[7])
        cb  = ContinuousCallback(AstroEOMs.crtbpMEMFSwitchingCallbackCondition,
                                AstroEOMs.crtbpMEMFSwitchingCallbackAffect!;
                                interp_points = 25)
        p   = ODEProblem(ODEFunction{false}(AstroEOMs.crtbpMinimumFuelOptimalControlMEMFEOMs),
                    x0, ts, (spaceCraft, crtbpParams, sw))

        # Solve 
        sol = solve(p, Vern9(), reltol = 1e-14, abstol = 1e-14, callback = cb)

        # Grab interpolated states
        n   = 5000 # Number of points 
        tv  = range(start = ts[1], stop = ts[2], length = n)
        data = zeros(n, 15)
        for i in eachindex(tv)
            data[i,1]       = tv[i]
            data[i,2:end]   .= sol(tv[i])
        end

        # Write to file
        open(datadir("solution.txt"), "w") do io; writedlm(io, data); end

        # ===== Integrate CRTBP EOMs with Constrained Final State to Obtain Halo
        tsh = (0.0, 12.0*86400.0 / crtbpParams.TU)
        p   = ODEProblem(ODEFunction{false}(AstroEOMs.crtbpEomNoControl),
                SVector(fcs[1:6]...), tsh, crtbpParams)
        sol = solve(p, Vern9(), reltol = 1e-14, abstol = 1e-14)

        # Grab interpolated states
        n   = 500
        tv  = range(start = tsh[1], stop = tsh[2], length = n)
        hdata = zeros(n, 6)
        for i in eachindex(tv)
            hdata[i,:] .= sol(tv[i])
        end

        # Write to file
        open(datadir("halo.txt"), "w") do io; writedlm(io, hdata); end
    end
end

main()