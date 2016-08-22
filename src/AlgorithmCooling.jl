export DemonCooling,cool!,heat!,trotter,propagate

immutable DemonCooling <: QuDynamics.QuPropagatorMethod
    t::Real
    gamma::Real
end

# auto-set parameters

function DemonCooling(H::AbstractMatrix)
    maxEigen = eigs(H;nev=1,which=:LR)[1][1]|>real
    minEigen = eigs(H;nev=1,which=:SR)[1][1]|>real

    m = Model(solver=IpoptSolver())

    @variable(m, t >= 0)
    @variable(m, γ >= 0)
    # 0 <= γ <= (minEigen/(maxEigen-minEigen)+0.5)*π

    μ = (maxEigen+minEigen)/2
    @NLobjective(m, Max, ((1-sin(minEigen*t-γ))*(1-sin(γ)))/(1-sin(μ*t-γ)))
    @constraint(m, -π/2<=γ-minEigen*t<=π/2)
    @constraint(m, -π/2<=γ-maxEigen*t<=π/2)
    status = solve(m)

    @show getvalue(t)
    @show getvalue(γ)
    DemonCooling(getvalue(t),getvalue(γ))
end

#trotter expansion
function trotter(A::AbstractMatrix,B::AbstractMatrix,P::Int64)
    return (expm(full(A/(2*P)))*expm(full(B/P))*expm(full(A/(2*P))))^P
end


function propagate(prob::DemonCooling, eq::AQCShrodingerEq, t, current_t, current_qustate)
    s = eq.t/eq.T

    # Use trotter expansion to simulate experimental implementation
    # P = 3 recommended!
    U = trotter(-im*prob.t*(1-s)*coeffs(eq.HB),
                    -im*prob.t*s*coeffs(eq.HP),3)

    # @show "pass"
    #
    # U =expm((-im*((1-s)*coeffs(eq.HB)+s*coeffs(eq.HP))*prob.t)|>full)

    p,next_state = cooling!(current_qustate,U,prob.t,prob.gamma)
    eq.p *= p
    return next_state
end

# for general use
function cooling!(state::AbstractQuVector,U::AbstractMatrix,t::Real,gamma::Real)
    dice = rand()
    if dice <= 0.5*(1+sin(gamma))
        # println("heating")
        heat!(state,U,t,gamma)
        return 0.5*(1+sin(gamma)),state
    else
        # println("cooling")
        cool!(state,U,t,gamma)
        return 0.5*(1-sin(gamma)),state
    end
end

function heat!(state::AbstractQuVector,U::AbstractMatrix,t::Real,gamma::Real)
    next_state = 0.5*( coeffs(state)+im*exp(im*gamma)*U*coeffs(state) ) |> Base.normalize!
    state.coeffs = next_state
    return state
end

function cool!(state::AbstractQuVector,U::AbstractMatrix,t::Real,gamma::Real)
    next_state = 0.5*( coeffs(state)-im*exp(im*gamma)*U*coeffs(state) ) |> Base.normalize!
    state.coeffs = next_state
    return state
end

function set_cooling_parameters(HB::AbstractMatrix,BP::AbstractMatrix,s::Real)
    H = (1-s)*HB+s*HP
    maxEigen = eigs(H;nev=1,which=:LR)[1][1]|>real
    minEigen = eigs(H;nev=1,which=:SR)[1][1]|>real

    m = Model(solver=IpoptSolver())

    @variable(m, t >= 0)
    @variable(m, γ >= 0)
    # 0 <= γ <= (minEigen/(maxEigen-minEigen)+0.5)*π

    μ = (maxEigen+minEigen)/2
    @NLobjective(m, Max, ((1-sin(minEigen*t-γ))*(1-sin(γ)))/(1-sin(μ*t-γ)))
    @constraint(m, -π/2<=γ-minEigen*t<=π/2)
    @constraint(m, -π/2<=γ-maxEigen*t<=π/2)
    status = solve(m)

    @show getvalue(γ)

    return getvalue(t),getvalue(γ)
end
