export DemonCooling,cool!,heat!,trotter,propagate

immutable DemonCooling <: QuDynamics.QuPropagatorMethod
    t::Real
    gamma::Real
end

# auto-set parameters

function DemonCooling(H::AbstractMatrix)
    maxEigen = eigs(H;nev=1,which=:LR)[1][1]|>real
    # minEigen = eigs(H;nev=1,which=:SR)[1][1]|>real

    # m = Model(solver=IpoptSolver())
    #
    # @variable(m, t >= 0)
    # @variable(m, γ >= 0)
    # 0 <= γ <= (minEigen/(maxEigen-minEigen)+0.5)*π

    # μ = (maxEigen+minEigen)/2
    # @NLobjective(m, Max, ((1-sin(minEigen*t-γ)))/(1-sin(μ*t-γ)))
    # @constraint(m, -π/2<=γ-minEigen*t<=π/2)
    # @constraint(m, -π/2<=γ-maxEigen*t<=π/2)
    # status = solve(m)
    #
    # @show getvalue(t)
    # @show getvalue(γ)
    DemonCooling(π/(2*maxEigen),0)
end

#trotter expansion
function trotter(A::AbstractMatrix,B::AbstractMatrix,P::Int64)
    return (expm(full(A/(2*P)))*expm(full(B/P))*expm(full(A/(2*P))))^P
end


function propagate(prob::DemonCooling, eq::AQCShrodingerEq, t, current_t, current_qustate)
    s = eq.t/eq.T

    # Use trotter expansion to simulate experimental implementation
    # P = 3 recommended!
    # U = trotter(-im*prob.t*(1-s)*coeffs(eq.HB),
    #                 -im*prob.t*s*coeffs(eq.HP),3)

    H = (1-s)*eq.HB+s*eq.HP|>coeffs

    # @show "pass"
    #
    # U =expm((-im*((1-s)*coeffs(eq.HB)+s*coeffs(eq.HP))*prob.t)|>full)

    p,next_state = cooling!(current_qustate,H,prob.t,prob.gamma)
    eq.p *= p
    return next_state
end

function probability(H::AbstractMatrix,state::AbstractVector,t::Real,gamma::Real)
    c_ret = 0
    eigens = (H|>full|>eigfact)

    for k = 1:length(state)
        c_ret += abs2(state⋅eigens[:vectors][:,k])*(1-sin((eigens[:values][k]*t-gamma)|>abs))
    end

    h_ret = 0

    for k = 1:length(state)
        h_ret += abs2(state⋅eigens[:vectors][:,k])*(1+sin((eigens[:values][k]*t-gamma)|>abs))
    end

    sum = c_ret+h_ret
    c_ret = c_ret/sum
    h_ret = h_ret/sum
    return c_ret,h_ret
end

# for general use
function cooling!(state::AbstractQuVector,H::AbstractMatrix,t::Real,gamma::Real)
    dice = rand()

    U = expm(-im*t*H|>full)

    p_cool,p_heat = probability(H,coeffs(state),t,gamma)

    if dice <= p_heat
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
    minEigen = eigs(H;nev=1,which=:SR)[1][1]|>real

    eigens = full(H)|>eigfact

    m = Model(solver=IpoptSolver())

    @variable(m, t >= 0)
    @variable(m, γ >= 0)
    # 0 <= γ <= (minEigen/(maxEigen-minEigen)+0.5)*π

    μ = (maxEigen+minEigen)/2
    @NLobjective(m, Max, (1-sin(minEigen*t-γ))/(sum(x->(1-sin(x*t-γ)))))
    # @NLobjective(m, Max, ((1-sin(minEigen*t-γ))*(1-sin(γ)))/(1-sin(μ*t-γ)))
    @constraint(m, -π/2<=γ-minEigen*t<=π/2)
    @constraint(m, -π/2<=γ-maxEigen*t<=π/2)
    status = solve(m)

    @show getvalue(γ)

    return getvalue(t),getvalue(γ)
end
