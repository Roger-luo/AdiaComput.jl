function cooler!(Hs::AdiaComputer,gamma::Real,t::Real)
    HB = (-im*t*(1-Hs.location)/6)*Hs.HB
    HP = (-im*t*Hs.location/3)*Hs.HP
    temp_mat = Hs.P*diagexp(HB)*Hs.invP
    temp_mat = temp_mat*diagexp(HP)*temp_mat
    temp_mat = temp_mat^3


    Hs.state.-((im*exp(im*gamma))*temp_mat*Hs.state)
    Hs.state.*Complex128(0.5)
    normalize!(Hs.state)
end

function heater!(Hs::AdiaComputer,gamma::Real,t::Real)
    HB = (-im*t*(1-Hs.location)/6)*Hs.HB
    HP = (-im*t*Hs.location/3)*Hs.HP
    temp_mat = Hs.P*diagexp(HB)*Hs.invP
    temp_mat = temp_mat*diagexp(HP)*temp_mat
    temp_mat = temp_mat^3

    Hs.state.+((im*exp(im*gamma))*temp_mat*Hs.state)
    Hs.state.*Complex128(0.5)
    normalize!(Hs.state)
end

function daemon!(
    Hs::AdiaComputer,
    gamma::Real,
    t::Real
    )
    dice = rand()
    if dice <= 0.5*(1+sin(gamma))
        #get |1> (higher energy)
        heater!(Hs,gamma,t)
        return 0.5*(1-sin(gamma))
    else
        #get |0> (lower energy)
        cooler!(Hs,gamma,t)
        return 0.5*(1+sin(gamma))
    end
end

function cooling!(
    Hs::AdiaComputer;
    n=5)

    count = 0
    gamma,t = CoolingPara(Hs)

    while count<n
        @show count
        daemon!(Hs,gamma,t)
        count += 1
    end
    return Hs
end

# function cost_func(t,gamma)

function CoolingPara(Hs::AdiaComputer)
    H = Hamiltonian(Hs)
    maxEigen = eigs(Hamiltonian(Hs);nev=1,which=:LR)[1][1]|>real
    minEigen = eigs(Hamiltonian(Hs);nev=1,which=:SR)[1][1]|>real

    gamma = (maxEigen+minEigen)/(maxEigen-minEigen) * pi/2 * 0.1

    if (gamma-pi/2)>0
        t = 0.5*((gamma-pi/2)/minEigen + (gamma+pi/2)/maxEigen)
    else
        t = 0.5*((gamma-pi/2)/maxEigen + (gamma+pi/2)/maxEigen)
    end
    return gamma,t
end

export cooling!
