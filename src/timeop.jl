function realtimeop!(Hs::AdiaComputer)
    HP = (-im*Hs.location/3)*Hs.HP
    HB = (-im*(1-Hs.location)/6)*Hs.HB

    timemat = Hs.P*diagexp(HB)*Hs.invP
    timemat = timemat*diagexp(HP)*timemat

    timemat = timemat^3
    Hs.state = timemat*Hs.state
    Hs.location += Hs.dt/Hs.maxtime
end

function next_timestep!(Hs::AdiaComputer;evopercentage::Real=1/3,nev=6)
    @assert (0<=Hs.location+evopercentage)&&( Hs.location+evopercentage-1<0.1) "evolutoin percentage out of bounds(should be in [0,1])"

    const evotime = Hs.maxtime
    const dt      = Hs.dt

    for i=Hs.location*evotime:dt:(Hs.location+evopercentage)*evotime
        realtimeop!(Hs)
    end
    return Hs
end

export next_timestep!,SimilarMatrix
