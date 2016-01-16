# const T             = 1e3;
# const initial_state = [1/sqrt(2),1/sqrt(2)]

# using PyPlot

# function Hamilton(s)
#     @assert 0<=s<=1
#     return 0.5*[1+s s-1;s-1 1-s]
# end

# function next_timestep(state::AbstractVector,step::Int64,dt=1e-2)
#     return expm(-im*dt*Hamilton(step*dt/T))*state
# end

# function evolution(dt=1)
#     state = initial_state
#     eigenvalue = collect(complex(conj(state).'*Hamilton(0)*state))
#     for i=1:1e3
#         state = next_timestep(state,Int64(i),dt)
#         append!(eigenvalue, conj(state).'*Hamilton(i*dt/T)*state)
#     end
#     return eigenvalue,state
# end

# # @show Hamilton(1)*Hamilton(1)

# res,state = evolution()
# @show state


# plot(real(res))
# show()

blas_set_num_threads(4)
A=rand(1000,1000)
B=rand(1000,1000)
@time for i=1:100 A*B end