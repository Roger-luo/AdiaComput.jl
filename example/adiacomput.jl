using AdiaComput,CUDArt

blas_set_num_threads(4)


function adia(n::Int)
    ins,ans = generate(n)
    q = AdiaComputer(ins,n,1;GPU=true)
    @time q|>next_timestep!|>next_timestep!|>next_timestep!
    q.prob*=norm(to_host(q.state)[ans[1]+1])^2
    print("\033[34;1m normal $n bits: \033[31;1m $(q.prob)\n \033[0m")

    # prob = 0
    cq = AdiaComputer(ins,n,1;GPU=true)

    prob = @time @parallel (+) for i=1:100
        print("\r\033[34;1m cooling $n Bits:\033[31;1m $i%\033[0m")
        cq = next_timestep!(cq)
        cq = cooling!(cq)
        cq = next_timestep!(cq)

        cq = cooling!(cq)
        cq = next_timestep!(cq)
        cq.prob*=norm(to_host(cq.state)[ans[1]+1])^2
        res = cq.prob
        reset!(cq)
        res
    end
    print("\033[34;1m cooled $n bits: \033[31;1m $(prob/100)\033[0m\n")
end

adia(5)
