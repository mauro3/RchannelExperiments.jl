export calculate_press_grad, calculate_press_grad_all

"""
    calculate_press_grad(res, port1, port2)

Calculates the pressure gradient over the ports:


"""
calculate_press_grad(res::RcExp, port1, port2) =
    calculate_press_grad(res.ex, res.resl)
function calculate_press_grad(ex::RcExp, resl::ExpLabViewResults, port1, port2)
    @unpack portdists, port2channels = ex.setup
    @unpack press, press_std = resl
    #press, press_std = resl.press_raw, resl.press_std_raw
    ports = keys(port2channels)
    slot1, slot2 = find(x->x==port1, ports)[1], find(x->x==port2, ports)[1]
    deltaL = sum(portdists[slot1:slot2-1])
    grad = (press[port2]-press[port1])./deltaL
    grad_std = sqrt.(press_std[port2].^2 + press_std[port1].^2)./deltaL
    return grad, grad_std, (port1=>port2, deltaL)
end
calculate_press_grad(res::RcExp) = calculate_press_grad(res.ex, res.resl)
function calculate_press_grad(ex::RcExp, resl::ExpLabViewResults)
    @unpack port2channels = ex.setup
    # select ice1, ice2 or in4
    nums = [8,9,4:-1:1;]
    port1 = 0
    for n in nums
        if haskey(port2channels, Port.PressPorts(n))
            port1 = Port.PressPorts(n)
            break
        end
    end
    # select ice4, ice3 or out1
    nums = [11,10,4:7;]
    port2 = 0
    for n in nums
        if haskey(port2channels, Port.PressPorts(n))
            port2 = Port.PressPorts(n)
            break
        end
    end
    calculate_press_grad(res, port1, port2)
end

calculate_press_grad_all(res::RcRes) = calculate_press_grad_all(res.ex, res.resl)
function calculate_press_grad_all(ex::RcExp, resl::ExpLabViewResults,
                                  ports1 = collect(keys(ex.setup.port2channels))[1:end-1],
                                  ports2 = collect(keys(ex.setup.port2channels))[2:end])
    grads = Float64[]
    grads_std = Float64[]
    ports = []
    for (p1,p2) in zip(ports1, ports2)
        pd, spd, ps = calculate_press_grad(ex, resl, p1, p2)
        append!(grads, pd)
        append!(grads_std, spd)
        push!(ports, ps)
    end
    nl = length(ports1)
    return reshape(grads, length(grads)÷nl, nl),
           reshape(grads_std, length(grads_std)÷nl, nl),
           ports
end


calc_f(grad_p, D, Q) =  (grad_p * pi^2 * D^5)/(8*Q^2*rho)
function calc_f(res::RcRes, port1, port2, inds=nothing)
    @unpack resl, resi = res
    @unpack trun, irun = resl
    if inds!=nothing
        irun = irun[inds]
        trun = trun[inds]
    end
    grad_p, grad_p_std, ports = calculate_press_grad(res, port1, port2)
    grad_p, grad_p_std = grad_p[irun], grad_p_std[irun]
    Q = resl.Q[irun]
    D, D_std = get_time_series(resi, trun)
    calc_f.(grad_p, D, Q)
end


function Re_nr(D,Q)
    A = pi/4*D^2
    D*Q/(nu*A)
end
function Re_nr(res::RcRes, inds=nothing)
    @unpack resl, resi = res
    @unpack trun, irun = resl
    if inds!=nothing
        irun = irun[inds]
        trun = trun[inds]
    end
    D, D_std = get_time_series(resi, trun)
    Q = resl.Q[irun]
    Re_nr.(D,Q)
end
