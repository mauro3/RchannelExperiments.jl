export calculate_press_grad, calculate_press_grad_all

"""
    port_dist(port1, port2, port2channels, portdists)
    port_dist(port1, port2, res)

Distance between two pressure ports
"""
port_dist(port1, port2, res::RcRes) = port_dist(port1, port2, res.ex)
function port_dist(port1, port2, ex::RcExp)
    @unpack portdists, port2channels = ex.setup
    port_dist(port1, port2, port2channels, portdists)
end
function port_dist(port1, port2, port2channels, portdists)
    ports = keys(port2channels)
    slot1, slot2 = find(x->x==port1, ports)[1], find(x->x==port2, ports)[1]
    deltaL = sum(portdists[slot1:slot2-1])
    return deltaL
end


"""
    calculate_press_grad(res, port1, port2)

Calculates the pressure gradient over the ports:


"""
calculate_press_grad(res::RcRes, port1, port2) =
    calculate_press_grad(res.ex, res.resl)
function calculate_press_grad(ex::RcExp, resl::ExpLabViewResults, port1, port2)
    @unpack press, press_std = resl
    #press, press_std = resl.press_raw, resl.press_std_raw
    deltaL = port_dist(port1, port2, ex)
    grad, grad_std = calculate_press_grad(press[port1], press[port2], press_std[port1], press_std[port2], deltaL)
    return grad, grad_std, (port1=>port2, deltaL)
end
function calculate_press_grad(press1, press2, press1_std, press2_std, deltaL)
    grad = (press1-press2)./deltaL
    grad_std = sqrt.(press1_std.^2 + press2_std.^2)./deltaL
    return grad, grad_std
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


"""
    f_haaland(Re, epsilon, Dh)
    f_haaland(Re, rel_roughness)

D-W f factor by Haaland 1983.  Epsilon is the size of the roughness,
Dh hyd-diameter.
"""
f_haaland(Re, epsilon, Dh) = f_haaland(Re, epsilon/Dh)
function f_haaland(Re, rel_roughness)
    tmp = -1.8 * log( ((rel_roughness/3.7)^1.11 + 6.9/Re) )
    1/tmp^2
end

### stacking
function stack_it(vec, step)
    halfstep = step÷2
    out = eltype(vec)[]
    for i = halfstep+1:step:length(vec)-halfstep
        push!(out, mean(vec[i-halfstep:i+halfstep]))
    end
    out
end
function stack_it_time(time, step)
    halfstep = step÷2
    out = eltype(time)[]
    for i = halfstep+1:step:length(time)-halfstep
        push!(out, time[i])
    end
    out
end


### purge high std
"""
    find_low_std(press_std, ports, irun=1:length(first(press_std)), thresh=100)
    -> inds

Find indices of low std.
"""
function find_low_std(press_std, ports, irun=1:length(first(press_std)), thresh=100)
    out = press_std[ports[1]][irun].<thresh
    for i=2:length(ports)
        out = (press_std[ports[i]][irun].<thresh) .& out
    end
    find(out)
end
