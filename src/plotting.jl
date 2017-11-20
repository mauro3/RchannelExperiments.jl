
function plot_grad(res::RcRes; fig_axs=P.subplots(1,1),
                   t0=[:exp, :trun, :time, :ns][2])
    ex = res.ex
    fig, axs = fig_axs



end

function plot_res(res::RcRes; fig_axs=P.subplots(6,1,sharex="all"),
                  t0=[:exp, :trun, :time, :ns][2],
                  kwargs...)
    ex = res.ex
    fig, axs = fig_axs

    plot_res(ex, res.resl, fig_axs=(fig,axs[1:2]); t0=t0, kwargs...)

    plot_res(ex, res.resi, fig_axs=(fig,axs[4:6]); t0=res.resl.trun[1:100:end], kwargs...)
end

function plot_moody_chart(Re=10.^(log10(4000):0.01:5), rel_roughness=10.^(-4:0.5:0.1);
                          ax=P.subplots(1,1)[2],  setticks=true, loglog=true)
    P.axes(ax)
    minmax = [Inf,0.0]
    for r in rel_roughness
        v = f_haaland.(Re,r)
        minmax[1] = min(minmax[1], minimum(v))
        minmax[2] = max(minmax[2], maximum(v))
        if loglog
            P.loglog(Re, v, "k:", alpha=0.5)
        else
            P.plot(Re, v, "k:", alpha=0.5)
        end
    end
    if setticks
        xt = [1000, 4000, 10000, 40000, 10^5, 4*10^5, 10^6]
        xti = find(x->Re[1]<=x+10, xt)[1]:find(x->Re[end]<=x+10, xt)[1]
        yt = [0.002:0.002:0.01, 0.02:0.02:0.1, 0.2:0.2:1;]
        yti = find(y->y>=minmax[1], yt)[1]:find(y->y>minmax[2], yt)[1]
        P.xticks(xt[xti],xt[xti])
        P.yticks(yt[yti],yt[yti])
    end
    P.xlabel("Re")
    P.ylabel("f")
end


function proc_for_f(rr, inds)
    mm = 1e-3
    l = 1e-3
    D_rel = 0.02
    D_abs = 1mm
    Q_rel = 0.01
    Q_abs = 0.1l
    p_abs = RcLabView.m_h2o_to_Pa(0.00)
    p_std_fac = 1
    p_rel = 0.01
    d_deltaL = 5mm
    stacksize = 100


    port1, port2, ind, remark = inds
    @unpack Q, press, press_std, trun, irun  = rr.resl
    ind4orig = irun[ind]
    ii = find_low_std(press_std, [port1, port2], ind4orig);
    ind4orig = ind4orig[ii]
    trun = rr.resl.t[ind4orig]

    dia_mean, dia_std, dia_med, dia_10, dia_90, scalloping, sym =
        get_time_series(rr.resi, rr.resl.t[ind4orig])

    D, Q, p1, p2, grad_p, tt = add_uncertainty_stack_first(rr, port1, port2,
                                                           D_rel, D_abs,
                                                           Q_rel, Q_abs,
                                                           p_abs, p_std_fac, p_rel,
                                                           d_deltaL,
                                                           stacksize,
                                                           ind4orig)
    scalloping = stack_it(scalloping, stacksize)


    f = calc_f.(grad_p, D, Q);
    Re = Re_nr.(D,Q);
    Rev = Measurements.value.(Re)
    fv = Measurements.value.(f)

    return Re, f, Rev, fv, scalloping, tt
end
function plot_res_moody(rr, inds; ax=P.subplots(1,1)[2], colorbar=true,
                        colors=[:scalopping,:time, :none][2],
                        label=rr.ex.name,
                        fmin = 0.002,
                        fmax = 0.15,
                        Remin = 35000,
                        Remax = 80000)


    Re, f, Rev, fv, scalloping = proc_for_f(rr, inds)

    colvar,clab = if colors==:scalopping
        (scalloping, "Scalloping ()")
    elseif colors==:time
        (tt/60, "Run-time (min)")
    elseif colors==:none
        (1,1)
    else
        error()
    end



    plot_moody_chart(logspace(log10(Remin), log10(Remax), 1000),
                       10.^(-5:0.5:-0.5),
                       ax=ax, setticks=false, loglog=false)

    P.axes(ax)
    ax[:semilogy]()
    P.errorbar(Re,f,fmt="k.")
    if colors==:none
        P.scatter(Rev, fv, 200*ones(fv), label=label)
    else
        P.scatter(Rev, fv, 200*ones(fv), colvar)
    end
    if colorbar && colors!=:none
        cb = P.colorbar()
        cb[:set_label](clab)
    end
    P.xlabel("Re")
    P.ylabel("f")
    P.title(rr.ex.name)


    P.ylim((fmin,fmax))
    P.xlim((Remin,Remax))
    xt = [1000, 4000, 10000, 40000, 10^5, 4*10^5, 10^6]
    xti = find(x->Remin<=x+10, xt)[1]:find(x->Remax<=x+10, xt)[1]
    yt = [0.002:0.002:0.01, 0.02:0.02:0.1, 0.2:0.2:1;]
    yti = find(y->y>=fmin, yt)[1]:find(y->y>fmax, yt)[1]
    # ax[:set_xticks](xt[xti], string.(xt[xti]))
    # ax[:set_yticks](yt[yti], string.(yt[yti]))
    # P.xticks(xt[xti], string.(xt[xti]))
    P.yticks(yt[yti], string.(yt[yti]))

end

function plot_res_scatter_f_scalops(rr, inds; ax=P.subplots(1,1)[2],
                                    colorbar=true, label="")
    Re, f, Rev, fv, scalloping, tt = proc_for_f(rr, inds)
    P.axes(ax)
    P.semilogy(scalloping, fv, ".", label=label)
    # P.xlabel("Scalloping ()")
    # P.ylabel("f ()")
    P.xlabel("Scalloping")
    P.ylabel("f")
    #P.scatter(scalloping, fv, 500*ones(fv), tt/60)
    # if colorbar
    #     cb = P.colorbar()
    #     cb[:set_label]("Time (min)")
    # end
end
