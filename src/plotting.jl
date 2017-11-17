
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
                          ax=P.subplots(1,1)[2],  setticks=true)
    P.axes(ax)
    minmax = [Inf,0.0]
    for r in rel_roughness
        v = f_haaland.(Re,r)
        minmax[1] = min(minmax[1], minimum(v))
        minmax[2] = max(minmax[2], maximum(v))
        P.loglog(Re, v, "k:", alpha=0.5)
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
