
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
