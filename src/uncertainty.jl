using Measurements

function add_uncertainty_stack_first(res, port1, port2,
                                     D_rel, D_abs,
                                     Q_rel, Q_abs,
                                     p_abs, p_std_fac, p_rel,
                                     d_deltaL,
                                     stacksize, inds=res.resl.irun)
    @unpack Q, press, press_std, trun  = res.resl

    ii = R.find_low_std(press_std, [port1, port2], inds);
    trun = trun[ii];
    inds = inds[ii];

    D = R.get_time_series(res.resi, trun)[1]
    Q = Q[inds]
    p1, p2 = press[port1][inds], press[port2][inds]
    p1_std, p2_std = press_std[port1][inds], press[port2][inds]
    deltaL = R.port_dist(port1, port2, res)

    D, Q, p1, p2, grad_p = add_uncertainty(D, Q, p1, p2, deltaL,
                                           D_rel, D_abs,
                                           Q_rel, Q_abs,
                                           p_abs,
                                           p1_std*p_std_fac, p2_std*p_std_fac,
                                           p_rel)
    D, Q, p1, p2, p1_std, p2_std = [R.stack_it(vec,stacksize) for vec in (D,Q,p1,p2,p1_std,p2_std)];
    tstac = R.stack_it_time(trun, stacksize)
    grad_p = R.calculate_press_grad(p1, p2, 0, 0, deltaL)[1]
    return D, Q, p1, p2, grad_p, tstac
end


function add_uncertainty(res::RcRes, port1, port2,
                         D_rel, D_abs,
                         Q_rel, Q_abs,
                         p_abs, p_std_fac, p_rel,
                         d_deltaL)
    @unpack Q, press, press_std, trun, irun  = res.resl
    p1, p2 = press[port1][irun], press[port2][irun]
    p1_std, p2_std = press_std[port1][irun], press[port2][irun]
    Q = Q[irun]
    D = R.get_time_series(res.resi, trun)[1]

    deltaL = R.port_dist(port1, port2, res) ± d_deltaL
    add_uncertainty(D, Q, p1, p2, deltaL,
                    D_rel, D_abs,
                    Q_rel, Q_abs,
                    p_abs, p1_std*p_std_fac, p2_std*p_std_fac,
                    p_rel)
end
function add_uncertainty(D, Q, p1, p2, deltaL,
                         D_rel, D_abs,
                         Q_rel, Q_abs,
                         p_abs, p1_std, p2_std, p_rel)
    D = relabs.(D, D_rel, D_abs)
    Q = relabs.(Q, Q_rel, Q_abs)
    p1 = relabs.(p1, p_rel, p_abs+p1_std)
    p2 = relabs.(p2, p_rel, p_abs+p2_std)
    grad_p = R.calculate_press_grad(p1, p2, 0, 0, deltaL)[1]
    return D, Q, p1, p2, grad_p
end

relabs(v, rel, abs) = v±(abs + v*rel)

function P.plot(t, v::Vector{<:Measurements.Measurement}; kwargs...)
    @show "here"
    va = Measurements.value.(v)
    err = Measurements.uncertainty.(v)
    P.fill_between(t, va-err, va+err, alpha=0.3)
    P.plot(t, va; kwargs...)
end
P.plot(v::Vector{<:Measurements.Measurement}; kwargs...) = P.plot(0:length(v)-1, v; kwargs...)

# https://matplotlib.org/examples/pylab_examples/ellipse_demo.html
function P.scatter(x::Vector{<:Measurements.Measurement}, y::Vector{<:Measurements.Measurement}, col::Vector;
                   scale_size=10, kwargs...)
    xx = Measurements.value.(x)
    xe = Measurements.uncertainty.(x)
    yy = Measurements.value.(y)
    ye = Measurements.uncertainty.(y)
    sizes = max.(xe,ye)*scale_size
    P.scatter(xx,yy,sizes,col; kwargs...)
end
function P.scatter(x::Vector{<:Measurements.Measurement}, y::Vector{<:Measurements.Measurement}; scale_size=10, kwargs...)
    xx = Measurements.value.(x)
    xe = Measurements.uncertainty.(x)
    yy = Measurements.value.(y)
    ye = Measurements.uncertainty.(y)
    sizes = max.(xe,ye)*scale_size
    P.scatter(xx,yy,sizes; kwargs...)
end


function P.errorbar(x::Vector{<:Measurements.Measurement}, y::Vector{<:Measurements.Measurement}; kwargs...)
    xx = Measurements.value.(x)
    xe = Measurements.uncertainty.(x)
    yy = Measurements.value.(y)
    ye = Measurements.uncertainty.(y)
    P.errorbar(xx,yy,ye,xe; kwargs...)
end
