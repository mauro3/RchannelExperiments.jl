module RcLabView
import ..RchannelExperiments: plot_res
using ..RchannelExperiments

using Parameters
using LsqFit, VAWTools
import PyPlot
const P=PyPlot
import DataStructures: OrderedDict
import DSP

using ..RcSetup
import ..RcImages: median_filter
import Base.Dates: Time, Millisecond, millisecond

export ExpLabView, ExpLabViewResults

@with_kw struct ExpLabView
    file::String # output file of Lab computer
    zero_press_window::Tuple{Time,Time}
    equal_press_windows::Vector{Tuple{Time,Time}}
end
@with_kw struct ExpLabViewResults @deftype Vector{Float64}
    exl::ExpLabView
    time::Vector{Dates.Time} # Time of day
    t     # time since start of experiment (s)
    trun  # total time the water was running since start of experiment
          # note length(t)!=length(trun)==length(t[irun]) (s)
    irun::Vector{Int} # time[irun] where times with flowing water
    @assert length(t)==length(time)==length(Q)
    @assert length(trun)==length(irun)==length(t[irun])
    Q # discharge m^3/s
    press::OrderedDict{Port.PressPorts,Vector{Float64}} # pressure at different ports Pa
    press_std::OrderedDict{Port.PressPorts,Vector{Float64}} # std of press (if several sensors)

    press_raw::OrderedDict{Port.PressPorts,Vector{Float64}} # pressure at different ports Pa
    press_std_raw::OrderedDict{Port.PressPorts,Vector{Float64}} # std of press (if several sensors)
end

function ExpLabViewResults(ex)
    @unpack exl, experiment_start, setup, experiment_running = ex

    time, Q, press0 = prep(exl.file)
    ch = RcSetup.used_channels(setup)
    press0 = press0[:,ch]
    # time in secs:
    t = time2seconds(time, experiment_start)
    # set zero
    press0 = set_to_val(time, press0, exl.zero_press_window, 0)
    # correct
    press_out, means, vals, params, errs = correct_offsets(
        time, press0, exl.equal_press_windows)
    # assign columns of press_out to ports
    press = OrderedDict{Port.PressPorts,Vector{Float64}}()
    press_std = OrderedDict{Port.PressPorts,Vector{Float64}}()
    press_raw = OrderedDict{Port.PressPorts,Vector{Float64}}()
    press_std_raw = OrderedDict{Port.PressPorts,Vector{Float64}}()
    for (p,cs) in setup.port2channels
        length(cs)==0 && continue
        press_raw[p] = vec(mean(press_out[:,cs],2))
        press_std_raw[p] = vec(std(press_out[:,cs],2))
        # median filter to remove spikes
        tmp, tmp_std = press_raw[p], press_std_raw[p]
        n = 5
        if n>0
            tmp = median_filter(tmp, (n,))
            tmp_std = median_filter(tmp_std, (n,))
        end

        nn = 5
        if nn>2
            tmp = lowpass_filter2(tmp, nn)
            tmp_std = lowpass_filter2(tmp_std, nn)
        end

        press[p] = tmp
        press_std[p] = tmp_std
    end
    # figure out when water was running
    #irun = find((Q.>0.5/1000) .& (t.>=0))
    trun, irun = times_of_flow(time, ex)

    ExpLabViewResults(exl, time, t, trun, irun, Q, press, press_std,
                      press_raw, press_std_raw)
end

"""
    prep(exl::ExpLabView)

Returns
- time of day
- discharge
- pressures from all sensors
"""
function prep(file)
    # parse file-name to get time
    start = file2time(file)
    out = readdlm(file)
    time = start + Millisecond.(1000*out[:,1])
    Q = out[:,2]/1000 # m^3/s
    # the rest depends a bit on the setup
    press = out[:,3:end]
    press .= mBar_to_Pa.(press)
    return time, Q, press
end

"Transform pressure from mBar to Pascal"
mBar_to_Pa(p) = p*100

"Transform pressure from mBar to m H2O"
function mBar_to_m_h2o(p)
    # mbar -> Pa -> m
    p*100 /rho/g
end
Pa_to_m_h2o(p) = p/rho/g

"Transforms the labview file into a DateTime of the logging start"
file2datetime(file) = DateTime(splitext(splitdir(file)[2])[1], "yyyymmdd-HHMMSS")
file2time(file) = Time(file2datetime(file))


"""
    set_to_val(t, press, window, val=nothing)

Corrects the offsets in the pressure transducers.  Presumably to set
the zero.

- t -- time
- press -- pressure time series
- window -- time window which to use as reference of the 0,
            i.e. subtract the mean pressure during this time window
- val -- value to which to set the mean to.  Defaults to the mean of the mean.

"""
function set_to_val(t, press, window, val=nothing)
    window = time_window2index_window(t, window)
    means = mean(press[window,:],1)
    if val==nothing
        val = mean(means)
    end
    press .- means + val
end

"""
    correct_offsets(t, press, windows)

Correct the pressure with an offset which is dependent on the absolute
pressure.  This is needed for the pressure sensors.

Input
- t -- times as Vector{Time}
- press -- pressure time series
- windows -- time windows when all sensors should read the same pressure

For each window a mean pressure is calculated for each sensor `pm` and also the mean
of all means `mm`.  A line is then fitted between `pm` and `mm` and used to correct
each sensor.  (This means that the weight of each window is equal, i.e. independent of
its width.)

Return
- corrected pressure time series
- mean pressure value of each sensor at all windows
- mean pressure value of all sensors at all windows
- parameters (a,b) of line of best fit a*p+b
- errors (95%) for each sensor and each window

TODO: Maybe add ability to have the offset time dependent too.
"""
function correct_offsets(time, press, windows; verbose=false)
    nw = length(windows)
    ns = size(press,2)
    means = zeros(nw,ns)
    for (i,window) in enumerate(windows)
        window = time_window2index_window(time, window)
        @assert length(window)>0 "Window has zero length"
        means[i,:] =  mean(press[window,:],1)
    end
    @assert !any(isnan.(means)) "Found NaNs in mean"
    meanmeans = mean(means,2)[:]
    # assume a linear dependence
    offset(press, p) = p[1]*press + p[2]
    offsetinv(press, p) = (press -p[2])/p[1]
    function offsetJ(press,p)
        # Jacobian of offset
        J = Array{Float64}(length(press),length(p))
        J[:,1] = press    #doffset/dp[1]
        J[:,2] = ones(press)  #doffset/dp[2]
        J
    end

    # fit a and b such that error is minimized
    params = zeros(2,ns)
    errs = zeros(2,ns)
    press_out = similar(press)
    for j=1:ns
        xdata = means[:,j]
        ydata = meanmeans
        fit = curve_fit(offset, offsetJ, xdata, ydata, [1.0,0.0])
        if verbose
            @show j
            @show fit.converged
            @show fit.param
            @show estimate_errors(fit, 0.95)
            @show xdata
            @show ydata
            P.clf()
            P.scatter(xdata, ydata)
            P.plot(sort(xdata), offset(sort(xdata), fit.param))
            P.plot(sort(xdata), offset(sort(xdata), [1.0,0.0]))
            sleep(5)
        end
        @assert fit.converged "Could not find a line of best fit to correct pressure offsets"
        params[:,j] = fit.param
        # We can estimate errors on the fit parameters,
        # to get 95% confidence error bars:
        errs[:,j] = estimate_errors(fit, 0.95)
        # transform press
        ofs = p -> offset(p, fit.param)
        press_out[:,j] = ofs.(press[:,j])
    end
    return press_out, means, meanmeans, params, errs
end

#################
# Plotting

function plot_file(file, inds=1:10)
    time, Q, press0 = prep(file)
    start = file2time(file)
    t = Dates.value.(Dates.Millisecond.(time-start))/1000/60
    P.figure()
    ax = P.subplot(2,1,1)
    P.plot(t, Q*1000)
    ax = P.subplot(2,1,2,sharex=ax)
    P.plot(t, mBar_to_m_h2o.(press0[:,inds]))
end

function plot_res(ex, resl::ExpLabViewResults; ports=RcSetup.used_ports(ex.setup),
                  tscale=60, t0=[:exp,:logg,:trun, :time][3],
                  fig_axs=P.subplots(3,1,sharex="all"),
                  units=[:mH2O, :Pa, :Bar][1])
    fig,axs = fig_axs

    # value transform and y-label
    trans, yl1, yl2 = if units==:mH2O
        (Pa_to_m_h2o, "p (m H2O)", "Δp ()")
    elseif units==:Pa
        (identity, "p (Pa)", "Δp (Pa/m)")
    else
        error()
    end

    # time and indices
    tt,ii = if t0==:logg
        (resl.t-resl.t[1],
         1:length(resl.t))
    elseif t0==:exp
        (resl.t,
         1:length(resl.t))
    elseif t0==:trun
        (resl.trun,
         resl.irun)
    elseif t0==:time
        (Dates.value.(resl.time)/1000/60/60/24, # Convert time from milliseconds from day 0 to days from day 0
         1:length(resl.t))
    else
        error()
    end
    @show ports[1:end-1],ports[2:end]
    grad, grad_std, grad_ports = calculate_press_grad_all(ex, resl, ports[1:end-1], ports[2:end]);

    P.axes(axs[1])
    P.grid("on")
    P.plot(tt/tscale, resl.Q[ii]*1000)
    P.ylabel("Q (l/s)")
    P.title(ex.name)

    P.axes(axs[2])
    P.grid("on")
    for (k,v,vs) in zip(ports, [resl.press[p] for p in ports], [resl.press_std[p] for p in ports])
        t0==:time || P.fill_between(tt/tscale, trans.(v[ii]-vs[ii]), trans.(v[ii]+vs[ii]), alpha=0.3)
        t0==:time || P.plot(tt/tscale, trans.(v[ii]), label=string(k))
        t0==:time && error() # does not work P.plot_date(tt, trans.(v[ii]), label=string(k))
    end

    P.legend()
    P.ylabel(yl1)

    P.axes(axs[3])
    P.grid("on")

    for (i,p) in enumerate(grad_ports)
        P.fill_between(tt/tscale, trans.(grad[ii,i]-grad_std[ii,i]), trans.(grad[ii,i]+grad_std[ii,i]), alpha=0.3)
        P.plot(tt/tscale, trans.(grad[ii,i]), label=string(p))
    end

    P.ylabel(yl2)
    P.xlabel("Time")

    # P.xlim((93,97))
    # P.ylim((1.55,2.5))
end

end
