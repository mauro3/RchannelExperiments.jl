module RcLabView
using Parameters
using LsqFit, VAWTools
import PyPlot
const P=PyPlot

using ..RcSetup
import Base.Dates: Time, Millisecond, millisecond

export ExpLabView, ExpLabViewResults

const rho = 1000
const g = 9.81


@with_kw struct ExpLabView
    file::String # output file of Lab computer
    zero_press_window::Tuple{Time,Time}
    equal_press_windows::Vector{Tuple{Time,Time}}
end
@with_kw struct ExpLabViewResults @deftype Vector{Float64}
    time::Vector{Dates.Time} # Time of day
    t     # time since start of experiment
    trun  # total time the water was running since start of experiment
          # note length(t)!=length(trun)==length(t[irun])
    irun::Vector{Int} # time[irun] where times with flowing water
    @assert length(t)==length(time)==length(Q)
    @assert length(trun)==length(irun)==length(t[irun])
    Q # discharge m^3/s
    press::Dict{Port.PressPorts,Vector{Float64}} # pressure at different ports Pa
    press_std::Dict{Port.PressPorts,Vector{Float64}} # std of press (if several sensors)
end

function ExpLabViewResults(rcexp)
    @unpack explv, experiment_start, setup= rcexp

    time, Q, press0 = prep(explv)
    ch = RcSetup.used_channels(setup)
    press0 = press0[:,ch]
    t = Dates.value.(Dates.Nanosecond.((time - Time(experiment_start))))/1e9 # time in secs
    # set zero
    press0 = set_to_val(time, press0, explv.zero_press_window, 0)
    # correct
    press_out, means, vals, params, errs = correct_offsets(time, press0, explv.equal_press_windows)
    # assign columns of press_out to ports
    press = Dict{Port.PressPorts,Vector{Float64}}()
    press_std = Dict{Port.PressPorts,Vector{Float64}}()
    for (p,cs) in setup.port2channels
        length(cs)==0 && continue
        press[p] = vec(mean(press_out[:,cs],2))
        press_std[p] = vec(std(press_out[:,cs],2))
    end
    # figure out when water was running
    irun = find(Q.>0.0005)
    dt = [0;diff(t)]
    trun = cumsum(dt[irun])

    ExpLabViewResults(time, t, trun, irun, Q, press, press_std)
end

"""
    prep(explv::ExpLabView)

Returns
- time of day
- discharge
- pressures from all sensors
"""
function prep(explv::ExpLabView)
    # parse file-name to get time
    start = Time(DateTime(splitext(splitdir(explv.file)[2])[1], "yyyymmdd-HHMMSS"))
    out = readdlm(explv.file)
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

time_window2index_window(t, window) =
    findfirst(t, window[1]):findfirst(t, window[end])

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
- t -- time
- press -- pressure time series
- windows -- time windows when all sensors should read the same pressure

Return
- corrected pressure time series
- mean pressure value of each sensor at all windows
- mean pressure value of all sensors at all windows
- parameters (a,b) of line of best fit a*p+b
- errors (95%) for each sensor and each window

TODO: Maybe add ability to have the offset time dependent too.
"""
function correct_offsets(t, press, windows) # don't work:, weights=ones(length(windows)))
    nw = length(windows)
    ns = size(press,2)
    means = zeros(nw,ns)
    for (i,window) in enumerate(windows)
        window = time_window2index_window(t, window)
        means[i,:] =  mean(press[window,:],1)
    end
    vals = mean(means,2)[:]
    # assume a linear dependence
    offset(press, p) = p[1]*press + p[2]
    offsetinv(press, p) = (press -p[2])/p[1]
    function offsetJ(press,p)
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
        ydata = vals
        fit = curve_fit(offset, offsetJ, xdata, ydata, [1.0,0.0])
        @assert fit.converged
        params[:,j] = fit.param
        # We can estimate errors on the fit parameters,
        # to get 95% confidence error bars:
        errs[:,j] = estimate_errors(fit, 0.95)
        # transform press
        ofs = p -> offset(p, fit.param)
        press_out[:,j] = ofs.(press[:,j])
    end
    return press_out, means, vals, params, errs
end

tscale=60;
t0=[:exp,:logg][1];
figsize=(15,15);
units=[:mH2O, :Pa, :Bar][1]
function plot_res(reslv::ExpLabViewResults, ports=Port.all_ports; tscale=60, t0=[:exp,:logg,:trun][3], figsize=nothing, units=[:mH2O, :Pa, :Bar][1])
    if units==:mH2O
        trans = Pa_to_m_h2o
        yl = "p (m H2O)"
    elseif units==:Pa
        trans = identity
        yl = "p (Pa)"
    else
        error()
    end
    (tt,ii) = if t0==:logg
        (reslv.t-reslv.t[1],
         1:length(reslv.t))
    elseif t0==:exp
        (reslv.t,
         1:length(reslv.t))
    elseif t0==:trun
        (reslv.trun,
         reslv.irun)
    else
        error()
    end

    fig = figsize==nothing ? P.figure() : P.figure(figsize=figsize)
    map((k,v)->P.plot(tt/tscale, trans.(v[ii]), label=string(k)), keys(reslv.press), values(reslv.press))
    P.legend()
    P.xlabel("Time")
    P.ylabel(yl)
    fig
end

end
