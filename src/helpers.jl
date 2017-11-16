

export time_window2index_window, time2seconds, times_of_flow,
    plot_res, quantile_filter!, lowpass_filter2,
    lowpass_filter, median_filter


"Plot results"
function plot_res end

time2seconds(time,ex) = time2seconds(time, ex.experiment_start)
time2seconds(time,experiment_start::DateTime) =
    Dates.value.(Dates.Nanosecond.((time - Time(experiment_start))))/1e9

"""
    time_window2index_window(times, window)

Return the indices of all t∈times for which window[1]<=t<window[2].
"""
time_window2index_window(time, window) =
    findfirst(t->t>=window[1],time):findfirst(t->t>window[end], time)

"""
    times_of_flow(time, ex)
    times_of_flow(time, experiment_running)

Return times (in sec since start) and indices when water was flowing.
"""
times_of_flow(time, ex) = times_of_flow(time, ex.experiment_start, ex.experiment_running)
function times_of_flow(time, experiment_start::DateTime,
                       experiment_running::AbstractVector)
    irun = vcat([time_window2index_window(time, w) for w in experiment_running]...)
    t = time2seconds(time, experiment_start)
    dt = [0;diff(t)]
    trun = cumsum(dt[irun])
    return trun, irun
end

###############
# Filtering

"Image median filter"
median_filter(img, window=(3,3)) = Images.mapwindow(median!, img, window)

"""
Lowpass filter with a cutoff at `wavelength`,
expressed as ratio of wavelength/length(lines).
"""
function lowpass_filter(lines, wavelength::Float64=0.03)
    wavelength_pix = round(Int, wavelength*size(lines,1))
    frq = 2/wavelength_pix
    responsetype = DSP.Lowpass(frq)
    designmethod = DSP.Butterworth(4)
    round.(Int, DSP.filtfilt(DSP.digitalfilter(responsetype, designmethod), lines))
end

"""
Lowpass filter with a cutoff at `wavelength`,
expressed as ratio of wavelength/length(lines).
"""
function lowpass_filter2(v, wavelength_samples)
    frq = 2/wavelength_samples
    responsetype = DSP.Lowpass(frq)
    designmethod = DSP.Butterworth(4)
    DSP.filtfilt(DSP.digitalfilter(responsetype, designmethod), v)
end


"Filters outliers above certain quantiles"
function quantile_filter!(vec,
                          quant_cut=0.1, quant_offset=0, quant_factor=1,
                          pre_window=length(vec)÷10, post_window = pre_window+1,
                          niter = 10,
                          verbose=false)
    lenu = lenl = -1
    for i=1:niter
        upper, medu = find_spikes_quantile(vec, >,
                                           quant_cut, quant_offset, quant_factor,
                                           pre_window, post_window;
                                           verbose=verbose)
        lower, medl = find_spikes_quantile(vec, <,
                                           1-quant_cut, quant_offset, quant_factor,
                                           pre_window, post_window;
                                           verbose=verbose)

        lenu==length(upper) && break
        lenl==length(lower) && break
        lenu = length(upper)
        lenl = length(lower)
        vec[upper] = medu
        vec[lower] = medl
    end
    vec
end
