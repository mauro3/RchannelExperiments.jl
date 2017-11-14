"""
    elapsed_time(res, i)

Elapsed time in seconds for i-th picture
"""
elapsed_time(res::ExpImgsResults, i)::Float64 =
    (Dates.Second(res.capture_times[i] - res.ep.experiment_start)).value

"""
    get_time_series(res::ExpImgsResults)

Get time series data as matrix:
- time
- mean diameter
- median diameter
- 10% quantile
- 90% quantile
- scalloping
- symmetry
"""
function get_time_series(res::ExpImgsResults)
    @unpack dia_mean, dia_quant, scalloping, tb_cor = res
    t = elapsed_time.(res, 1:length(res))
    return t, dia_mean, dia_quant', scalloping, tb_cor
end

"Image median filter"
median_filter(img, window=(3,3)) = mapwindow(median!, img, window)

#################
# findpeaks
# put into VAWTools

"""
    findpeaks(vec; min_dominance=0.0, min_prominence=0.0)

Return

- peaks

- index of peaks

- prominence (difference between peak height and the highest low-point to
  traverse to get to a higher point)

- dominance (separation to next higher point, or maximal distance to the edge)

TODO: use troughs vector
"""
function findpeaks(vec; min_dominance=0.0, min_prominence=0.0)
    ## find all peaks and all troughs
    peaks = eltype(vec)[]
    ips = Int[]
    troughs = eltype(vec)[]
    its = Int[]
    for i=2:length(vec)-1
        if vec[i-1]<vec[i] && vec[i+1]<=vec[i]
            push!(peaks, vec[i])
            push!(ips, i)
        elseif vec[i-1]>vec[i] && vec[i+1]>=vec[i]
            push!(troughs, vec[i])
            push!(its, i)
        end
    end

    # calculate dominance and prominence
    doms = Int[]
    proms = eltype(vec)[]
    for (i,ip) in enumerate(ips)
        ele = peaks[i]
        left = @view vec[ip-1:-1:1]
        right = @view vec[ip+1:end]
        dom1 = findfirst(x->x>ele, left)
        dom2 = findfirst(x->x>ele, right)
        if dom1==0 && dom2==0
            # highest peak
            max_dist_to_margin = max(ip, length(vec)-ip)
            push!(doms, max_dist_to_margin)
            # push!(doms, length(vec))
            push!(proms, ele-minimum(vec))
        else
            if dom1==0 # no higher peak on left
                push!(doms,dom2)
                #min1 = minimum(@view vec[ip-1:-1:1])
                min2 = minimum(@view vec[ip+1:ip+dom2])
                min1 = min2
            elseif dom2==0 # no higher peak on right
                push!(doms,dom1)
                min1 = minimum(@view vec[ip-1:-1:ip-dom1])
                #min2 = minimum(@view vec[ip+1:end])
                min2 = min1
            else
                push!(doms,min(dom1,dom2))
                min1 = minimum(@view vec[ip-1:-1:ip-dom1])
                min2 = minimum(@view vec[ip+1:ip+dom2])
            end
            push!(proms, ele - max(min1,min2))
        end
    end

    inds = (doms.>=min_dominance) .& (proms.>=min_prominence)
    return peaks[inds], ips[inds], proms[inds], doms[inds]
end

"""
    findpeaks_onesided(vec, dir::Symbol; min_dominance=0.0, min_prominence=0.0)

As findpeaks but only scans prominence and dominance in direction `dir` (:left,:right)

Return

- peaks

- index of peaks

- prominence (difference between peak height and lowest point to
  traverse to get to a higher point)

- dominance (separation to next higher point, or minimal distance to the edge)

The API inspired by
https://ch.mathworks.com/help/signal/ref/findpeaks.html.
"""
function findpeaks_onesided(vec, dir; min_dominance=0.0, min_prominence=0.0)
    @assert dir in [:left,:right]
    ## find all peaks and all troughs
    peaks = eltype(vec)[]
    ips = Int[]
    troughs = eltype(vec)[]
    its = Int[]
    for i=2:length(vec)-1
        if vec[i-1]<vec[i] && vec[i+1]<=vec[i]
            push!(peaks, vec[i])
            push!(ips, i)
        elseif vec[i-1]>vec[i] && vec[i+1]>=vec[i]
            push!(troughs, vec[i])
            push!(its, i)
        end
    end

    # calculate dominance and prominence
    doms = Int[]
    proms = eltype(vec)[]
    for (i,ip) in enumerate(ips)
        ele = peaks[i]
        side = dir==:left ? @view(vec[ip-1:-1:1]) : @view(vec[ip+1:end])
        dom = findfirst(x->x>ele, side)
        if dom==0
            # highest peak to border
            push!(doms, length(side))
            push!(proms, ele-minimum(side))
        else
            push!(doms,dom)
            side2 = dir==:left ? @view(vec[ip-1:-1:ip-dom]) : @view(vec[ip+1:ip+dom])
            push!(proms, ele - minimum(side2))
        end
    end

    inds = (doms.>=min_dominance) .& (proms.>=min_prominence)
    return peaks[inds], ips[inds], proms[inds], doms[inds]
end



find_first_big_peak(vec, dir=:left) = find_first_big_peak(findpeaks_onesided(vec, dir==:left ? :right : :left)[1:3]..., dir)
function find_first_big_peak(peaks, locs, proms, dir=:left)
    # run with a median filter over it
end

## This does not work well:
# using OnlineStats
# function find_outlier_median(vec, quant=0.1, exp_weight=1/4,
#                              what=[:below,:above,:both][3];
#                              verbose=false)
#     @assert 0<exp_weight<=1
#     # mean for the first 1/0.2 observations, then Exponential.
#     # e-fold of importance given by exp_weight

#     w = ExponentialWeight(exp_weight)

#     # none of them seem quite kosher:
#     Qu = [QuantileMM, QuantileSGD, QuantileMSPI][3]

#     s = Series(w, Qu([quant, 0.5, 1-quant]))

#     outliers = Int[]
#     q1 = Float64[]
#     q2 = Float64[]
#     q3 = Float64[]

#     fit!(s, vec[1], 1.0) # use 1.0 weight, see https://github.com/joshday/OnlineStats.jl/issues/93
#     for i = 2:length(vec)
#         v = vec[i]
#         sv = value(s)[1]
#         push!(q1,sv[1])
#         push!(q2,sv[2])
#         push!(q3,sv[3])
#         if what==:both
#             sv[1]<v<sv[3] || push!(outliers, i)
#         elseif what==:below
#             sv[1]<v || push!(outliers, i)
#         elseif what==:above
#             v<sv[3] || push!(outliers, i)
#         end
#         fit!(s, v)
#     end
#     sv = value(s)[1]
#     push!(q1,sv[1])
#     push!(q2,sv[2])
#     push!(q3,sv[3])

#     if verbose
#         P.figure()
#         P.plot(1:length(vec),vec)
#         map(x->P.plot(2:length(vec)+1,x) ,(q1,q2,q3))
#         P.plot(outliers, vec[outliers], ".")
#     end

#     w = Bounded(EqualWeight(), exp_weight)
#         w = ExponentialWeight(exp_weight)
#     @show Series(vec, w, QuantileMSPI([quant, 0.5, 1-quant]))
#     outliers
# end


"""
    find_spikes_quantile(vec, comp_vec_q=>, quant=0.9,
                          quant_offset=mean(vec)/2, quant_factor=1,
                          pre_window=5, post_window = 2;
                          verbose=false)

Find index of spikes/outliers by comparing each value with a certain quantile
of a running window.

    vec -- input vector
    comp_vec_q=> -- comparison to use: > means remove outliers above
    quant=0.9 -- what quantile to compare against
    quant_offset=mean(vec)/2
    quant_factor=1
    pre_window=5, post_window=0 -- window size over which to look at values

The comparison is then `comp_vec_q(v[i], q*quant_factor + quant_offset)`, i.e.
the quantile value `q` is multiplied by the `quant_factor` and `quant_offset` is added.

The window for the i-th value is `i-pre_window:i+post_window-1`, i.e. if
`post_window==0` then the value itself is not considered.

Keyword:
    verbose=false -- if true plot

Return

- index of spikes/outliers
- value of median over window at outliers

"""
function find_spikes_quantile(vec, comp_vec_q= >,
                               quant=0.9, quant_offset=median(vec), quant_factor=1,
                               pre_window=length(vec)÷10, post_window = 0;
                              verbose=false)
    ex = 3
    pre_w = 1./(pre_window+1:-1:2)
    post_w = 1./(1:post_window).^ex
    w = StatsBase.Weights([pre_w; post_w])
    outliers = Int[]
    qs = Float64[]
    med = Float64[]
    for i = pre_window+1:length(vec)-post_window
        vi = vec[i]
        v = @view vec[i-pre_window:i+post_window-1]
        if length(v)>0
            q = quant_factor * quantile(v, w, quant) + quant_offset
            push!(qs, q)
            push!(med, median(v))
            if comp_vec_q(vi, q)
                push!(outliers, i)
            end
        end
    end
    if verbose
        P.plot(1:length(vec),vec)
        map(x->P.plot(pre_window+1:length(vec)-post_window,x) ,(qs,))
        P.plot(outliers, vec[outliers], ".")
    end
    outliers, med[outliers-pre_window]
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


#############
# Saving and Loading
function filename(ep::ExpImgs)
    @unpack dir, thin_num, ns, algo = ep
    "$(dir)-thin$(thin_num)-ns$(length(ns))-$(algo)"
end
filename(res::ExpImgsResults) = filename(res.ep)
"""
    save_result(res::ExpImgsResults; overwrite=false, store_imgs=false)


"""
function save_result(res::ExpImgsResults; overwrite=false, store_imgs=false)
    fln = "output/"*filename(res)*".jld"
    if isfile(fln) && !overwrite
        warn("File exists and overwrite==false, not saving!")
        return nothing
    end
    if !store_imgs && length(res.imgs)>0
        res = ExpImgsResults(res,
                             imgs=[])
    end
    JLD.jldopen(fln, "w") do file
        file["res"] = res
    end
    nothing
end
load_result(fln) = JLD.load(fln)
load_result(ep::ExpImgs) = load_result(filename(ep)*".jld")


###############
# Filtering

"""
Lowpass filter with a cutoff at `wavelength`,
expressed as ratio of wavelength/image width.
"""
function lowpass_filter(lines, wavelength=0.03)
    wavelength_pix = wavelength*size(lines,1)
    @show frq = 2/wavelength_pix
    responsetype = DSP.Lowpass(frq)
    designmethod = DSP.Butterworth(4)
    round.(Int, DSP.filtfilt(DSP.digitalfilter(responsetype, designmethod), lines))
end
