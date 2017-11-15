#########
# Width
getmid(img) = fld1(size(img,1),2)
getmid(ep::ExpImgs) = fld1(size(ep)[1],2)

top2ind(top, img_or_ep) = getmid(img_or_ep) - top + 1
bottom2ind(bottom, img_or_ep) = bottom + getmid(img_or_ep)
top_inds(img_or_ep) = getmid(img_or_ep):-1:1
bottom_inds(img_or_ep) = getmid(img_or_ep)+1:size(img,1)

################
# width via color tresholding


function findtrough(img, title; verbose=false)
    nbins = 100
    hi = fit(Histogram, img[:], closed=:left, nbins=nbins)
    peaks, ips, proms, doms = findpeaks(-hi.weights[1:(end-end÷4)])
    i = findmax(proms[ips.<(nbins÷2)])[2] # look only if first half of histogram
    if proms[i]<length(img)/500
        warn("No color boundary found in $title")
        return 0
    end
    verbose && P.plot(1:length(hi.weights), -hi.weights)
    #verbose && P.bar(hi.edges[1][1:end-1], hi.weights)
    verbose && @show ips[i]
    ips[i]
end


"""

"""
function keep_only_big_regions(img, loc, siz=round(Int,(sqrt(length(img))/20)^2),
                               k=10; verbose=false)
    segments = felzenszwalb(img, k, siz)
    verbose && @show segments

    # map(i->segment_mean(segments,i), labels_map(segments))
    tmp = labels_map(segments)
    # set region at `loc` to true, rest to false.
    tmp.==tmp[loc...]
end


function channel_width_thresholding(img_, loc, median_filter_region, title; verbose=false)
    trough = findtrough(img_, title, verbose=verbose)
    trough==0 && return zeros(Int,size(img_,2)),zeros(Int,size(img_,2)),zeros(Int,size(img_,2))
    img_ = img_.>trough
    if median_filter_region!=(0,0)
        img_ = median_filter(img_, median_filter_region)
    end
    img_3 = keep_only_big_regions(img_, loc,
                                  round(Int,(sqrt(length(img_))/50)^2),
                                  verbose=verbose);
    verbose && imshow(img_);
    verbose && imshow(img_3);

    mid = getmid(img_)
    # top = squeeze(sum(img_3[1:mid,:],1),1)
    # bottom = squeeze(sum(img_3[mid+1:end,:],1),1)
    top = [findlast(img_3[mid:-1:1,i]) for i=1:size(img_,2)]
    bottom = [findlast(img_3[mid+1:end,i]) for i=1:size(img_,2)]

    # remove outliers
    quantile_filter!(top, 0.05, 5, 1, 5, 6)
    quantile_filter!(bottom, 0.05, 5, 1, 5, 6)

    return top, bottom, top+bottom
end


###############
# Edge detection
"""
    channel_width_edgedetection(img_, minhalfwidth, gauss_w = 10, quant=0.8, gap=11; verbose=false)

Determine width using edge detection.
"""
channel_width_edgedetection(img_, ep; verbose=false) =
    channel_width_edgedetection(img_, ep.minhalfwidth, ep.gauss_w, ep.quant, ep.gap, verbose=verbose)
function channel_width_edgedetection(img_, minhalfwidth, gauss_w,
                                     quant, gap, morphgrad=false; verbose=false)
    @assert !iseven(gap) "gap must be a odd number"
    # works on img_
    verbose && imshow(img_);

    mid = getmid(img_)
    ncols = size(img_,2)

    # Gauss filter:
    img_f = imfilter(img_, Kernel.gaussian(gauss_w));

    # eges:
    if morphgrad
        mag = morphogradient(img_f)
    else
        grad_x, grad_y, mag, orient = imedge(img_f);
        # Only look at gradient in x-direction
        # and only when going from more green to less green
        # in outward direction.
        grad_x[mid:end,:] = -grad_x[mid:end,:]
        grad_x = -grad_x
        grad_x[grad_x.<0] = 0
        mag = grad_x
    end
    verbose && imshow(mag);
    #verbose && imshow(grad_x);


    # # Find significant peaks on each column:
    # out = falses(mag);
    # out_ = zeros(mag);
    # for i=1:ncols
    #     peaks,locs,proms,doms = findpeaks(mag[:,i]);
    #     big = find(x->x>quantile(peaks, quant) && x>quantile(proms, quant), peaks)
    #     out[locs[big],i] = true
    #     out_[locs,i] = proms
    # end
    # verbose && imshow(out);

    # Find significant peaks on each column, version 2:
    out = falses(img_);
    vverbose = false;
    for i=1:ncols
        # right (lower) border
        # peaks,locs,proms,doms = findpeaks(mag[:,i]);
        peaks_t,locs_t,proms_t,doms_t = findpeaks(mag[mid:-1:1,i]);
        peaks_b,locs_b,proms_b,doms_b = findpeaks(mag[mid+1:end,i]);
        locs_tt = mid +1 - locs_t
        locs_bb = locs_b-1+mid

        big_t = find_spikes_quantile(proms_t, >,
                                     0.9, mean(proms_t), 1,
                                     length(proms_t)÷10, 0,
                                     verbose=vverbose)[1]
        big_b = find_spikes_quantile(proms_b, >,
                                     0.9, mean(proms_b), 1,
                                     length(proms_b)÷10, 0,
                                     verbose=vverbose)[1]

        out[locs_tt[big_t],i] = true
        out[locs_bb[big_b],i] = true
    end
    verbose && imshow(out);



    # # find most prominent peak on each half:
    # out_p = falses(mag);
    # for i=1:ncols
    #     # right (lower) border
    #     peaks,locs,proms,doms = findpeaks(mag[:,i]);
    #     big = find(x->x>quantile(peaks, quant) && x>quantile(proms, quant), peaks)
    #     out[locs[big],i] = true
    #     out_[locs,i] = proms
    # end
    # verbose && imshow(out);

    # # Just find the two biggest peaks in both halves
    # top = zeros(Int, ncols);
    # bottom = zeros(Int, ncols);
    # for i=1:ncols
    #     top[i] = findmax(mag[1:mid-minhalfwidth,i])[2]
    #     bottom[i] = findmax(mag[mid+minhalfwidth+1:end,i])[2] + minhalfwidth
    # end
    # return top, bottom, top+bottom

    ## remove small bits
    # labs0 = label_components(out, trues(3,3));
    # inds0 = component_indices(labs0);
    # len0 = component_lengths(labs0);
    # # drop single & double pixels
    # for (i,l) in enumerate(len0)
    #     if l<3
    #         labs0[inds0[i]] = 0
    #     end
    # end
    # out = labs0.>0;
    # verbose && imshow(out);

    ## join bigger ones and remove left-over smaller ones
    labs = label_components(out, trues(gap,gap));
    inds = component_indices(labs);
    len = component_lengths(labs);
    # drop all components with are below certain length
    tlen = 50
    for (i,l) in enumerate(len)
        if l<tlen
            labs[inds[i]] = 0
        end
    end
    out = labs.>0;
    verbose && imshow(labs.>0);


    ## TODO: actually join the big pieces

    # now determine distance from middle+minhalfwidth to first peak
    top = zeros(Int, size(mag,2))
    bottom = zeros(Int, size(mag,2))
    for i=1:size(mag,2)
        tmp = findfirst(@view out[mid-minhalfwidth:-1:1,i])
        top[i] = tmp==0 ? -1 : tmp + minhalfwidth
        tmp = findfirst(@view out[mid+minhalfwidth+1:end, i])
        bottom[i] = tmp==0 ? -1 : tmp + minhalfwidth
    end
    verbose && PyPlot.figure()
    verbose && PyPlot.plot(top, ":")
    verbose && PyPlot.plot(bottom, ":")

    # remove outliers
    quantile_filter!(top,
                     0.05, 5, 1, 5, 6)
    quantile_filter!(bottom,
                     0.05, 5, 1, 5, 6)

    # top[top.==-1] = NaN
    # bottom[bottom.==-1] = NaN

    verbose && PyPlot.plot(top)
    verbose && PyPlot.plot(bottom)
    return top, bottom, top+bottom
end

"""
    filter1d_outliers!(vec::Vector{<:Integer}, halfwindow=1)

Filter outliers using median.  Probably better ways to do this.
"""
function filter1d_outliers!(vec::Vector{<:Integer}, halfwindow=1)
    len = length(vec)
    for i=eachindex(vec)
        vec[i] = round.(eltype(vec), median(@view vec[max(1,i-halfwindow):min(len,i+halfwindow)]))
    end
    vec
end

# function plottop()
#     for i=1:size(mag,2)
#         PyPlot.plot(top[i]*0+i, top[i],".")
#         PyPlot.plot(bottom[i]*0+i, bottom[i],".")
#     end
# end

###############
#

"""
    channel_width(ep::ExpImgs; verbose=false)

Returns the best of both outlines for all images in a ep.  Returns

- tops, bottoms: pixel distance to top and bottom outline in terms of
  pixels of the cropped image (i.e. not corrected with center_dist)
- t_dist, b_dist: distance in meters to top and bottom from a line
  through the original image center and parallel to (ep.p1, ep.p2).

"""

function channel_width(ep::ExpImgs; verbose=false, vverbose=false, saveit=true,
                       overwrite=false, store_imgs=false)
    @unpack dir, ns, p1, p2, thin_num, algo,
            minhalfwidth_orig, gauss_w, quant, gap, median_filter_region = ep

    println("\n Figuring out R-channel width for pictures in dir: $dir\n")
    minhalfwidth = minhalfwidth_orig÷thin_num

    tops = Int[]
    bottoms = Int[]
    last_top = [minhalfwidth]
    last_bottom = [minhalfwidth]
    ncol = size(ep)[2]
    center_dists = Int[]
    thumbs = []
    imgs = []
    capture_times = DateTime[]
    @showprogress for (n,fl) in zip(ns, image_files(ep))
        img, center_dist, capture_time = prep_img(fl, ep; verbose=vverbose)
        t, b = _channel_width(img, last_top, last_bottom, ep, "$n:  $fl",
                             algo,
                             verbose=verbose,
                             vverbose=vverbose)
        last_top = t
        last_bottom = b
        append!(tops, t)
        append!(bottoms,b)
        push!(center_dists, center_dist)
        push!(thumbs, thin(img, 16))
        store_imgs && push!(imgs, img)
        push!(capture_times, capture_time)
    end
    tops = reshape(tops, size(ep)[2], length(tops)÷size(ep)[2])
    bottoms = reshape(bottoms, size(ep)[2], length(tops)÷size(ep)[2])
    ts_dist = pixel2meter.(tops.+center_dists', ep)
    bs_dist = pixel2meter.(bottoms.-center_dists', ep)
    res = ExpImgsResults(ep,
                         capture_times,
                         center_dists,
                         tops,
                         bottoms,
                         ts_dist,
                         bs_dist,
                         rchannel_stats(ts_dist,bs_dist)...,
                         thumbs,
                         imgs,
                         Dict()
                         )
    saveit && save_result(res, overwrite=overwrite, store_imgs=store_imgs)
    return res
end


function _channel_width(img, last_top, last_bottom, ep::ExpImgs, title, algo;
                       verbose=false, vverbose=false)
    @unpack minhalfwidth_orig, gauss_w, quant, gap, median_filter_region, thin_num = ep
    minhalfwidth = minhalfwidth_orig÷thin_num

    if length(last_bottom)==1
        last_bottom = zeros(Int, size(img,2)) .+ last_bottom
    end
    if length(last_top)==1
        last_top = zeros(Int, size(img,2)) .+ last_top
    end

    # minhalfwidth_cur = max(min(minimum(last_top), minimum(last_bottom)) - minhalfwidth÷5, minhalfwidth)
    minhalfwidth_cur = minhalfwidth
    if algo in [:both, :thresh]
        top_t, bottom_t, tot_t = channel_width_thresholding(img,
                                                            (size(img,1)÷2, size(img,2)÷2),
                                                            median_filter_region, title;
                                                            verbose=vverbose)
        if algo==:thresh
            new_top, new_bottom = top_t, bottom_t
        end
    end
    if algo in [:both, :edge]
        top_e, bottom_e, tot_e = channel_width_edgedetection(img, minhalfwidth_cur,
                                                             gauss_w, quant, gap;
                                                             verbose=vverbose)
        if algo==:edge
            new_top, new_bottom = top_e, bottom_e
        end
    end
    if !(algo in [:both, :edge, :thresh])
        error("algo=$algo not recognized")
    end

    if algo==:both
        # now produce a best vector somehow:

        # Go through to find points of agreement which are farther than v_old.
        # Of those pick the closest to v_old, this will serve as starting point
        tv, bv = -1, -1
        ta, ba = 10^10, 10^10
        ti, bi = 0, 0
        for i=1:length(top_e)
            tmp_,vv = how_close(top_t[i], top_e[i], last_top[i])
            if tmp_<ta
                ta = tmp_
                ti = i
                tv = vv
            end
            tmp_,vv = how_close(bottom_t[i], bottom_e[i], last_bottom[i])
            if tmp_<ba
                ba = tmp_
                bi = i
                bv = vv
            end
        end
        if tv<1
            # no starting point found!
            error("no starting point found for tops in: $title")
        end
        if bv<1
            error("no starting point found for bottoms in: $title")
        end

        new_top = zeros(last_top)
        new_top[ti] = tv
        if ti<length(new_top)
            new_top[ti+1] = pick_next_point(top_t[ti+1], top_e[ti+1], new_top[ti], last_top[ti])
            for i=ti+2:length(new_top)
                new_top[i] = pick_next_point(top_t[i], top_e[i], new_top[i-1], new_top[i-2], last_top[i])
            end
        end
        if ti>1
            new_top[ti-1] = pick_next_point(top_t[ti-1], top_e[ti-1], new_top[ti], last_top[ti])
            for i=ti-2:-1:1
                new_top[i] = pick_next_point(top_t[i], top_e[i], new_top[i+1], new_top[i+2], last_top[i])
            end
        end

        new_bottom = zeros(last_bottom)
        new_bottom[bi] = bv
        if bi<length(new_bottom)
            # start
            new_bottom[bi+1] = pick_next_point(bottom_t[bi+1], bottom_e[bi+1], new_bottom[bi], last_bottom[bi])
            for i=bi+2:length(new_bottom)
                new_bottom[i] = pick_next_point(bottom_t[i], bottom_e[i], new_bottom[i-1], new_bottom[i-2], last_bottom[i])
            end
        end
        if bi>1
            # start
            new_bottom[bi-1] = pick_next_point(bottom_t[bi-1], bottom_e[bi-1], new_bottom[bi], last_bottom[bi])
            for i=bi-2:-1:1
                new_bottom[i] = pick_next_point(bottom_t[i], bottom_e[i], new_bottom[i+1], new_bottom[i+2], last_bottom[i])
            end
        end
    end

    window = 10
    niter = 30
    quantile_filter!(new_top, 0.1, 2, 1, window, window+1, niter)
    quantile_filter!(new_bottom, 0.1, 2, 1, window, window+1, niter)

    window = 11
    new_top = round.(Int,mapwindow(mean, new_top, window))
    new_bottom = round.(Int,mapwindow(mean, new_bottom, window))

    if verbose
        if algo==:both
            plot_all_n_new(img, new_top, new_bottom, top_t, bottom_t, top_e, bottom_e;
                           col="r", label="", title=title,
                           legend=false)
        elseif algo==:edge
            plot_all_n_new(img, new_top, new_bottom, top_e, bottom_e, top_e, bottom_e;
                           col="r", label="", title=title,
                           legend=false)
        elseif algo==:thresh
            plot_all_n_new(img, new_top, new_bottom, top_t, bottom_t, top_t, bottom_t;
                           col="r", label="", title=title,
                           legend=false)
        else
            error()
        end
        P.draw()
    end
    return new_top, new_bottom
end

"""
    how_close(v1, v2, v_old, weightdiv=5)

Squared distance of three points, the distance to the last has only
1/weightdiv the weight.  Returns also the v1 or v1 which is closer to
v_old.
"""
function how_close(v1, v2, v_old, weightdiv=5)
    if v_old<0
        return 10^6, -1
    else
        fit = ((v1-v_old)^2 + (v2-v_old)^2)÷weightdiv + (v2-v1)^2
        tmp = (v1-v_old)^2>(v2-v_old)^2 ? v2 : v1
        return fit, tmp
        # if v1==v2
        #     # weight closeness to v_old
        #     return sqrt((v1-v_old)^2)/5000, v1
        # else
        #     tmp = (v1-v_old)^2>(v2-v_old)^2 ? v2 : v1
        #     return (v1-v2)^2.0, tmp
        # end
    end
end

"""
    pick_next_point(v1, v2, vp, v_old)
    pick_next_point(v1, v2, vp, vpp, v_old)

Try to pick the better of v1 and v2.
"""
function pick_next_point(v1, v2, vp, v_old)
    fits = zeros(2)
    for (i,v) in enumerate((v1,v2))
        # oldness = v>v_old ? (v-v_old)^2 : 2*(v-v_old)^2
        oldness = (v-v_old)^2
        fits[i] = (v-vp)^2 + 2*oldness÷3
    end
    ii = findmin(fits)[2]
    if ii==1
        return v1
    elseif ii==2
        return v2
    else # mean is best -> try harder to find which v is better
        error()
    end
    ## old:
    # if v1>=v_old && v2>=v_old
    #     # both good.  Pick the one closer to
    #     # previous point:
    #     return (v1-vp)^2>(v2-vp)^2 ? v2 : v1
    # elseif v1<v_old && v2<v_old
    #     # oh no!
    #     return v_old
    # elseif v1<v_old
    #     return v2
    # elseif v2<v_old
    #     return v1
    # else
    #     error()
    # end
end

function pick_next_point(v1, v2, vp, vpp, v_old)
    fits = zeros(2)
    dv = vpp-vp
    for (i,v) in enumerate((v1,v2))
        # oldness = v>v_old ? (v-v_old)^2 : 2*(v-v_old)^2
        oldness = (v-v_old)^2
        dvv = v-vp
        fits[i] = (dv-dvv)^2 + oldness÷5
    end
    ii = findmin(fits)[2]
    if ii==1
        return v1
    elseif ii==2
        return v2
    else
        error()
    end
end


function rchannel_stats(tdist, bdist)
    d = tdist+bdist
    dia_mean = mean(d)
    dia_std = std(d)
    dia_quant = quantile(d, [0.1, 0.5, 0.9])
    scalopping = (dia_quant[3]-dia_quant[1])/dia_quant[2]
    scalopping2 = dia_std./dia_mean
    tb_cor = cor(tdist, bdist)
    return dia_mean, dia_quant, scalopping, tb_cor
end
function rchannel_stats(tdists::Matrix, bdists::Matrix)::Tuple{
    Vector{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    tmp = hcat([[rchannel_stats(tdists[:,i],bdists[:,i])...] for i=1:size(bdists,2)]...)
    return tmp[1,:], hcat(tmp[2,:]...), tmp[3,:], tmp[4,:]
end
