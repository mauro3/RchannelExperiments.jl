
function imshow_file(fl, ep)
    @unpack p1,p2 = ep
    img = Images.load(fl);
    guidict = ImageView.imshow(img);
    ImageView.annotate!(guidict, ImageView.AnnotationPoint(p1[2], p1[1], shape='x', size=20, linewidth=2));
    ImageView.annotate!(guidict, ImageView.AnnotationPoint(p2[2], p2[1], shape='x', size=20, linewidth=2));
end


function plot_lines_on_img(img, top, bottom;
                          col="r", label="", ax=nothing,
                          legend=true)
    ax==nothing || P.axes(ax)
    mid = getmid(img)
    greys =  P.get_cmap("Greys")
    P.imshow(img, aspect="equal", extent=(1,size(img,2),size(img,1), 1), cmap=greys)
    top = top2ind(top, img)
    bottom = bottom2ind(bottom, img)
    P.plot(1:length(top), top, c=col, lw=1, label="$label top")
    P.plot(1:length(top), bottom, c=col, lw=1, label="$label bottom")
    legend && P.legend()
end

"""
Two subplots with one made to show all lines and the other just the new one.
Both overlain by the current image.
"""
function plot_all_n_new(img, top, bottom, top_t=top, bottom_t=bottom, top_e=top_t, bottom_e=bottom_t;
                        col="r", label="", title="",
                        legend=false)
    ax = P.subplot(2,1,1)
    plot_lines_on_img(img, top, bottom,
                     col=col, label=label, ax=ax,
                     legend=legend)
    P.title(title)
    ax = P.subplot(2,1,2)
    P.cla()
    plot_two_lines(img, top_e, bottom_e, top_t, bottom_t)
    plot_lines_on_img(img, top, bottom,
                      col=col, label=label, ax=ax,
                      legend=legend)
end
function plot_all_n_new(ep::ExpImgs, tops, bottoms;
                        col="r", label="", title="",
                        legend=false)
    P.figure()
    for i=1:length(ep.ns)
        img, _, _ = prep_img(i, ep)
        plot_all_n_new(img, tops[:,i], bottoms[:,i],
                       col=col, label=label, title=title)
        P.draw()
        sleep(0.01)
    end
end

function plot_two_lines2(img, top_e, bottom_e, top_t, bottom_t; title="")
    #P.figure()
    ax = P.subplot(2,1,1)
    plot_lines_on_img(img, top_e, bottom_e)
    P.ylabel("edge")
    P.title(t)
    P.subplot(2,1,2, sharex=ax, sharey=ax)
    plot_lines_on_img(img, top_t, bottom_t)
    P.ylabel("thresh")
end

function plot_two_lines(img, top_e, bottom_e, top_t, bottom_t; title="", ax=nothing)
    ax==nothing || P.axes(ax)
    plot_lines_on_img(img, top_e, bottom_e, col="b", label="edge")
    plot_lines_on_img(img, top_t, bottom_t, col="g", label="thresh")
    P.title(title)
end

function plot_proms(mag, i)
    peaks,locs,proms,doms = findpeaks(mag[:,i]);
    peaks2,locs2,proms2,doms2 = findpeaks(peaks);
    peaks3,locs3,proms3,doms3 = findpeaks(peaks2);
    peaks4,locs4,proms4,doms4 = findpeaks(peaks3);
    P.figure()
    P.plot(1:size(mag,1), mag[:,i])
    P.plot(locs,peaks)
    P.plot(locs[locs2], peaks2)
    P.plot(locs[locs2[locs3]], peaks3)
    P.plot(locs[locs2[locs3[locs4]]], peaks4)
end



# https://brushingupscience.wordpress.com/2016/06/21/matplotlib-animations-the-easy-way/
# https://genkuroki.github.io/documents/Jupyter/20170624%20Examples%20of%20animations%20in%20Julia%20by%20PyPlot%20and%20matplotlib.animation.html
function animate_res(res::ExpImgsResults; image=[:none,:thumb,:img][1], save=false, interval=200)
    image==:img && error("Not implemented")
    @unpack ep, ts, bs, thumbs, imgs = res
    scale = ep.thin_num/ep.p2m.ppm *100 # cm
    min_alpha = 0.2
    fac_alpha = 0.8
    c = "b"
    ts = top2ind(ts, ep)*scale
    bs = bottom2ind(bs, ep)*scale
    fig, ax = P.subplots(figsize=(5, 3), dpi=100)
    greys =  P.get_cmap("Greys")
    ax[:set](xlim=(0, size(ep)[2]*scale), ylim=(0, size(ep)[1]*scale))
    ax[:set](xlabel="Along channel (cm)", ylabel="Across channel (cm)")
    ax[:invert_yaxis]()

    if image==:thumb
        im = ax[:imshow](thumbs[1], extent=(0,size(ep)[2]*scale,size(ep)[1]*scale,0),
                         aspect="equal", cmap=greys)
    end
    pts = (0:size(ep)[2]-1)*scale
    lines_t = ax[:plot](pts,ts[:,1], c=c)
    lines_t[end][:set_alpha](1)
    lines_b = ax[:plot](pts,bs[:,1], c=c)
    lines_b[end][:set_alpha](1)
    ii = 1
    function ani_frame(i)
        if image==:thumb
            im[:set_array](thumbs[i])
        end
        for l in lines_t
            a = l[:get_alpha]()
            l[:set_alpha](max(min_alpha, a*fac_alpha))
        end
        for l in lines_b
            a = l[:get_alpha]()
            l[:set_alpha](max(min_alpha, a*fac_alpha))
        end
        append!(lines_t, ax[:plot](pts,ts[:,i], c=c))
        append!(lines_b, ax[:plot](pts,bs[:,i], c=c))
        lines_t[end][:set_alpha](1)
        lines_b[end][:set_alpha](1)
        ii +=1
        ii>1000 && error()
        [lines_t;lines_b]
    end
    myanim = anim.FuncAnimation(fig, ani_frame, frames=2:size(ts,2), interval=interval, repeat=false)
    fln = filename(ep)*".mp4"
    save && myanim[:save](fln, bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end

function plot_res(res::ExpImgsResults; tscale=60)
    P.figure()
    t, dia_mean, dia_quant, scalloping, tb_cor = get_time_series(res)
    dia_std = std(res.ts_dist+res.bs_dist,1)'
    scalloping2 = dia_std./dia_mean
    t ./= tscale
    ax = P.subplot(2,1,1)
    P.plot(t, dia_mean, label="mean")
    P.plot(t, dia_quant[:,2], label="median")
    P.plot(t, dia_quant[:,1], label="10%")
    P.plot(t, dia_quant[:,3], label="90%")
    P.legend()
    P.ylabel("diameter (m)")
    P.subplot(2,1,2,sharex=ax)
    P.plot(t, scalloping, label="scallop quantile (0..2)")
    P.plot(t, scalloping2, label="scallop std/mean (0..2)")
    P.plot(t, tb_cor, label="symmetry (correlation) (-1..1)")
    P.legend()
    st = res.ep.experiment_start
    P.xlabel("Time period since $(Dates.Time(st)) on $(Dates.Date(st)) ($(tscale)s)")
    P.ylabel("factor ()")
end
