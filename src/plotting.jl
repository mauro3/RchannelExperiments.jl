
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

function plot_all_n_new(img, top, bottom, top_t=top, bottom_t=bottom, top_e=top_t, bottom_e=bottom_t;
                        col="r", label="", ax=nothing, title="",
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
function animate_line(tops, bottoms, ep)
    min_alpha = 0.2
    fac_alpha = 0.8
    c = "b"
    tops = top2ind(tops, ep)
    bottoms = bottom2ind(bottoms, ep)
    fig, ax = P.subplots(figsize=(5, 3))
    ax[:set](xlim=(1, ep.siz[2]), ylim=(1, ep.siz[1]))
    ax[:invert_yaxis]()

    pts = 1:ep.siz[2]
    lines_t = ax[:plot](pts,tops[:,1], c=c)
    lines_t[end][:set_alpha](1)
    lines_b = ax[:plot](pts,bottoms[:,1], c=c)
    lines_b[end][:set_alpha](1)
    ii = 1
    function ani_frame(i)
        for l in lines_t
            a = l[:get_alpha]()
            l[:set_alpha](max(min_alpha, a*fac_alpha))
        end
        for l in lines_b
            a = l[:get_alpha]()
            l[:set_alpha](max(min_alpha, a*fac_alpha))
        end
        append!(lines_t, ax[:plot](pts,tops[:,i], c=c))
        append!(lines_b, ax[:plot](pts,bottoms[:,i], c=c))
        lines_t[end][:set_alpha](1)
        lines_b[end][:set_alpha](1)
        ii +=1
        ii>1000 && error()
        [lines_t;lines_b]
    end
    myanim = anim.FuncAnimation(fig, ani_frame, frames=2:size(tops,2), interval=200, repeat=false)
    fln = filename(ep)*".mp4"
    myanim[:save](fln, bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end
