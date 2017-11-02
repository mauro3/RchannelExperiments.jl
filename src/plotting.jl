
function imshow_file(fl, ep)
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

function plot_all_n_new(img, top, bottom, top_t, bottom_t, top_e, bottom_e;
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
