###########
# Image preparation
#
# TODO
# - scale
# - co-registration

"""
    crop(img, topleft, bottomright)

Crops an image.
"""
function crop(img, topleft, bottomright)
    i1, i2 = indices(img)
    img[max(topleft[1],i1[1]):min(bottomright[1],i1[end]),
        max(topleft[2],i2[1]):min(bottomright[2],i2[end])]
end
crop(img, left, right, down, up) = img[1+left:end-right,1+down:end-up]
"""
    thin(img,step)

Thin an image
"""
thin(img,step) = img[1:step:end,1:step:end]

"""
    rotate_n_crop(img, ep::ExpImgs; verbose=false)
    rotate_n_crop(img, p1::Tuple{Int,Int}, p2::Tuple{Int,Int},
                  halfheight::Integer; verbose=false)

Rotates the image such that the line through p1 and p2 is horizontal.
It then crops it such that the width spans p1 and p2 and the
height is 2x halfheight.

Besides the processed image it also returns a list of tuples to compute
the offset in pixels between the new image center line and the old one. I.e.
to correct the tops and bottoms:

    tops .+ round.(Int,midoffset)
    bottomright .- round.(Int,midoffset)

Notes: needs to be in sync with ExpImgs defaults definition.
"""
rotate_n_crop(img, ep::ExpImgs; verbose=false) =
    rotate_n_crop(img, ep.p1, ep.p2, ep.halfheight_crop, verbose=verbose)
function rotate_n_crop(img, p1::Tuple{Int,Int}, p2::Tuple{Int,Int}, halfheight::Integer; verbose=false)
    mid_orig = getmid(img)
    if verbose
        guidict = imshow(img);
        annotate!(guidict, AnnotationPoint(p1[2], p1[1], shape='x', size=20, linewidth=2));
        annotate!(guidict, AnnotationPoint(p2[2], p2[1], shape='x', size=20, linewidth=2));
    end
    @assert p1[2]<=p2[2] "point p1 must be on the left of p2"
    angle = atan2(p2[1]-p1[1], p2[2]-p1[2])
    len = round(Int, sqrt( (p2[1]-p1[1])^2 + (p2[2]-p1[2])^2 ))
    tfm = recenter(RotMatrix(-angle), [p1...])
    img2 = warp(img, tfm);
    topleft = (p1[1]-halfheight, p1[2])
    bottomright = (p1[1]+halfheight, p1[2]+len)
    midoffset = ((p2[1]-p1[1])/(p2[2]-p1[2]), p1[1]-mid_orig)
    if verbose
        guidict = imshow(img2);
        annotate!(guidict, AnnotationPoint(p1[2], p1[1], shape='x', size=20, linewidth=2));
        annotate!(guidict, AnnotationBox(topleft[2], topleft[1], bottomright[2], bottomright[1], linewidth=2, color=RGB(0,0,1)))
    end
    crop(img2, topleft, bottomright), midoffset
end

"""
    prep_img(path::String, ep::ExpImgs; verbose=false)
    prep_img(img_color, ep; verbose=false)

Prepare image by:
- rotate and crop
- colordiff it
"""
function prep_img(img_num::Int, ep::ExpImgs; verbose=false)
    @unpack dir = ep
    path = ["$dir/$f" for f in readdir(dir)][ep.ns[img_num]];
    prep_img(path, ep; verbose=verbose)
end
function prep_img(path::String, ep::ExpImgs; verbose=false)
    verbose && println(path)
    img_color = load(path);
    prep_img(img_color, ep; verbose=verbose)
end
function prep_img(img_color::AbstractArray, ep::ExpImgs; verbose=false)
    @unpack p1, p2, halfheight_crop, thin_num = ep
    img_color, midoffset = rotate_n_crop(img_color, p1, p2, halfheight_crop, verbose=verbose)
    img_color = thin(img_color, thin_num);
    @assert size(img_color)==ep.siz
    # calculate the difference in color
    img = colordiffit(img_color, ep, verbose=verbose);
    return img, midoffset[1]*(1:ep.siz[2]) + midoffset[2]
end

"""
    colordiffit(img_color, ep::ExpImgs; verbose=false)

Use `colordiff` make an image of perceived color difference
"""
function colordiffit(img_color, ep::ExpImgs; verbose=false)
    loc = ep.color_loc
    c0 = img_color[loc...]
    if verbose
        guidict = imshow(img_color)
        dp = size(img_color,1)÷50
        annotate!(guidict, AnnotationBox(loc[2]+dp, loc[1]+dp, loc[2]-dp, loc[1]-dp, linewidth=2, color=RGB(0,0,1)))
    end
    colordiff.(img_color, c0)
end
