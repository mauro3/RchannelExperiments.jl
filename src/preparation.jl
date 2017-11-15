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
function thin(img,step, alg=[:imresize, :restrict, :subsampling][3])
    if alg==:imresize
        ## using filtering: a bit slow
        sz = (length(1:step:size(img,1)), length(1:step:size(img,2)))
        σ = maximum(map((o,n)->0.25*o/n, size(img), sz))
        kern = KernelFactors.gaussian((σ,σ))   # from ImageFiltering
        return imresize(imfilter(img, kern, NA()), sz)
    elseif alg==:restrict
        ### using restrict, slightly faster
        l2 = log2(step)
        @assert isinteger(l2) "Can only thin by power of 2."
        for i=1:l2
            s = size(img)
            img = ImageFiltering.padarray(img, ImageFiltering.Pad(:reflect,1,1)) # does not work
            img = ImageTransformations.restrict(img)[1:fld1(s[1],2),1:fld1(s[2],2)]
        end
        return img
    elseif alg==:subsampling
        ### just sub-sampling: causes aliasing
        return img[1:step:end,1:step:end]
    end
    error()
end

"""
    rotate_n_crop(img, ep::ExpImgs; verbose=false)
    rotate_n_crop(img, p1::Tuple{Int,Int}, p2::Tuple{Int,Int},
                  halfheight::Integer; verbose=false)

Rotates the image such that the line through p1 and p2 is horizontal.
It then crops it such that the width spans p1 and p2 and the
height is 2x halfheight.

Besides the processed image it also returns the distance of the (p1,p2) line to
the center of the original image.  This is needed to calculate real distances:
`top+center_dist` and `bottom-center_dist` gives the number of pixels from a line
through the center-point.

Notes: needs to be in sync with ExpImgs defaults definition.
"""
rotate_n_crop(img, ep::ExpImgs; verbose=false) =
    rotate_n_crop(img, ep.p1, ep.p2, ep.halfheight_crop, verbose=verbose)
function rotate_n_crop(img, p1::Tuple{Int,Int}, p2::Tuple{Int,Int}, halfheight::Integer; verbose=false)
    mid_point = map(x->fld1(x,2), size(img))
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
    if verbose
        guidict = imshow(img2);
        annotate!(guidict, AnnotationPoint(p1[2], p1[1], shape='x', size=20, linewidth=2));
        annotate!(guidict, AnnotationBox(topleft[2], topleft[1], bottomright[2], bottomright[1], linewidth=2, color=RGB(0,0,1)))
    end
    center_dist = cross([(mid_point.-p1)...,0], [(p2.-p1)...,0])[3]/norm([(p1.-p2)...])
    crop(img2, topleft, bottomright), round(Int, center_dist)
end

"""
    prep_img(path::String, ep::ExpImgs; verbose=false)
    prep_img(img_color, ep; verbose=false)

Load and prepare image by:
- rotate and crop
- colordiff it

Also returns:
- center_dist -- distance to previous center
- capture_time -- time of image capture
"""
function prep_img(path::String, ep::ExpImgs; verbose=false)
    verbose && println(path)
    img_color = load(path);

    dt = ImageMagick.magickinfo(path, "exif:DateTimeOriginal")["exif:DateTimeOriginal"]
    _prep_img(img_color, ep; verbose=verbose)..., DateTime(dt, "yyyy:mm:dd HH:MM:SS") - ep.time_correction
end
function _prep_img(img_color::AbstractArray, ep::ExpImgs; verbose=false)::Tuple{Matrix{Float64}, Int}
    @unpack p1, p2, halfheight_crop, thin_num = ep
    img_color, center_dist = rotate_n_crop(img_color, p1, p2, halfheight_crop, verbose=verbose)
    img_color = thin(img_color, thin_num);
    @assert size(img_color)==size(ep) "size(img_color)=$(size(img_color)) not equal size(ep)=$(size(ep))"
    # calculate the difference in color
    img = colordiffit(img_color, ep, verbose=verbose);
    # find EXIF time and correct it
    return img, center_dist
end

"""
    colordiffit(img_color, ep::ExpImgs; verbose=false)

Use `colordiff` make an image of perceived color difference
"""
function colordiffit(img_color, ep::ExpImgs; verbose=false)
    loc1, loc2 = round.(Int, (size(img_color,1)*ep.color_loc[1], size(img_color,2)*ep.color_loc[2]))
    c0 = img_color[loc1,loc2]
    if verbose
        guidict = imshow(img_color)
        dp = size(img_color,1)÷50
        annotate!(guidict, AnnotationBox(loc2+dp, loc1+dp, loc2-dp, loc1-dp, linewidth=2, color=RGB(0,0,1)))
    end
    colordiff.(img_color, c0)
end
