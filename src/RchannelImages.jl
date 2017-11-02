#__precompile__()

module RchannelImages

@time using Images, ImageView, ImageSegmentation
@time using CoordinateTransformations, OffsetArrays
@time using ProgressMeter
@time import PyPlot
const P = PyPlot
@time using StatsBase
@time using Parameters

include("helpers.jl")

export ExpImgs, channel_width, prep_img

"""
This holds all data needed to process pictures from one Experiment.
"""
@with_kw struct ExpImgs @deftype Int
    dir::String
    p1::Tuple{Int,Int}
    p2::Tuple{Int,Int}
    halfheight_crop = 1200
    thin_num = 4
    siz::Tuple{Int,Int} = (fld1(2*halfheight_crop+1,thin_num),
                           fld1(round(Int, sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2)) +1, thin_num) )
    color_loc::Tuple{Int,Int} = (fld1(siz[1],2), fld1(siz[2],2))
    ns::Vector{Int}
    minhalfwidth_orig  # minimal half-width of R-channel in pixels of original image
    gauss_w = 1
    quant::Float64=0.6
    gap = 19
    median_filter_region::Tuple{Int,Int} = (3,7) # (top-bottom, left-right)
end

include("preparation.jl")
include("width.jl")
include("plotting.jl")
end # module
