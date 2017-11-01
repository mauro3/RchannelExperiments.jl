#__precompile__()

module RchannelImages

using Images, ImageView, ImageSegmentation
using CoordinateTransformations, OffsetArrays
import PyPlot
const P = PyPlot
using StatsBase
using Parameters

include("helpers.jl")

export ExpPics, channel_width

"""
This holds all data needed to process pictures from one Experiment.
"""
@with_kw struct ExpPics @deftype Int
    dir::String
    p1::Tuple{Int,Int}
    p2::Tuple{Int,Int}
    halfheight = 1200
    thin_num = 4
    siz::Tuple{Int,Int} = (fld1(2*halfheight+1,thin_num),
                           fld1(round(Int, sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2)) +1, thin_num) )
    color_loc::Tuple{Int,Int} = (fld1(siz[1],2), fld1(siz[2],2))
    ns::Vector{Int}
    minhalfwidth
    gauss_w = 1
    quant::Float64=0.6
    gap = 19
    median_filter_region::Tuple{Int,Int} = (0,0)
end

include("preparation.jl")
include("width.jl")
include("plotting.jl")
end # module
