#__precompile__()

module RchannelImages

using Images, ImageView, ImageSegmentation
using CoordinateTransformations, OffsetArrays
using ProgressMeter
import PyPlot
import FileIO
import JLD
# const JLD = JLD2
const P = PyPlot
import PyCall
PyCall.@pyimport matplotlib.animation as anim
using StatsBase
using Parameters
import DSP

export ExpImgs, channel_width, prep_img


"""
All parameters needed to turn pixels into meters, assuming
an image without distortion.
"""
@with_kw struct Pixel2Meter @deftype Float64
    h1 # distance surface to pupil of camera
    h2 # thickness of water-pool or perspex
    h3 # depth from bottom of pool to r-channel center
    n1 = 1.0
    n2 = [1.339, 1.48][1] # [water, perspex]
    n3 = 1.31 # ice
    ppm::Float64 # pixel (of full resolution) per meter at surface
end


"""
This holds all data needed to process the automatic pictures
from one Experiment.
"""
@with_kw struct ExpImgs @deftype Int
    dir::String
    p1::Tuple{Int,Int}
    p2::Tuple{Int,Int}
    halfheight_crop = 1200
    thin_num = 4
    color_loc::Tuple{Float64,Float64} = (0.5,0.5) # location of where to pic the color (in original image, coords 0..1)
    ns::Vector{Int}
    minhalfwidth_orig  # minimal half-width of R-channel in pixels of original image
    algo::Symbol=[:both, :thresh, :edge][1]
    gauss_w = 1
    quant::Float64=0.6
    gap = 19
    median_filter_region::Tuple{Int,Int} = (3,7) # (top-bottom, left-right)
    p2m::Pixel2Meter
end
function Base.size(ep::ExpImgs)
    @unpack halfheight_crop, thin_num, p1, p2 = ep
    return fld1(2*halfheight_crop+1,thin_num), fld1(round(Int, sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2)) +1, thin_num)
end

function image_files(ep::ExpImgs)
    @unpack dir, ns = ep
    ["$dir/$f" for f in readdir(dir)][ns]
end
"""
Hold the result from the processing
"""
@with_kw struct ExpImgsResults
    ep::ExpImgs
    center_dist::Vector{Int}
    ts::Matrix{Int} # tops
    bs::Matrix{Int} # bottoms
    ts_dist::Matrix{Float64}
    bs_dist::Matrix{Float64}
    # products
    dia_mean::Vector{Float64} # mean diameter
    dia_quant::Matrix{Float64} # 10%, 50%, 90% quantiles
    scalloping::Vector{Float64} # scalloping factor
    tb_cor::Vector{Float64}
    # images
    thumbs::Vector{Any} # thumbnail images
    imgs::Vector{Any} # full images
end

include("optics.jl")
include("helpers.jl")
include("preparation.jl")
include("width.jl")
include("plotting.jl")
end # module
