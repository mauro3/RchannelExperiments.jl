module RchannelExperiments
using Reexport
using Parameters

# The experiment setup, such as geometry, measurements, etc.
include("setup/RcSetup.jl")
@reexport using .RcSetup
# Image processing to get diameter, scalloping, etc
include("images/RcImages.jl")
@reexport using .RcImages
# All measurements which were recorded by the computer
include("labview/RcLabView.jl")
@reexport using .RcLabView
# hand-recorded temperature
include("temperature/RcTemperature.jl")
@reexport using .RcTemperature

"""
This is one Experiment.
"""
@with_kw struct RcExp
    name::String # probably just the ISO date: 2017-10-06
    # start-time of experiment.  Plotting will be relative to this time
    experiment_start::Dates.DateTime
    setup::ExpSetup
    # sub exps
    expimg::ExpImgs
    explv::ExpLabView
    exptemp::ExpTemp=ExpTemp(notebook_pdf="")
    extras::Dict{Symbol}=Dict{Symbol,Any}()
end

"Holds all results"
@with_kw struct RcRes
    rex::RcExp
    resi::ExpImgsResults
    reslv::ExpLabViewResults
    rest::ExpTempResults=ExpTempResults()
    extras::Dict{Symbol}=Dict{Symbol,Any}()
end
# get relevant data
function get_Q(res::RcRes)
    return t, Q
end
function get_D(res::RcRes)
    return t, D
end
function get_p(res::RcRes)
    return t, p
end
function calc_f(res::RcRes)
end
function calc_heat_flux(res::RcRes)
end


end
