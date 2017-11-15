module RchannelExperiments
using Reexport
using Parameters
import Base.Dates: Time

export RcExp, RcRes

# The experiment setup, such as geometry, measurements, etc.
include("setup/RcSetup.jl")
@reexport using .RcSetup
# Image processing to get diameter, scalloping, etc
include("images/RcImages.jl")
@reexport using .RcImages
# All measurements which were recorded by the computer
include("labview/RcLabView.jl")
@reexport using .RcLabView
# volume
include("volume/RcVolume.jl")
@reexport using .RcVolume
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
    # The time periods when water is flowing
    experiment_running::Vector{Tuple{Time,Time}}=Tuple{Time,Time}[]
    # Time of events
    events::Dict{Dates.DateTime,Any}=Dict{Dates.DateTime,Any}()
    setup::ExpSetup
    # sub exps
    expimg::Union{ExpImgs,Void}
    explv::Union{ExpLabView,Void}
    expvol::ExpVol=ExpVol()
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
