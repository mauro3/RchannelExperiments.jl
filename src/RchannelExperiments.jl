module RchannelExperiments
using Reexport, Parameters
import Images
import Base.Dates: Time
import DSP, PyPlot
const P=PyPlot

export RcExp, RcRes, rho, g
const rho = 1000
const g = 9.81
const mu = 0.001792 # dynamic viscosity of water at 0C https://www.thermexcel.com/english/tables/eau_atm.htm
const nu = mu/rho  # kinematic viscosity


include("helpers.jl")
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
    experiment_start::DateTime
    # The time periods when water is flowing
    experiment_running::Vector{Tuple{Time,Time}}=Tuple{Time,Time}[]
    # Time of events
    events::Dict{Any,Any}=Dict{Any,Any}()
    setup::Setup
    # sub exps
    exi::Union{ExpImgs,Void}
    exl::Union{ExpLabView,Void}
    exv::ExpVol=ExpVol()
    ext::ExpTemp=ExpTemp(notebook_pdf="")
    extras::Dict{Symbol}=Dict{Symbol,Any}()
end

"Holds all results"
@with_kw struct RcRes
    ex::RcExp
    resi::Union{ExpImgsResults,Void}
    resl::Union{ExpLabViewResults,Void}
    rest::ExpTempResults=ExpTempResults()
    resv::ExpVolResults=ExpVolResults()
    extras::Dict{Symbol}=Dict{Symbol,Any}()
end

include("plotting.jl")

include("proc.jl")
include("uncertainty.jl")
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
