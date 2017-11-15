module RcSetup
using Parameters
export ExpSetup, Port, all_ports

module Port
export PressPorts, all_ports
@enum PressPorts in1 in2 in3 in4 out1 out2 out3 out4 ice1 ice2 ice3 ice4
const all_ports = [in1, in2, in3, in4, out1, out2, out3, out4, ice1, ice2, ice3, ice4]
end

# TODO: figure something out for pressure ports in ice
@with_kw struct ExpSetup
    drawing_pdf::String # path to pdf
    port2channels::Dict{Port.PressPorts,Vector{Int}} # mapping physical port to channels
    @assert length(used_channels(port2channels))==length(unique(used_channels(port2channels)))
end
used_channels(di::Dict{Port.PressPorts,Vector{Int}}) = sort(vcat(values(di)...))
used_channels(ex::ExpSetup) = used_channels(ex.port2channels)

end
