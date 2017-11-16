module RcSetup
using Parameters
import DataStructures: OrderedDict
export Setup, Port, all_ports, OrderedDict

module Port
export PressPorts, all_ports
@enum PressPorts in1 in2 in3 in4 out1 out2 out3 out4 ice1 ice2 ice3 ice4
const all_ports = [in1, in2, in3, in4, ice1, ice2, ice3, ice4, out1, out2, out3, out4]
end
"Distance between fixed installed ports"
const def_portdists = [0.88, 0.44, 0.44, 1.68, 0.34, 0.34, 0.68]

# TODO: figure something out for pressure ports in ice:  The change depending on time

const P2C = OrderedDict{Port.PressPorts,Vector{Int}}
@with_kw struct Setup
    drawing_pdf::String # path to pdf
    port2channels::P2C # mapping physical port to channels
                                                            # keep the ports in order
    @assert length(used_channels(port2channels))==length(unique(used_channels(port2channels)))
    portdists::Vector{Float64}=def_portdists # needs to be in order of keys(port2channels)
end
"Return the used measurement channels"
used_channels(di::P2C) = sort(vcat(values(di)...))
used_channels(setup::Setup) = used_channels(setup.port2channels)
"Return the used measurement ports"
used_ports(di::P2C) = collect(keys(di))
used_ports(setup::Setup) = used_ports(setup.port2channels)

end
