module RcLabView
using Parameters
export ExpLabView, ExpLabViewResults

@with_kw struct ExpLabView
    file::String # output file of Lab computer
end
struct ExpLabViewResults
end

end
