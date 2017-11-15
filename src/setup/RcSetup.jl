module RcSetup
using Parameters
export ExpSetup

@with_kw struct ExpSetup
    drawing_pdf::String # path to pdf
end

end
