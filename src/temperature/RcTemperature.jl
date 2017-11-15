module RcTemperature
using Parameters
export ExpTemp, ExpTempResults

@with_kw struct ExpTemp
    notebook_pdf::String # path to pdf-scan of notebook containing the measurements
end

struct ExpTempResults
end

end
