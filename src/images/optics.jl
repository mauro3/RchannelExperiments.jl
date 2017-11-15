# want:
# - truew -- true width of Rchannel in m
#
# have
# - ppm: pix per m at surface
#

"""
    pixel2meter(pixels, exi::ExpImgs)

This gives the pixel to meter conversion for the plane in which the
R-channel lies.
"""
function pixel2meter(pixels, exi::ExpImgs)
    r1 = pixels/pixels_per_m(exi)
    r1_to_R(r1, exi)
end

pixels_per_m(exi::ExpImgs) = exi.p2m.ppm/exi.thin_num


"""
    r1_to_R(r1, p2m::Pixel2Meter)
    r1_to_R(r1, exi::ExpImgs)

Transforms a distance as measured at the surface of the block into a
distance inside the ice block.

Notes:
- see Yuri's thesis
- green is about 500nm wavelength
- water refractive index for green at 0C is about 1.339
  http://optics.sgu.ru/_media/optics/staff/bashkatov/bashkatov_spie_03_5068_393.pdf
- perspex is about 1.48
"""
function r1_to_R(r1, exi::ExpImgs)
    p2m = exi.p2m
    r1_to_R(r1, p2m.h1, p2m.h2, p2m.h3, p2m.n1, p2m.n2, p2m.n3)
end
function r1_to_R(r1, h1, h2, h3, n1=1.0, n2=[1.339,1.48][1], n3=1.31)
    α = atan(r1/h1)
    β = asin(sin(α)*n1/n2)
    γ = asin(n2/n3*sin(β))
    Δr2 = tan(β)*h2
    r3 = (tan(γ)*h3 + r1 + Δr2) / (1+(tan(γ))^2)
    return r3/cos(γ)
end
