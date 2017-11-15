using RchannelExperiments
using Base.Test
const R=RchannelExperiments
const RI=RchannelExperiments.RcImages
# write your own tests here

x = linspace(-5pi,5pi,156);
vec = sin.(x).*cos.(x/5)+0.005x.^2;
peaks, ips, proms, doms = RI.findpeaks(vec);
@test peaks ≈ [1.96984, 0.497992, 0.709033, 0.964469, 0.362738, 1.24738] atol=1e-5
@test maximum(peaks)==maximum(vec)
@test all(doms.==[148, 17, 27, 45, 10, 121])
@test proms ≈ [2.9112, 0.369982, 0.694124, 1.47311, 0.10339, 2.18874] atol=1e-5

x = linspace(-5pi,5pi,156);
vec = sin.(x).*cos.(x/5)+0.005x.^2;
peaks, ips, proms, doms = RI.findpeaks(vec);
peaks_l, ips_l, proms_l, doms_l = RI.findpeaks_onesided(vec, :left);
peaks_r, ips_r, proms_r, doms_r = RI.findpeaks_onesided(vec, :right);
@test peaks ≈ peaks_r
@test peaks ≈ peaks_l
rs_d = [1,2,3,4,5]
ls_d = [6]
rs_p = [1,2,4,5]
ls_p = [3,6]
@test proms[rs_p] ≈ proms_r[rs_p]
@test proms[ls_p] ≈ proms_l[ls_p]
@test doms[rs_d] ≈ doms_r[rs_d]
@test doms[ls_d] ≈ doms_l[ls_d]
