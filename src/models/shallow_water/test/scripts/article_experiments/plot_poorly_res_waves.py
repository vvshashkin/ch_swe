import Ngl
import numpy as np
from sys import argv
from add_cubed_sphere import add_cubed_sphere


ref_file = argv[1]

schemes = ["Ch_21", "Ch_42", "Ch_63", "Ah_63"]
dts = {16:"1800",32:"0900", 48: "0600", 64:"0450", 96:"0300", 128:"0225",192:"0150"}

map_N = 64
map_Nt = 25
avg_Nt = 600
#plot map plots

Nlon = 4*map_N
Nlat = 2*map_N+1

wks = Ngl.open_wks("pdf","poorly_resolved_map")
Ngl.define_colormap(wks,"GMT_polar")

cn_res = Ngl.Resources()
cn_res.cnFillOn = True
cn_res.cnLinesOn = False
cn_res.cnLineLabelsOn = False
cn_res.mpCenterLonF = 180.0
cn_res.mpGridAndLimbOn = False
cn_res.lbOrientation="Horizontal"
cn_res.lbLabelFontHeightF = 0.015
cn_res.tiMainFontHeightF = 0.02
cn_res.mpGeophysicalLineColor = "Transparent"
cn_res.mpGreatCircleLinesOn = True
cn_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
cn_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
cn_res.cnLevelSelectionMode = "ExplicitLevels"
cn_res.cnLevels = np.linspace(-0.15,0.15,16)
cn_res.lbLabelBarOn = False

cn_res.nglDraw = False
cn_res.nglFrame = False

# ref_file = ref_sol_path+"gauss_h_{:03d}".format(map_N)+corner+".dat"
# h_exact = np.fromfile(ref_file, count=Nlon*Nlat, dtype=np.float32, offset=4*Nlon*Nlat*map_Nt).reshape((Nlat,Nlon))

# print(np.min(h_exact), np.max(h_exact))

files=["h_N{:03d}_dt".format(map_N)+dts[map_N]+"_"+scheme+".dat" for scheme in schemes]
files.append("h_N064_dt0450_Ch_63_orig.dat")
# files.append("h_N064_dt0450_Ch_63.dat")
files.append(ref_file)

titles = [scheme for scheme in schemes]
titles.append("Ch_63_orig")
titles.append("Exact")

plots = []
for i in range(len(files)):
    h_file = files[i]
    h = np.fromfile(h_file, count=Nlon*Nlat, dtype=np.float32, offset=4*Nlon*Nlat*map_Nt).reshape((Nlat,Nlon))

    cn_res.tiMainString = titles[i]
    plots.append(Ngl.contour_map(wks,h,cn_res))
    add_cubed_sphere(wks,plots[-1])

pres = Ngl.Resources()
pres.nglPanelLabelBar = True
Ngl.panel(wks,plots,(3,2),pres)
Ngl.delete_wks(wks)

wks = Ngl.open_wks("pdf","poorly_resolved_map_stationary")
Ngl.define_colormap(wks,"WhiteYellowOrangeRed")

cn_res.cnLevelSelectionMode = "AutomaticLevels"
cn_res.lbLabelBarOn = True

plots = []
for i in range(len(files)):
    h_file = files[i]
    h = np.fromfile(h_file, count=avg_Nt*Nlon*Nlat, dtype=np.float32).reshape((avg_Nt,Nlat,Nlon))
    h = np.average(h,axis=0)

    cn_res.tiMainString = titles[i]
    plots.append(Ngl.contour_map(wks,h,cn_res))
    add_cubed_sphere(wks,plots[-1])

pres = Ngl.Resources()
Ngl.panel(wks,plots,(3,2),pres)
Ngl.delete_wks(wks)

Ngl.end()
