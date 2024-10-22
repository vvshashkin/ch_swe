import numpy as np
import re
import Ngl
from sys import argv
from add_cubed_sphere import add_cubed_sphere

path = argv[1]

schemes = ["Ch_21","Ch_42","Ch_63", "Ah_63"]
#non-linear test:
#Ns = [20,40,80,160]
#dts = {20 : 800, 40 : 400, 80 : 200, 160 : 100}
#map_Ns = [20, 40]
#linear test
Ns = [48,64,96,128,192]
dts = {32: 900, 48: 600, 64: 450, 96: 300, 128: 225, 192: 150}
map_Ns = [48, 64]

resolutions=["N{:03d}_dt{:04d}".format(N, dts[N]) for N in Ns]

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

cn_res.nglDraw = False
cn_res.nglFrame = False

overlay_res = Ngl.Resources()
overlay_res.nglDraw = False
overlay_res.nglFrame = False
overlay_res.cnLineThicknessF = 2.0
overlay_res.cnInfoLabelOn = False

for N in map_Ns:
    Nlon = 4*N
    Nlat = 2*N+1
    cn_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    cn_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
    overlay_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
    overlay_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
    wksres = Ngl.Resources()
    wksres.wkOrientation="portrait"
    wks = Ngl.open_wks("pdf", "ts2_h_error_N{:03d}".format(N),wksres)
    Ngl.define_colormap(wks,"GMT_polar")
    plots = []
    for scheme in schemes:
        fd = open("h_N"+"{:03d}_".format(N)+"dt{:04d}_".format(dts[N])+scheme+".dat","rb")
        h0 = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape(Nlat,Nlon)
        fd.seek(4*Nlon*Nlat*10,0)
        h1 = np.fromfile(fd,count=Nlon*Nlat,dtype=np.float32).reshape(Nlat,Nlon)
        fd.close()
        cn_res.tiMainString = scheme
        plots.append(Ngl.contour_map(wks,h1-h0,cn_res))
        overlay_plot = Ngl.contour(wks, h0, overlay_res)
        Ngl.overlay(plots[-1], overlay_plot)
        add_cubed_sphere(wks,plots[-1])

    textres = Ngl.Resources()
    textres.txFontHeightF = 0.020
    #Ngl.text_ndc(wks,"Height error at day 10 [m]",0.5,.97,textres)
    pres = Ngl.Resources()
    Ngl.panel(wks,plots,(2,2),pres)
    Ngl.delete_wks(wks)

def read_err(fname):
    fd = open(fname,"r")
    l2_h = []
    linf_h = []
    for line in fd.readlines():
        if re.search("Errors,",line):
            errs = re.findall("\d+\.\d+E.\d+|\d+\.\d+",line)
            l2_h.append(float(errs[1]))
            linf_h.append(float(errs[2]))
    return np.array(l2_h), np.array(linf_h)

# l221, linf = read_err(path+"/errors_N160_dt100_Ch_21.txt")
# l242, linf = read_err(path+"/errors_N160_dt100_Ch_42.txt")
# l263, linf = read_err(path+"/errors_N160_dt100_Ch_63.txt")


# wks = Ngl.open_wks("pdf", "l2_h_t")

# res = Ngl.Resources()
# res.trYLog = True

# plot = Ngl.xy(wks,range(1,241),np.array([l221,l242,l263]),res)

# Ngl.delete_wks(wks)

l2_conv = []
linf_conv = []
path = path+"/swm_"

for scheme in schemes:
    for res in resolutions:
        fname = path+res+"_"+scheme+".out"
        l2, linf = read_err(fname)
        l2_conv.append(np.max(l2))
        linf_conv.append(np.max(linf))
l2_order = []
linf_order=[]
x = np.array(Ns)
l2_0 = {2: 1e1, 3: 2e0, 4: 1e-1}
linf_0 = {2: 2e1, 3: 2e0, 4: 4e-1}
for order in [2,3,4]:
    for i in range(len(x)):
        l2_order.append(l2_0[order]*(1.0*x[0]/x[i])**order)
        linf_order.append(linf_0[order]*(1.0*x[0]/x[i])**order)

l2 = np.array(l2_conv).reshape((len(schemes),len(resolutions)))
linf = np.array(linf_conv).reshape((len(schemes),len(resolutions)))
l2_ord = np.array(l2_order).reshape((3,len(resolutions)))
linf_ord = np.array(linf_order).reshape((3,len(resolutions)))

wks = Ngl.open_wks("pdf", "ts2_conv")
res = Ngl.Resources()
res.trXLog = True
res.trYLog = True
res.tiXAxisString="N~B~c"
res.tmXBMode = "Explicit"
res.tmXBValues = x
res.tmXBLabels = x
res.trXMinF = 32
res.trXMaxF = 256
res.nglDraw = False
res.nglFrame = False
res.xyLineThicknessF  = 2
res.xyLineColors=["red","green","blue","orange"]
res.xyMarkLineMode="MarkLines"
res.xyMarker = 16

res.tiYAxisString="~F10~l~B~2~N~ ~F21~error"
plot1 =  Ngl.xy(wks,x,l2,res)
res.tiYAxisString="~F10~l~B~~F34~%~N~~F21~ error"
plot2 =  Ngl.xy(wks,x,linf,res)

lres = Ngl.Resources()
text_label_orders=['2nd order', '3rd order', '4th order']
for i in range(3):
    lres.gsLineLabelString = text_label_orders[i]
    Ngl.add_polyline(wks, plot1, x, l2_ord[i,:],lres)
    Ngl.add_polyline(wks, plot2, x, linf_ord[i,:],lres)

lg_res = Ngl.Resources()
lg_res.vpWidthF = 0.12
lg_res.vpHeightF = 0.11
lg_res.lgLineThicknessF  = 2.0
lg_res.lgLineColors = ["red","green","blue", "orange"][::-1]
lg_res.lgDashIndexes = [0,0,0,0]
lg_res.lgLabelFontHeightF = 0.013
lg_res.lgPerimOn = True
lg1 = Ngl.legend_ndc(wks,len(schemes),schemes[::-1],0.1,0.45,lg_res)
lg2 = Ngl.legend_ndc(wks,len(schemes),schemes[::-1],0.6,0.45,lg_res)

Ngl.panel(wks,(plot1,plot2),(1,2))

Ngl.delete_wks(wks)

def linreg(f,x):
    a = (np.sum(f*x)-np.sum(f)*np.sum(x)/len(x)) / (np.sum(x*x)-np.sum(x)**2/len(x))
    c = np.sum(f-a*x) / len(x)
    return a,c

x = np.log(Ns)
f2 = np.log(l2)
finf = np.log(linf)
for i in range(len(schemes)):
    a2,c2 = linreg(f2[i,:],x)
    ai,ci = linreg(finf[i,:],x)
    print(schemes[i], "conv order:", a2, ai)

Ngl.end()
