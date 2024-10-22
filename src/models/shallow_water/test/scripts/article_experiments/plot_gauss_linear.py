import Ngl
import numpy as np
from sys import argv
from add_cubed_sphere import add_cubed_sphere


Nt = 601
corner=""
ref_sol_path = argv[1]
if len(argv)>2: corner = argv[2]

if(corner != ""): corner = "_"+corner

schemes = ["Ch_21", "Ch_42", "Ch_63", "Ah_63"]
# dts = {32:"1800", 48: "1200", 64:"0900", 96:"0600", 128:"0450"}
dts = {16:"1800",32:"0900", 48: "0600", 64:"0450", 96:"0300", 128:"0225",192:"0150"}
# Ns = [32,48,64,128]
Ns = [48,64,96,192]

map_N = 64
map_Nt = 600

def get_max_err(h_file,N,ref_sol_path,corner):
    Nlon = 4*N
    Nlat = 2*N+1

    w = np.cos(np.linspace(-0.5*np.pi,0.5*np.pi,Nlat)).reshape((Nlat,1))

    l2_max = 0.0
    linf_max = 0.0

    h = np.fromfile(h_file, count=Nlon*Nlat*Nt, dtype=np.float32).reshape((Nt,Nlat,Nlon))
    ref_file = ref_sol_path+"gauss_h_{:03d}".format(N)+corner+".dat"
    h_exact = np.fromfile(ref_file, count=Nlon*Nlat*Nt, dtype=np.float32).reshape((Nt,Nlat,Nlon))
    for i in range(Nt):
        e = h[i,:,:]-h_exact[i,:,:]

        linf = np.max(np.abs(e))# / np.max(h_exact[i,:,:])
        linf_max = max(linf,linf_max)

        l2 = np.sqrt(np.sum(e**2*w) / np.sum(w) / Nlon )# / \
             #np.sqrt(np.sum(h_exact[i,:,:]**2*w) / np.sum(w) / Nlon ) 
        l2_max = max(l2,l2_max)
    return l2_max, linf_max

def get_max_err2(out_file):
    import re
    l2_max = 0.0
    linf_max = 0.0

    fd = open(out_file,"r")
    for ln in fd.readlines():
        if re.search("Errors,",ln):
            errs = re.findall("\d+\.\d+E.\d+|\d+\.\d+",ln)
            l2_max = max(l2_max,float(errs[1]))
            linf_max = max(linf_max,float(errs[2]))

    return l2_max, linf_max

l2_all = []
linf_all = []

for scheme in schemes:
    for N in Ns:
        # h_file = "h_N{:03d}_dt".format(N)+dts[N]+"_"+scheme+corner+".dat"
        # l2, linf = get_max_err(h_file, N, ref_sol_path, corner)
        l2, linf = get_max_err2("swm_N{:03d}_dt".format(N)+dts[N]+"_"+scheme+corner+".out")
        print(scheme,N,l2)
        l2_all.append(l2)
        linf_all.append(linf)

l2_all = np.array(l2_all).reshape(len(schemes),len(Ns))
linf_all = np.array(linf_all).reshape(len(schemes),len(Ns))


def linreg(f,x):
    a = (np.sum(f*x)-np.sum(f)*np.sum(x)/len(x)) / (np.sum(x*x)-np.sum(x)**2/len(x))
    c = np.sum(f-a*x) / len(x)
    return a,c

x = np.log(Ns)
f2 = np.log(l2_all)
finf = np.log(linf_all)

a2 = []
c2 = []
ainf = []
cinf = []

for i in range(len(schemes)):
    a,c = linreg(f2[i,:],x)
    a2.append(a)
    c2.append(c)
    a,c = linreg(finf[i,:],x)
    ainf.append(a)
    cinf.append(c)
    print(schemes[i], "conv order:", a2[-1], ainf[-1])

wks = Ngl.open_wks("pdf","gauss_linear"+corner)
res = Ngl.Resources()
res.trXLog = True
res.trYLog = True
res.tiXAxisString="N~B~c"
res.tmXBMode = "Explicit"
res.tmXBValues = Ns
res.tmXBLabels = Ns
res.trXMinF = 32
res.trXMaxF = 256
res.nglDraw = False
res.nglFrame = False
res.xyLineThicknessF  = 2
res.xyLineColors=["red","green","blue","orange"]
res.xyMarkLineMode="MarkLines"
res.xyMarker = 16

res.tiYAxisString="~F10~l~B~2~N~ ~F21~error"
plot1 =  Ngl.xy(wks,Ns,l2_all,res)
res.tiYAxisString="~F10~l~B~~F34~%~N~~F21~ error"
plot2 =  Ngl.xy(wks,Ns,linf_all,res)

orders = [2,4,5]
order_labels = ['2nd order', '4th order', '5th order']
l2_order = []
linf_order=[]
x = np.array(Ns)
l2_0 = {2: 2e-2, 3: 0.9e-2, 4: 2e-3, 5:1e-4}
linf_0 = {2: 2e-1, 3: 9e-2, 4: 2e-2, 5:1e-3}
for order in orders:
    for i in range(len(x)):
        l2_order.append(l2_0[order]*(1.0*x[0]/x[i])**order)
        linf_order.append(linf_0[order]*(1.0*x[0]/x[i])**order)

l2_ord = np.array(l2_order).reshape((len(orders),len(x)))
linf_ord = np.array(linf_order).reshape((len(orders),len(x)))

lres = Ngl.Resources()
text_label_orders=order_labels
for i in range(3):
    lres.gsLineLabelString = text_label_orders[i]
    Ngl.add_polyline(wks, plot1, x, l2_ord[i,:],lres)
    Ngl.add_polyline(wks, plot2, x, linf_ord[i,:],lres)

lg_res = Ngl.Resources()
lg_res.vpWidthF = 0.12
lg_res.vpHeightF = 0.11
lg_res.lgLineThicknessF  = 2.0
lg_res.lgLineColors = ["red","green","blue","orange"][::-1]
lg_res.lgDashIndexes = [0,0,0,0]
lg_res.lgLabelFontHeightF = 0.013
lg_res.lgPerimOn = True
leg = ["Ch_21","Ch_42", "Ch_63", "Ah_63"]
lg1 = Ngl.legend_ndc(wks,len(leg),leg[::-1],0.1,0.45,lg_res)
lg2 = Ngl.legend_ndc(wks,len(leg),leg[::-1],0.6,0.45,lg_res)

Ngl.panel(wks,(plot1,plot2),(1,2))

Ngl.delete_wks(wks)

#plot map plots

Nlon = 4*map_N
Nlat = 2*map_N+1

wks = Ngl.open_wks("pdf","gauss_linear_map"+corner)
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

cn_res.nglDraw = False
cn_res.nglFrame = False

overlay_res = Ngl.Resources()
overlay_res.nglDraw = False
overlay_res.nglFrame = False
overlay_res.cnLineThicknessF = 2.0
overlay_res.cnInfoLabelOn = False
overlay_res.sfXArray = np.linspace(0.0,360.0,Nlon,endpoint=False)
overlay_res.sfYArray = np.linspace(-90.0,90.0,Nlat)
overlay_res.cnLevelSelectionMode = "ExplicitLevels"
overlay_res.cnLevels = np.arange(-1.0,1.0,0.1)
overlay_res.cnMonoLineDashPattern = False
overlay_res.cnLineDashPatterns = np.where(overlay_res.cnLevels < 0.0, 1,0)

plots = []

ref_file = ref_sol_path+"gauss_h_{:03d}".format(map_N)+corner+".dat"
h_exact = np.fromfile(ref_file, count=Nlon*Nlat, dtype=np.float32, offset=4*Nlon*Nlat*map_Nt).reshape((Nlat,Nlon))

print(np.min(h_exact), np.max(h_exact))

for scheme in schemes:
    h_file = "h_N{:03d}_dt".format(map_N)+dts[map_N]+"_"+scheme+corner+".dat"
    h = np.fromfile(h_file, count=Nlon*Nlat, dtype=np.float32, offset=4*Nlon*Nlat*map_Nt).reshape((Nlat,Nlon))
    err = h-h_exact

    cn_res.tiMainString = scheme
    plots.append(Ngl.contour_map(wks,err,cn_res))
    overlay_plot = Ngl.contour(wks, h_exact, overlay_res)
    Ngl.overlay(plots[-1], overlay_plot)
    add_cubed_sphere(wks,plots[-1])

pres = Ngl.Resources()
Ngl.panel(wks,plots,(2,2),pres)
Ngl.delete_wks(wks)

Ngl.end()
