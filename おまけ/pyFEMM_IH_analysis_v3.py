import femm
import matplotlib.pyplot as plt
import numpy as np
import csv

#ワークの外側の表面のみ、メッシュを細かく刻むようにしたもの.
#導電率, 温度係数: http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/rstiv.html

#---------------------config----------------------

sigma_coil=58#MS/m , 銅の導電率.
sigma_work0=37.7#MS/m , 加熱対象の導電率. @20 degree.
Temp0=20#基準温度. 20 degree.
alpha_work=4.29e-3#1/degree, 加熱対象の温度係数.
mu=4*np.pi*1e-7#透磁率.
I_coil=1#コイルの電流.変更不可.

scale=1e-3#長さのスケール.

r_wire=2.5#wire radius.
mesh_size_min=0.1
mesh_size_max=1
mesh_size_ratio=2 #mesh_size = 表皮厚さ / mesh_size_ratio
mesh_size_large=1 #[mm] 加熱対象内部の荒いメッシュサイズ.
k=2#r_work_mid=r_work-delta*k

r_coil=42#coil radius.
h_coil=75#coil hight.
r_work=28#work radius.

y_work_bottom=0#加熱対象の底面のy座標.
y_work_top=h_coil#加熱対象の上底のy座標.

turn=8#コイルの巻き数.

Temp=20#温度. degree.
freq_min=1e2#周波数の下限.
freq_max=1e4#周波数の上限.
freq_num=2#周波数の分割数.

#結果ファイルの保存path.
file_path_csv_write=r"C:\Users\???\jhjh.csv"

#--------------------------------------------------

matf=np.logspace(np.log10(freq_min),np.log10(freq_max),freq_num)

sigma_work=sigma_work0/(1+alpha_work*(Temp-Temp0))

x_air=(r_coil+r_work)/2
y_air=h_coil/2

y_work=(y_work_top+y_work_bottom)/2

#https://www.femm.info/wiki/pyfemm

#周波数を振って何度も実行.
matans=[]
for freq in matf:
    # The package must be initialized with the openfemm command.
    femm.openfemm()
    # We need to create a new Magnetostatics document to work on.
    femm.newdocument(0)

    femm.mi_probdef(freq, 'millimeters', 'axi', 1.e-8, 0, 30)

    delta=np.sqrt(2/(2*np.pi*freq*sigma_work*1e6*mu))*1e3
    mesh_size=delta/mesh_size_ratio
    if mesh_size>mesh_size_max:
        mesh_size=mesh_size_max
    elif mesh_size<mesh_size_min:
        mesh_size=mesh_size_min

    r_work_mid=r_work-delta*k
    if r_work_mid<0:
        r_work_mid=mesh_size
    
    #コイルの描画.
    mat_y_coil=np.linspace(0,h_coil,turn)
    for y in mat_y_coil:
        femm.mi_drawarc(r_coil+r_wire,y-r_wire,r_coil+r_wire,y+r_wire,180,1)#x1,y1,x2,y2,angle,maxseg
        femm.mi_drawarc(r_coil+r_wire,y+r_wire,r_coil+r_wire,y-r_wire,180,1)
        femm.mi_addblocklabel(r_coil+r_wire,y)

    #加熱対象の描画.
    femm.mi_drawrectangle(0, y_work_bottom, r_work, y_work_top)
    femm.mi_drawline(r_work_mid,y_work_bottom,r_work_mid,y_work_top)
    femm.mi_addblocklabel(r_work_mid/2,y_work)
    femm.mi_addblocklabel((r_work_mid+r_work)/2,y_work)

    #コイル内の材質とCircuitを定義.
    for y in mat_y_coil:
        femm.mi_selectlabel(r_coil+r_wire,y)
        femm.mi_setblockprop('Coil', 0, mesh_size, 'icoil', 0, 0, 1)
        femm.mi_clearselected()
    
    femm.mi_selectlabel(r_work_mid/2,y_work)
    femm.mi_setblockprop('Work1', 0, mesh_size_large, '<None>', 0, 0, 0)
    femm.mi_clearselected()

    femm.mi_selectlabel((r_work_mid+r_work)/2,y_work)
    femm.mi_setblockprop('Work2', 0, mesh_size, '<None>', 0, 0, 0)
    femm.mi_clearselected()

    # Define an "open" boundary condition using the built-in function:
    femm.mi_makeABC()

    # Add block labels, one to each the steel, coil, and air regions.
    femm.mi_addblocklabel(x_air,y_air)

    # Add some block labels materials properties
    femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
    femm.mi_addmaterial('Coil', 1, 1, 0, 0, sigma_coil, 0, 0, 1, 0, 0, 0)
    femm.mi_addmaterial('Work1', 1, 1, 0, 0, sigma_work, 0, 0, 1, 0, 0, 0)
    femm.mi_addmaterial('Work2', 1, 1, 0, 0, sigma_work, 0, 0, 1, 0, 0, 0)

    # Add a "circuit property" so that we can calculate the properties of the
    femm.mi_addcircprop('icoil', I_coil, 1)

    # Apply the materials to the appropriate block labels
    femm.mi_selectlabel(x_air,y_air)
    femm.mi_setblockprop('Air', 0, 1, '<None>', 0, 0, 0)
    femm.mi_clearselected()

    # Now, the finished input geometry can be displayed.
    femm.mi_zoomnatural()

    # We have to give the geometry a name before we can analyze it.
    femm.mi_saveas('IH.fem')

    # Now,analyze the problem and load the solution when the analysis is finished
    femm.mi_analyze()
    femm.mi_loadsolution()

    vals=femm.mo_getcircuitproperties('icoil')
    current=vals[0]
    voltage=vals[1]
    P_coil=np.real(voltage)/2

    femm.mo_seteditmode('area')
    femm.mo_selectblock(r_work_mid/2,y_work)
    P_work1=femm.mo_blockintegral(4)
    P_work1=np.real(P_work1)
    
    femm.mo_selectblock((r_work_mid+r_work)/2,y_work)
    P_work2=femm.mo_blockintegral(4)
    P_work2=np.real(P_work2)

    P_work=P_work1+P_work2

    print(P_coil, P_work,np.imag(voltage))
    matans.append([P_coil, P_work,np.imag(voltage)])

    femm.mo_close()

    # When the analysis is completed, FEMM can be shut down.
    femm.closefemm()

mat_eta=[]
matR=[]
matX=[]
for y in matans:
    P_coil=y[0]
    P_work=y[1]
    mat_eta.append(P_work/P_coil)
    matR.append(P_coil*2)
    matX.append(y[2])

fig, ax = plt.subplots(nrows=1, ncols=1, squeeze=False, tight_layout=True, figsize=[8,6], sharex = "col")
ax[0,0].plot(matf,mat_eta,"k-")
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')

fig, ax = plt.subplots(nrows=1, ncols=1, squeeze=False, tight_layout=True, figsize=[8,6], sharex = "col")
ax[0,0].plot(matf,matR,"k-")
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')

matX=np.array(matX)
fig, ax = plt.subplots(nrows=1, ncols=1, squeeze=False, tight_layout=True, figsize=[8,6], sharex = "col")
ax[0,0].plot(matf,matX/(2*np.pi*matf),"k-")
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')
plt.show()

#csvに書き出し.
with open(file_path_csv_write, 'w', newline='') as fil:
    writer = csv.writer(fil)
    writer.writerow([sigma_coil,sigma_work0,Temp0,alpha_work,mu,I_coil,scale,r_wire,mesh_size,r_coil,r_work,h_coil,y_work_bottom,y_work_top,turn,Temp,freq_min,freq_max,freq_num])
    writer.writerow(["P_coil, P_work, X_coil"])
    for y in matans:
        writer.writerow(y)
