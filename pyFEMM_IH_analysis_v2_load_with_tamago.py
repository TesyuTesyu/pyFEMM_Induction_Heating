import matplotlib.pyplot as plt
import numpy as np
import csv

#たまごさんの式.
#パワーエレクトロニクス第4号, タマゴさん, 部分共振ZCS方式電流形並列共振プッシュプルインバータによる誘導炉の製作, 2024.
#https://pawami-sha.tumblr.com/
def tamago_calcu(rw,l,N,sigma,delta):
    R=np.pi*N**2/(2*l*sigma*delta*(np.exp(rw/delta)+np.exp(-rw/delta))**2)*(np.exp(2*rw/delta)*(2*rw-delta)-np.exp(-2*rw/delta)*(2*rw+delta)+2*delta-4*rw**2/delta)
    return R

#pyFEMM_IH_analysis_v2.py が書き出したファイルたちを読み出す.

folder_path=r"C:\Users\???"#フォルダパス.
matfile_path=["pyFEMM_IH_v1_20deg.csv","pyFEMM_IH_v1_20deg_half_radius_work.csv"]#ファイル名.いくつも可能.
Rdc_file_path="pyFEMM_IH_v1_wo_work.csv"#ワークなしにおける結果.
RDC=0.00196269#DC抵抗.

file_path=folder_path+"\\"+Rdc_file_path
matRcoil=[]
with open(file_path, 'r', encoding='utf-8') as file:
    csvreader = csv.reader(file, delimiter=',')
    i=0
    for row in csvreader:
        if i==0:
            values = [float(x) for x in row]
            [sigma_coil,sigma_work0,Temp0,alpha_work,mu,I_coil,scale,r_wire,mesh_size_coil,r_coil,r_work,h_coil,y_work_bottom,y_work_top,turn,Temp,freq_min,freq_max,freq_num]=values
        elif i==1:
            print(row)
            break
        i+=1
    for row in csvreader:
        values = [float(x) for x in row]
        matRcoil.append(values[0]*2)

matRcoil=np.array(matRcoil)

fig3, ax3 = plt.subplots(nrows=1, ncols=1, squeeze=False, tight_layout=True, figsize=[8,6], sharex = "col")

mat_eta=[]
matR=[]
matX=[]
for j in range(len(matfile_path)):
    file_path=folder_path+"\\"+matfile_path[j]
    mat_eta=[]
    matR=[]
    matX=[]
    with open(file_path, 'r', encoding='utf-8') as file:
        csvreader = csv.reader(file, delimiter=',')
        i=0
        for row in csvreader:
            if i==0:
                values = [float(x) for x in row]
                [sigma_coil,sigma_work0,Temp0,alpha_work,mu,I_coil,scale,r_wire,mesh_size_coil,r_coil,r_work,h_coil,y_work_bottom,y_work_top,turn,Temp,freq_min,freq_max,freq_num]=values
            elif i==1:
                print(row)
                break
            i+=1
        for row in csvreader:
            values = [float(x) for x in row]
            mat_eta.append(values[1]/values[0])
            matR.append(values[0]*2)
            matX.append(values[2])

    matf=np.logspace(np.log10(freq_min),np.log10(freq_max),int(freq_num))
    sigma_work=sigma_work0/(1+alpha_work*(Temp-Temp0))
    R_tamago=[]
    for f in matf:
        delta=np.sqrt(2/(2*np.pi*f*sigma_work*1e6*mu))
        R_tamago.append(tamago_calcu(r_work*scale,y_work_top*scale-y_work_bottom*scale,turn,sigma_work*1e6,delta))
    
    ax3[0,0].plot(matf,np.array(matR)-np.array(matRcoil),"k-")
    ax3[0,0].plot(matf,np.array(R_tamago),"r--")
    ax3[0,0].set_xscale('log')
    ax3[0,0].set_yscale('log')

plt.show()
