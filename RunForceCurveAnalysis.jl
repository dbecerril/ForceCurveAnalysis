push!(LOAD_PATH, pwd())

using Revise
using DelimitedFiles
using Plots
#pyplot()
plotly()

import ForceCurveAnalysis
fca = ForceCurveAnalysis

########################################################
#### Run this section to define the system to analyse
#########################################################

#folder_path = string(pwd(),"/testData" )
#folder_path =  "/media/david/389D-7C2D/Dic2021/prova1" #blackRed usb
folder_path = "/media/david/84E6-AD45/Jan2022/panc1"    #little usb

#folder_path = "/media/david/Daves/Pancreas/Pac1" #purple
file_list   = readdir(folder_path);

# read file parameters
rfp = [folder_path, "Jan26_PANC1",51, "2"]

#elastic cte.
elasticConst = 0.0002 #0.02

#detector sensitivity (V/nm)
Dsens = 0.000356

#Zpiezo sensitivity nm/step
Zsens = 0.084

# Outlier Tolarance 
tol1 = 0.07

# points in moving average window
nwin = 2

# Variables used to identify zero point part of curve
#[percentageInit,percentageFinl] 
zps = [1.0, 60.0]

# Variables used to make the derivative curve 
# [iterations,weigth,dx] ;     default: 200, 0.2, 1
dp = [150.0, 0.2, 2.0]

#xi,xf,yi,yf
mapsec = [1, rfp[3], 1, rfp[3]]

#indentatin curve part
icp = [20, 70]

system1 = fca.SystemParameters(rfp, elasticConst, Dsens, Zsens, tol1, nwin, zps, dp, mapsec, icp)
#####################################################################
###################### end of system definition #######################


######################################################
#################### OPtimizer  #######################
#######################################################
#@doc forceCurveOptimizer

#list0        =["glass1_Andata_14.txt","glass2_Andata_15.txt","cellfc2_Andata_17.txt","cellfc1_Andata_16.txt" ]
#list0 = readdir(folder_path)

y0 = 45
x0 = 4
points = [x0 12;
    x0 22;
    x0 32;
    x0 42];
    
a      = fca.forceCurveOptimizer( points,2, system1);
plot(a[1],a[2],a[3],a[4] ,layout = (2,2), fmt = :png)

####################  end optimizer#######################
#######################################################


#####################0.0007#.0007##################################
#################### Run  #######################
#######################################################


xi = 1
yi = 1


xf = rfp[3]
yf = 51

system1.AnalyzeSection = [xi,xf,yi,yf] 
@time slopeMap = fca.analyzeForceCurve(system1)

# write to files 
fca.writeToFiles(slopeMap,"jan26_PAC1",system1)


####################  end Run #######################
#######################################################