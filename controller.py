from Membrane import Membrane
from statistics import mean
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
import os
import shutil
import argparse



def plotTrace(trace, out_path):
    plt.figure(figsize=(5,15))
    plt.xlim(-400,300)
    plt.title(f"Frames: {int(trace.frames[0])}-{int(trace.frames[-1])} | Class: {trace.classification}")
    plt.plot(*zip(*memb.bounds['midline']), marker='o', color="blue",markersize=1)
    plt.plot(*zip(*memb.bounds['npcCyto']), marker='o', color="orange",markersize=1)
    plt.plot(*zip(*memb.bounds['bsk']), marker='o', color="green",markersize=1)
    plt.plot(*zip(*memb.bounds['nucChnl']), marker='o', color="red",markersize=1)
    plt.plot(*zip(*memb.bounds['cytChnl']), marker='o', color="purple",markersize=1)

    plt.plot(*zip(*trace.trajectory), marker='o', color="pink",markersize=3)
    for c,x in enumerate(trace.trajectory):
        plt.text(x[0],x[1],str(c+1),fontsize=8)

    plt.legend(list(memb.bounds.keys()))

    plt.savefig(f"{out_path}/traces/{trace.classification}/Trace{int(trace.frames[0])}-{int(trace.frames[-1])}.png")
    if trace.crossedMidline():
            plt.savefig(f"{out_path}/traces/midline/Trace{int(trace.frames[0])}-{int(trace.frames[-1])}.png")
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("membrane", help="Path to membrane csv file")
parser.add_argument("trace", help="Path to traces csv file")
parser.add_argument("pxl_size", help="PXL_SIZE used for membrane/trace data", type=int)
parser.add_argument("--plot","-p", help="Plot traces", action="store_true")
parser.add_argument("-i", help="Use this flag if all data is import", action="store_true")
parser.add_argument("-e", help="Use this flag if all data is export", action="store_true")
parser.add_argument("--normalize", '-n', help="Normalize trace points by given pixel size.", action="store_true")
parser.add_argument("--plot_incomplete","-pi", help="Plots incomplete traces", action="store_true")

args = parser.parse_args()

PXL_SIZE=args.pxl_size
NORMALIZE=args.normalize
IMPORT=args.i
EXPORT=args.e
DIRECTION=None

if IMPORT:
    DIRECTION=1
elif EXPORT:
    DIRECTION=2
else:
    DIRECTION=3

start = time.time()

df = pd.read_csv(args.membrane)
df2 = pd.read_csv(args.trace) 

completeDict={}
incompleteDict={}


memb = Membrane(df)

traceList = memb.traceParce(df2,DIRECTION,PXL_SIZE,NORMALIZE)
# for c,trace in enumerate(traceList):
#     print(f"{c}: {int(trace.frames[0])}-{int(trace.frames[-1])}")
# memb.locatePoints(traceList)
# print(traceList[62].regions)
# exit()

completeTraces = []
incompleteTraces = []

for c, trace in enumerate(traceList):
    if not trace.prune():
        completeTraces.append(trace)
        trace.classifyComplete()
        trace.setAvgDist()
        trace.regionDistance()
        trace.setDwellTime()
    else:
        incompleteTraces.append(trace)
        # print(f"Incomplete: {trace.regions}")


print(f"Out of {len(traceList)} traces there were: \n-{len(completeTraces)} complete traces\n-{len(incompleteTraces)} incomplete traces")
print(f"Total Time Spent: {round((time.time()-start),2)} seconds")

if DIRECTION == 1:
    out_path=f"./output/{args.trace.split('/')[1].split('.')[0]}_import"
elif DIRECTION == 2:
    out_path=f"./output/{args.trace.split('/')[1].split('.')[0]}_export"

if os.path.exists(f"{out_path}/traces"):
    shutil.rmtree(f"{out_path}/traces")
else:
    os.mkdir(f"{out_path}")
os.mkdir(f"{out_path}/traces")
os.mkdir(f"{out_path}/traces/successful")
os.mkdir(f"{out_path}/traces/docking_event")
os.mkdir(f"{out_path}/traces/abortive_central_channel")
os.mkdir(f"{out_path}/traces/abortive_cytoplasmic_fibril")
os.mkdir(f"{out_path}/traces/abortive_nuclear_basket")
os.mkdir(f"{out_path}/traces/midline")
os.mkdir(f"{out_path}/traces/incomplete")

regDistsTotal = {"nucleus":[],"nuclear_basket":[],"central_scaffold1":[],"central_scaffold2":[],"cytoplasmic_fibril":[],"cytoplasm":[]}
output=[]
for c1, trace in enumerate(completeTraces):
    regDistsTrace = {"nucleus":[],"nuclear_basket":[],"central_scaffold1":[],"central_scaffold2":[],"cytoplasmic_fibril":[],"cytoplasm":[]}
    print("------------------------------------------------------------------")
    print(f"Complete #{c1}: {trace.regions}")
    print(f"Frames: {int(trace.frames[0])}-{int(trace.frames[-1])} \nClass: {trace.classification}\nDwell Time: {trace.dwellTime} frame(s)")
    if len(trace.sameRegionDists) > 0 :
        print(f'Average Distances:{trace.sameRegionDists}')
        for dist in trace.sameRegionDists:
            regDistsTotal[dist[0]].append(dist[1])
            regDistsTrace[dist[0]].append(dist[1])

    nuc=round(mean(regDistsTrace['nucleus']),2) if len(regDistsTrace['nucleus']) > 0 else None
    nucBask=round(mean(regDistsTrace['nuclear_basket']),2) if len(regDistsTrace['nuclear_basket']) > 0 else None
    cs1=round(mean(regDistsTrace['central_scaffold1']),2) if len(regDistsTrace['central_scaffold1']) > 0 else None
    cs2=round(mean(regDistsTrace['central_scaffold2']),2) if len(regDistsTrace['central_scaffold2']) > 0 else None
    cf=round(mean(regDistsTrace['cytoplasmic_fibril']),2) if len(regDistsTrace['cytoplasmic_fibril'])> 0  else None
    cyto=round(mean(regDistsTrace['cytoplasm']),2) if len(regDistsTrace['cytoplasm']) > 0 else None
    
    outDict={'StartFrame':int(trace.frames[0]),'EndFrame':int(trace.frames[-1]),
                    'Class':trace.classification,'DwellTime':trace.dwellTime,
                    'Nucleus':nuc,
                    'NuclearBasket':nucBask,
                    'CentralScaffold1':cs1,
                    'CentralScaffold2':cs2,
                    'CytoplasmicFibril':cf,
                    'Cytoplasm':cyto
                }
    output.append(outDict)
    
    if args.plot:
        if c1 == 0:
            print("Plotting complete traces...")
        plotTrace(trace,out_path)
       
print("------------------------------------------------------------------\n")
#--------------------------PLOT ALL TRACES----------------------------------
plt.figure(figsize=(5,15))
for c1, trace in enumerate(completeTraces):
    plt.xlim(-400,300)
    if c1 == 0:
        plt.title(f"All Traces")
        plt.plot(*zip(*memb.bounds['midline']), marker='o', color="blue",markersize=1)
        plt.plot(*zip(*memb.bounds['npcCyto']), marker='o', color="orange",markersize=1)
        plt.plot(*zip(*memb.bounds['bsk']), marker='o', color="green",markersize=1)
        plt.plot(*zip(*memb.bounds['nucChnl']), marker='o', color="red",markersize=1)
        plt.plot(*zip(*memb.bounds['cytChnl']), marker='o', color="purple",markersize=1)
    plt.plot(*zip(*trace.trajectory), marker='o', color="pink",markersize=3)
plt.legend(list(memb.bounds.keys()))
plt.savefig(f"{out_path}/traces/AllTraces.png") 
plt.close()

if args.plot_incomplete:
    for c1, trace in enumerate(incompleteTraces):
        if len(trace.regions) <= 2:
            continue
        plotTrace(trace,out_path)

#-------------------------------------AVERAGE DISTS--------------------------
nucTot=round(mean(regDistsTotal['nucleus']),2) if len(regDistsTotal['nucleus']) > 0 else None
nucBaskTot=round(mean(regDistsTotal['nuclear_basket']),2) if len(regDistsTotal['nuclear_basket']) > 0 else None
cs1Tot=round(mean(regDistsTotal['central_scaffold1']),2) if len(regDistsTotal['central_scaffold1']) > 0 else None
cs2Tot=round(mean(regDistsTotal['central_scaffold2']),2) if len(regDistsTotal['central_scaffold2']) > 0 else None
cfTot=round(mean(regDistsTotal['cytoplasmic_fibril']),2) if len(regDistsTotal['cytoplasmic_fibril'])> 0  else None
cytoTot=round(mean(regDistsTotal['cytoplasm']),2) if len(regDistsTotal['cytoplasm']) > 0 else None

print(f"Regional Average Distances: Nucleus: {nucTot}nm, " +
        f"NuclearBasket: {nucBaskTot}nm, " +
        f"CentralScaffold1: {cs1Tot}nm, " +
        f"CentralScaffold2: {cs2Tot}nm, " +
        f"CytoplasmicFibril: {cfTot}nm, " +
        f"Cytoplasm: {cytoTot}nm")

pd.DataFrame(output).to_csv(f"{out_path}/out.csv")

