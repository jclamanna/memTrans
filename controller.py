from Membrane import Membrane
from statistics import mean
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import os
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("membrane", help="Path to membrane csv file")
parser.add_argument("trace", help="Path to traces csv file")
parser.add_argument("pxl_size", help="PXL_SIZE used for membrane/trace data", type=int)
parser.add_argument("--plot","-p", help="Plot traces", action="store_true")
args = parser.parse_args()

PXL_SIZE=args.pxl_size
DIRECTION=1

start = time.time()

df = pd.read_csv(args.membrane)
df2 = pd.read_csv(args.trace) 

# xMoves=[-30,-5,0,5,30]
xMoves=[0]
completeDict={}
incompleteDict={}

for x in xMoves:
    memb = Membrane(df)
    traceList = memb.traceParce(df2,DIRECTION,PXL_SIZE)
    memb.locatePoints(traceList)

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
    completeDict[x]= completeTraces
    incompleteDict[x] = incompleteTraces

maxX = 0
for key, val in completeDict.items():
    print (len(val))
    if len(val) > maxX:
        maxX = key

print(f"The membrane that resulted in the most complete traces was {maxX}")
completeTraces = completeDict[maxX]
incompleteTraces = incompleteDict[maxX]
print(f"Out of {len(traceList)} traces there were: \n-{len(completeTraces)} complete traces\n-{len(incompleteTraces)} incomplete traces")
print(f"Total Time Spent: {round((time.time()-start),2)} seconds")

out_path=f"./output/{args.trace.split('/')[1].split('.')[0]}"
if os.path.exists(f"{out_path}/traces"):
    shutil.rmtree(f"{out_path}/traces")
else:
    os.mkdir(f"{out_path}")
os.mkdir(f"{out_path}/traces")
os.mkdir(f"{out_path}/traces/successful")
os.mkdir(f"{out_path}/traces/docking_event")
os.mkdir(f"{out_path}/traces/failure_central_channel")
os.mkdir(f"{out_path}/traces/failure_cytoplasmic_fibril")
os.mkdir(f"{out_path}/traces/midline")

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
        plt.figure(figsize=(5,15))
        plt.xlim(-400,300)
        plt.title(f"Frames: {int(trace.frames[0])}-{int(trace.frames[-1])} | Class: {trace.classification}")
        plt.plot(*zip(*memb.bounds['midline']), marker='o')
        plt.plot(*zip(*memb.bounds['npcCyto']), marker='o')
        plt.plot(*zip(*memb.bounds['bsk']), marker='o')
        plt.plot(*zip(*memb.bounds['nucChnl']), marker='o')
        plt.plot(*zip(*memb.bounds['cytChnl']), marker='o')
        plt.plot(*zip(*trace.trajectory), marker='.', color="pink")
        for c,x in enumerate(trace.trajectory):
            plt.text(x[0],x[1],str(c+1),fontsize=8)
        plt.legend(list(memb.bounds.keys()))
        plt.savefig(f"{out_path}/traces/{trace.classification}/Trace{int(trace.frames[0])}-{int(trace.frames[-1])}.png")
        if trace.crossedMidline():
                plt.savefig(f"{out_path}/traces/midline/Trace{int(trace.frames[0])}-{int(trace.frames[-1])}.png")
        plt.close()

print("------------------------------------------------------------------\n")

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

