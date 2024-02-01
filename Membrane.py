import pandas as pd
import math
import scipy
from helpers import find_closest_lines,normalize
from Trace import Trace

class Membrane:
    bounds = None
    roi=500
    regions={
        (None, "bsk"): "nucleus",
        ("bsk", "nucChnl"): "nuclear_basket",
        ("nucChnl", "midline"): "central_scaffold1",
        ("midline", "cytChnl"): "central_scaffold2",
        ("cytChnl", "npcCyto"): "cytoplasmic_fibril",
        ("npcCyto", None): "cytoplasm"
    }

    def __init__(self, df):
        self.bounds = self.readDF(df)


    def readDF(self, df):
        lineDict={
            "midline": None,
            "npcCyto": None,
            "bsk": None,
            "nucChnl": None,
            "cytChnl": None
            }
        lines=[
            df.iloc[:,0:2],
            df.iloc[:,2:4],
            df.iloc[:,4:6],
            df.iloc[:,6:8],
            df.iloc[:,8:10]]
        for c,key in enumerate(lineDict.keys()):
            tupList=[]
            for i,j in lines[c].iterrows():
                tupList.append((j[0],j[1]))
            lineDict[key]=tupList
        return(lineDict)
    
    def checkROI(self, ptX, ptY):
        first = self.bounds['midline'][0][1] - self.roi
        last = self.bounds['midline'][-1][1] + self.roi
        if ptY < first or ptY > last:
            return False
        for pt in self.bounds['midline']:
            pt = pt[0]
            if pt > 0:
                pt = pt+self.roi
                if ptX < pt:
                    return True
                else:
                    return False
            else:
                pt = pt-self.roi
                if ptX > pt:
                    return True
                else:
                    return False
                
    def traceParce(self,df,direction,pxl_size):
        df = df[['Frame','X','Y']].copy()
        normalize(df,pxl_size)
        df = df.sort_values(by='Frame')
        prevFrame = None
        traceList=[]
        tracePts=[]
        frames=[]
        for key,row in df.iterrows():
            if self.checkROI(row["X"],row["Y"]):
                if prevFrame != None and row['Frame'] != prevFrame+1:
                    if len(tracePts) > 2:
                        trace = Trace(tracePts,frames,direction)
                        traceList.append(trace)
                    tracePts=[]
                    frames=[]
                tracePts.append(((row['X'],row['Y'])))
                frames.append(row['Frame'])
                prevFrame = row['Frame']
        return traceList

    def locatePoints(self, traceList):
        for trace in traceList:
            regionList=[]
            for point in trace.trajectory:
                left_line, right_line = find_closest_lines(self.bounds, point)
                regionList.append(self.getRegion(left_line,right_line))
            trace.setRegions(regionList)
    
    def getRegion(self, left_line,right_line):
        return self.regions.get((left_line, right_line), None)
    
    def moveXbounds(self, x):
        for line in self.bounds.keys:
            if line == "bsk" or line == "npcCyto":
                for c, pt in enumerate(self.bounds[line]):
                    self.bounds[line][c][0] = pt[0]+x
            else:
                continue
                    

    
