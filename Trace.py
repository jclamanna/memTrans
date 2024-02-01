from statistics import mean
import numpy as np

def distance(point1, point2):
    return np.linalg.norm(np.array(point1) - np.array(point2))

class Trace:
    trajectory=None
    frames = None
    direction = None
    classification = None
    avgDist = None
    dwellTime=None
    sameRegionDists = []
    regions = []

    def __init__(self, ptList,frames,direction):
        self.trajectory = ptList
        self.frames = frames
        self.direction = direction
    
    def __str__(self):
        return("".join(str(self.trajectory)))
    
    def setRegions(self, regions):
        self.regions = regions

    def prune(self):
        prev = ""
        start_index = None
        end_index = None
        target = "nucleus" if self.direction == 1 else "cytoplasm"
        opposite = "cytoplasm" if self.direction == 1 else "nucleus"

        if target not in self.regions:
            self.classification="incomplete"
            return True
        if len(set(self.regions)) == 1:
            self.classification="incomplete"
            return True

        for c, r in enumerate(self.regions):
            # Handle the first element in the list
            if c == 0 and r == target:
                if self.regions[1] == opposite:
                    self.classification="incomplete"
                    return True
                elif self.regions[1] != target:
                    start_index = 0
                    continue
            # Check if the previous region was the target (nucleus or cytoplasm)
            if prev == target:
                if r == opposite:
                    self.classification="incomplete"
                    return True
                elif r != target and start_index is None:
                    start_index = c - 1
                    continue
            if start_index is not None:
                # Determine the end_index based on the current region
                if r == opposite:
                    if c + 1 < len(self.regions) and self.regions[c + 1] == opposite:
                        end_index = c + 1
                        break
                    else:
                        end_index = c
                        break
                elif r == target and prev != target:
                    end_index = c
            prev = r 
            
        if start_index is not None and end_index is not None:
            self.regions = self.regions[start_index:end_index + 1]
            self.trajectory = self.trajectory[start_index:end_index + 1]
            self.frames = self.frames[start_index:end_index + 1]
        if self.regions[0] != target or (self.regions[-1] != target and self.regions[-1] != opposite):
            self.classification="incomplete"
            return True
            
        return False
    
    def isSuccess(self):
        target = "nucleus" if self.direction != 1 else "cytoplasm"
        if target in self.regions:
            return True
        return False
    
    def crossedMidline(self):
        if self.direction == 1:
            midZones = ["central_scaffold2","cytoplasmic_fibril","cytoplasm"]
        else:
            midZones = ["central_scaffold1","nuclear_basket","nucleus"]
        for zone in midZones:
            if zone in self.regions:
                return True
        return False
    
    def deepestRegion(self):
        order = ["nucleus","nuclear_basket","central_scaffold1","central_scaffold2","cytoplasmic_fibril","cytoplasm"] if self.direction == 1 else ["cytoplasm","cytoplasmic_fibril","central_scaffold2","central_scaffold1","nuclear_basket","nucleus"]
        regionIndex =[]
        for r in self.regions:
            regionIndex.append(order.index(r))
        return order[max(regionIndex)]

    def classifyComplete(self):
        if self.classification == "incomplete":
            print("Cannot classify incomplete trace!")
            return None
        
        if self.isSuccess():
            self.classification = "successful"
            return None
        
        deep = self.deepestRegion()
        if self.direction == 1:
            if deep == "nuclear_basket":
                self.classification="docking_event"
                return None
            elif deep == "central_scaffold1" or deep == "central_scaffold2":
                self.classification="failure_central_channel"
                return None
            elif deep == "cytoplasmic_fibril":
                self.classification = "failure_cytoplasmic_fibril"
                return None
    
    def setAvgDist(self):
        dists = set()
        pts = self.trajectory
        for c,pt in enumerate(pts):
            if c < len(pts)-1:
                dists.add(distance(pt,pts[c+1]))
        self.avgDist = mean(dists)
    
    def regionDistance(self):
        distList = []
        pts = self.trajectory
        for c,r in enumerate(self.regions):
            rDist = None
            if c < len(self.regions)-1:
                if r == self.regions[c+1]:
                    rDist = (r, round(distance(pts[c],pts[c+1]),2))
                if rDist != None:
                    # print(rDist)
                    distList.append(rDist)
        self.sameRegionDists = distList

    def setDwellTime(self):
        dt = 0
        for r in self.regions:
            if r != "nucleus" and r != "cytoplasm":
                dt+=1
        self.dwellTime = dt