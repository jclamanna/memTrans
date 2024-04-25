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

    def __init__(self, ptList,frames):
        self.trajectory = ptList
        self.frames = frames
    
    def __str__(self):
        return("".join(str(self.trajectory)))
    
    def setRegions(self, regions):
        self.regions = regions
    
    def setDirection(self, dir):
        self.direction = dir

    def prune(self):
        prev = ""
        start_index = None
        end_index = None

        if self.direction==1:
            origin="cytoplasm"
            target="nucleus"
        elif self.direction==2:
            origin="nucleus"
            target="cytoplasm"

        if origin not in self.regions:
            self.classification="incomplete"
            return True

        for c, r in enumerate(self.regions):
            # Handle the first element in the list
            if c == 0:
                if r == origin:
                    if self.regions[1] == target:
                        self.classification="incomplete"
                        return True
                    elif self.regions[1] != origin:
                        start_index = 0
                        continue
                else: 
                    if self.regions[1] == origin:
                        start_index = 1
                        continue
                    else:
                        self.classification="incomplete"
                        return True
 
            # Check if the previous region was the target (nucleus or cytoplasm)
            if prev == origin:
                if r == target:
                    self.classification="incomplete"
                    return True
                elif r != target and start_index is None:
                    start_index = c - 1
                    continue

            if start_index is not None:
                # Determine the end_index based on the current region
                if r == target:
                    if c + 1 < len(self.regions) and self.regions[c + 1] == target:
                        end_index = c + 1
                        break
                    else:
                        end_index = c
                        break
                elif r == origin and prev != origin:
                    end_index = c
            prev = r 
            
        if start_index is not None and end_index is not None:
            self.regions = self.regions[start_index:end_index + 1]
            self.trajectory = self.trajectory[start_index:end_index + 1]
            self.frames = self.frames[start_index:end_index + 1]

        if (self.regions[0] != origin) or (self.regions[-1] != origin and self.regions[-1] != target) or (len(self.regions) <= 1) or (len(set(self.regions))==1):
            self.classification="incomplete"
            return True
            
        return False
    
    def isSuccess(self):
        target = "nucleus" if self.direction == 1 else "cytoplasm"
        if target in self.regions:
            return True
        return False
    
    def crossedMidline(self):
        if self.direction == 1:
            midZones = ["central_scaffold1","nuclear_basket","nucleus"]
        else:
            midZones = ["central_scaffold2","cytoplasmic_fibril","cytoplasm"]
        for zone in midZones:
            if zone in self.regions:
                return True
        return False
    
    def deepestRegion(self):
        order = ["nucleus","nuclear_basket","central_scaffold1","central_scaffold2","cytoplasmic_fibril","cytoplasm"] if self.direction != 1 else ["cytoplasm","cytoplasmic_fibril","central_scaffold2","central_scaffold1","nuclear_basket","nucleus"]
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
        # print(deep)
        if self.direction != 1:
            if deep == "nuclear_basket":
                self.classification="docking_event"
                return None
            elif deep == "central_scaffold1" or deep == "central_scaffold2":
                self.classification="abortive_central_channel"
                return None
            elif deep == "cytoplasmic_fibril":
                self.classification = "abortive_cytoplasmic_fibril"
                return None
        else:
            if deep == "cytoplasmic_fibril":
                self.classification="docking_event"
                return None
            elif deep == "central_scaffold1" or deep == "central_scaffold2":
                self.classification="abortive_central_channel"
                return None
            elif deep == "nuclear_basket":
                self.classification = "abortive_nuclear_basket"
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