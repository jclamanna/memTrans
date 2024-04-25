import pandas as pd
from Trace import Trace, distance
import numpy as np

def normalize(df, PXL_SIZE):
    df["X"] = list(map(lambda x: x*PXL_SIZE, df.loc[:,"X"]))
    df["Y"] = list(map(lambda x: x*PXL_SIZE, df.loc[:,"Y"]))
    xAvg = df["X"].mean()
    yAvg = df["Y"].mean()
    df["X"] = list(map(lambda x: x-xAvg, df.loc[:,"X"]))
    df["Y"] = list(map(lambda x: x-yAvg, df.loc[:,"Y"]))
    return df

def closest_point_on_segment(segment, point):
    p1, p2 = segment
    p1, p2 = np.array(p1), np.array(p2)
    v = np.array(p2) - np.array(p1)
    u = np.array(point) - np.array(p1)
    t = np.dot(u, v) / np.dot(v, v)
    if t < 0:
        closest_point = p1
    elif t > 1:
        closest_point = p2
    else:
        closest_point = p1 + t * (p2 - p1)
    return closest_point

def find_closest_lines(lines, point):
    line_distances = {}
    
    for line_name, line_points in lines.items():
        min_dist = float('inf')
        closest_point_on_line = None
        for i in range(len(line_points) - 1):
            segment = (line_points[i], line_points[i + 1])
            closest_point = closest_point_on_segment(segment, point)
            dist = distance(point, closest_point)
            if dist < min_dist:
                min_dist = dist
                closest_point_on_line = closest_point
        line_distances[line_name] = (min_dist, closest_point_on_line[0])

    sorted_lines = sorted(line_distances.keys(), key=lambda line: line_distances[line][1])

    closest_left = None
    closest_right = None
    for i, line in enumerate(sorted_lines):
        if line_distances[line][1] < point[0]:
            closest_left = line
        elif line_distances[line][1] > point[0]:
            closest_right = line
            break

    return closest_left, closest_right

