'''
Author: Neeraj Sirdeshmukh
Delft University Of Technology, Netherlands 
'''
import math
from math import * 
from pyproj import *
from laspy.file import File
import numpy as np
import time
import pandas
import os
import psycopg2
from psycopg2.extras import *
from io import StringIO


class isea_pt(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
         
class isea_geo(object):
    # in radians
    def __init__(self, lon, lat):
        self.lon = lon 
        self.lat = lat 

# Set constants

M_PI = 3.14159265358979323846


E = 52.62263186
F = 10.81231696

DEG60 =1.04719755119659774614
DEG120= 2.09439510239319549229
DEG72 =1.25663706143591729537
DEG90 =1.57079632679489661922
DEG144= 2.51327412287183459075
DEG36 =0.62831853071795864768
DEG108 =1.88495559215387594306
DEG180= M_PI

ISEA_SCALE=0.8301572857837594396028083

# Degrees: 26.565051177065900134
V_LAT = 0.46364760899944494524

RAD2DEG = 180.0/M_PI
DEG2RAD = M_PI/180.0

# Icosahedron constants, radians

g =37.37736814 * DEG2RAD
G =36.0 * DEG2RAD
theta =30.0 * DEG2RAD

# vertices of triangles
vertex = [isea_geo(0.0,DEG90),
        isea_geo(DEG180,V_LAT),
        isea_geo(-DEG108,V_LAT),
        isea_geo(-DEG36,V_LAT),
        isea_geo(DEG36,V_LAT),
        isea_geo(DEG108,V_LAT),
        isea_geo(-DEG144,-V_LAT),
        isea_geo(-DEG72,-V_LAT),
        isea_geo(0.0,-V_LAT),
        isea_geo(DEG72,-V_LAT),
        isea_geo(DEG144,-V_LAT),
        isea_geo(0.0,-DEG90)
        ]

E_RAD =0.91843818702186776133
F_RAD =0.18871053072122403508


# icosahedron triangle centers 
icostriangles = [
    isea_geo(0.0, 0.0),
    isea_geo(-DEG144, E_RAD),
    isea_geo(-DEG72, E_RAD),
    isea_geo(0.0, E_RAD),
    isea_geo(DEG72, E_RAD),
    isea_geo(DEG144, E_RAD),
    isea_geo(-DEG144, F_RAD),
    isea_geo(-DEG72, F_RAD),
    isea_geo(0.0, F_RAD),
    isea_geo(DEG72, F_RAD),
    isea_geo(DEG144, F_RAD),
    isea_geo(-DEG108, -F_RAD),
    isea_geo(-DEG36, -F_RAD),
    isea_geo(DEG36, -F_RAD),
    isea_geo(DEG108, -F_RAD),
    isea_geo(DEG180, -F_RAD),
    isea_geo(-DEG108, -E_RAD),
    isea_geo(-DEG36, -E_RAD),
    isea_geo(DEG36, -E_RAD),
    isea_geo(DEG108, -E_RAD),
    isea_geo(DEG180, -E_RAD)
    
    ]

# Parameters taken from Snyder (1992)
TABLE_G = 0.6615845383
TABLE_H = 0.1909830056

RPRIME = 0.91038328153090290025


tri_v1 = [0, 0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 2, 3, 4, 5, 1, 11, 11, 11, 11, 11]

    
NEWORIGX = -0.6022955012659694 # TABLE_G * (-1) #old: 
NEWORIGY = -0.3477354703761901 # TABLE_H * (-2) #old: 

TANTHETA = tan(theta)
COTTHETA = 1.0 / TANTHETA

SINUPPERG = sin(G)
COSUPPERG = cos(G)

COSLOWERG = cos(g)
TANLOWERG = tan(g) 

# parameters for WGS84 ellipsoid

R = 6378137 # semi-major axis, meters 
b = 6356752.314245 # semi-minor axis, meters
flattening =  1/298.257223563


def az_adjustment(triangle):
    v = vertex[tri_v1[triangle]] # vertex ID of triangle
    c = icostriangles[triangle] # center of triangle


    # Azimuth from vertex to center of triangle
    adj = atan2(cos(v.lat) * sin(v.lon - c.lon),cos(c.lat) * sin(v.lat)- sin(c.lat) * cos(v.lat) * cos(v.lon - c.lon))

    return adj
  
def isea_snyder_forward(ll):
       
    for i in range(1,21):

        center = icostriangles[i];

        # step 1 , returns z(scaled meters) and Az (in radians) 
       
        z, Az = vincentyInverse(center.lat, center.lon, ll.lat, ll.lon)
        
        if (Az > M_PI):
            Az = Az - (2*M_PI)
       
        # not on this triangle                        
        if (z > g):
            continue
      
        
        # step 2 

        # This calculates a vertex coordinate whose azimuth is going to be assigned 0
        az_offset = az_adjustment(i)


        # This gives that vertex an azimuth of 0. For south pointing triangles the range of azimuths changes
        # from [-3.14 - 3.14] to [-6.28 0].
        # For north pointing triangles, makes no difference. 
        
        Az -= az_offset 
        
         
        # Adjust Az to fall between 0 and 120 degrees, record adjustment amount

        Az_adjust_multiples = 0;
        while Az < 0.0:
            Az += DEG120
            Az_adjust_multiples-=1

        while (Az > DEG120):
            Az -= DEG120
            Az_adjust_multiples+=1


        # Calculate q from eq 9. 
        
        q = atan(TANLOWERG/(cos(Az) + (sin(Az)*COTTHETA)))

        # not in this triangle 
        if (z > q):
            continue
        
        # Apply equations 5-8 and 10-12 in order
        # eq 6
        
        H = acos((sin(Az) * SINUPPERG * COSLOWERG) - (cos(Az) * COSUPPERG))


        #eq 7
        
        AG = (Az + G + H - DEG180)

        # eq 8 
        
        Azprime = atan2((2.0 * AG), ((RPRIME * RPRIME* TANLOWERG * TANLOWERG) - (2.0 * AG * COTTHETA)))

        # eq 10
        
        dprime = (RPRIME * TANLOWERG) / (cos(Azprime) + (sin(Azprime) * COTTHETA))
        
        # eq 11
        
        f = dprime / (2.0 * RPRIME * sin(q / 2.0))

        # eq 12
        
        rho = 2.0 * RPRIME * f * sin(z / 2.0)


        #add back the same 120 degree multiple adjustment from step 2 to Azprime
        
        Azprime += DEG120 * Az_adjust_multiples
        
        # calculate rectangular coordinates
        
        x = rho * sin(Azprime)
        y = rho * cos(Azprime)

        out = isea_pt(x,y)
        
        return out, i

def computeMorton2D(longitude, latitude, res):
    
    point = isea_geo(longitude * DEG2RAD,  latitude * DEG2RAD)
    
    # out contains coordinates from center of triangle 
    out, face = isea_snyder_forward(point)
    
    # Find new coordinates of point from lower left/upper left origin

    newPointX = 0
    newPointY= 0
    
    if ((face >= 1 and face <=5) or (face>=11 and face <=15)):
        newPointX = out.x - NEWORIGX
        newPointY = out.y - NEWORIGY
    else:
        newPointX = (out.x + NEWORIGX) * (-1)
        newPointY = (out.y - NEWORIGY) * (-1)


    # Rotate the axes, round down to nearest integer since addressing begins at 0
    # Scale coordinates of all dimensions to match resolution of DGGS

    origX = ((newPointX - ((1/(sqrt(3))) * newPointY)) /(NEWORIGX *(-2))) * totRange
    origY = ((newPointX + ((1/(sqrt(3))) * newPointY))/(NEWORIGX *(-2))) * totRange
         

    rotatedX = int(origX)
    rotatedY = int(origY) 
    
   
    # Convert to binary
    
    xBin = ('{0:0' + str(res) + 'b}').format(rotatedX)
    yBin = ('{0:0' + str(res) + 'b}').format(rotatedY)
    
    # Convert binary to int and use Morton formula
    
    morton = ('{0:0' + str(res) + '}').format(2 * int(yBin) + int(xBin))
    
    
    # Convert triangle face number to rhombus face number
    
    if face == 1 or face == 6:
        face = 0
    elif face == 11 or face == 16:
        face = 1
    elif face == 2 or face == 7:
        face = 2
    elif face == 12 or face == 17:
        face = 3
    elif face == 3 or face == 8:
        face = 4
    elif face == 13 or face == 18:
        face = 5
    elif face == 4 or face == 9:
        face = 6
    elif face == 14 or face == 19:
        face = 7
    elif face == 5 or face == 10:
        face = 8
    else:
        face = 9

    fullCode = str(face) + morton
    
    return fullCode


def computeMorton3D(longitude, latitude, height, hrange, res):
    
    point = isea_geo(longitude * DEG2RAD,  latitude * DEG2RAD)
    #out contains coordinates from center of triangle 
    out, face = isea_snyder_forward(point)
    
    # Find new coordinates of point from lower left/upper left origin

    newPointX = 0
    newPointY= 0
    
    if ((face >= 1 and face <=5) or (face>=11 and face <=15)):
        newPointX = out.x - NEWORIGX
        newPointY = out.y - NEWORIGY
    else:
        newPointX = (out.x + NEWORIGX) * (-1)
        newPointY = (out.y - NEWORIGY) * (-1)


    # Rotate the axes, round down to nearest integer since addressing begins at 0
    # Scale coordinates of all dimensions to match resolution of DGGS


    origX = ((newPointX - ((1/(sqrt(3))) * newPointY)) /(NEWORIGX *(-2))) * totRange
    origY = ((newPointX + ((1/(sqrt(3))) * newPointY))/(NEWORIGX *(-2))) * totRange
    Z = 0
    if height <=0:
        Z = ((hrange - ((-1) * height)) / (hrange * 2)) * totRange
    else:
        Z = ((hrange + height) / (hrange * 2)) * totRange
      

    rotatedX = int(origX)
    rotatedY = int(origY)
    intZ = int(Z)     
    
    
    # Convert to binary
    
    xBin = ('{0:0' + str(res) + 'b}').format(rotatedX)
    yBin = ('{0:0' + str(res) + 'b}').format(rotatedY)
    zBin = ('{0:0' + str(res) + 'b}').format(intZ)
    
    # Convert binary to int and use Morton formula
    
    morton = ('{0:0' + str(res) + '}').format(4 * int(zBin) + 2 * int(yBin) + int(xBin))
    
    
    # Convert triangle face number to rhombus face number
    
    if face == 1 or face == 6:
        face = 0
    elif face == 11 or face == 16:
        face = 1
    elif face == 2 or face == 7:
        face = 2
    elif face == 12 or face == 17:
        face = 3
    elif face == 3 or face == 8:
        face = 4
    elif face == 13 or face == 18:
        face = 5
    elif face == 4 or face == 9:
        face = 6
    elif face == 14 or face == 19:
        face = 7
    elif face == 5 or face == 10:
        face = 8
    else:
        face = 9

    fullCode =  str(face) + morton 
    
    return fullCode

def computeMorton4D(longitude, latitude, height, timeSec, hrange, trange, res):
    
    point = isea_geo(longitude * DEG2RAD,  latitude * DEG2RAD)
    #out contains coordinates from center of triangle 
    out, face = isea_snyder_forward(point)

    # Find new coordinates of point from lower left/upper left origin

    newPointX = 0
    newPointY= 0
    
    if ((face >= 1 and face <=5) or (face>=11 and face <=15)):
        newPointX = out.x - NEWORIGX
        newPointY = out.y - NEWORIGY
    else:
        newPointX = (out.x + NEWORIGX) * (-1)
        newPointY = (out.y - NEWORIGY) * (-1)


    # Rotate the axes, round down to nearest integer since addressing begins at 0
    # Scale coordinates of all dimensions to match resolution of DGGS

    origX = ((newPointX - ((1/(sqrt(3))) * newPointY)) /(NEWORIGX *(-2))) * totRange
    origY = ((newPointX + ((1/(sqrt(3))) * newPointY))/(NEWORIGX *(-2))) * totRange
    
    Z = 0
    if height <=0:
        Z = ((hrange - ((-1) * height)) / (hrange * 2)) * totRange
    else:
        Z = ((hrange + height) / (hrange * 2)) * totRange
        
    # Subtract 18 seconds to get UTC time
    T = ((timeSec - 18)/trange) * totRange



    rotatedX = int(origX)
    rotatedY = int(origY)
    intZ = int(Z)     
    intT = int(T)

   
    # Convert to binary
    
    xBin = ('{0:0' + str(res) + 'b}').format(rotatedX)
    yBin = ('{0:0' + str(res) + 'b}').format(rotatedY)
    zBin = ('{0:0' + str(res) + 'b}').format(intZ)
    tBin = ('{0:0' + str(res) + 'b}').format(intT)
    

    tVal = int('0'.join(str(8 * int(tBin))))
    zVal = int('0'.join(str(4 * int(zBin))))
    yVal = int('0'.join(str(2 * int(yBin))))
    xVal = int('0'.join(str(int(xBin))))
    
    total = tVal + zVal + yVal + xVal
    

    # Convert binary to int and use Morton formula
    # Drop 0 (if any) after face number to save 1 bit!
    
    morton3D = ('{0:0' + str(res) + '}').format(4 * int(zBin) + 2 * int(yBin) + int(xBin))
    
    morton4D = ('{0:0' + str((2* res))  + '}').format(total)


    # Convert triangle face number to rhombus face number
    
    if face == 1 or face == 6:
        face = 0
    elif face == 11 or face == 16:
        face = 1
    elif face == 2 or face == 7:
        face = 2
    elif face == 12 or face == 17:
        face = 3
    elif face == 3 or face == 8:
        face = 4
    elif face == 13 or face == 18:
        face = 5
    elif face == 4 or face == 9:
        face = 6
    elif face == 14 or face == 19:
        face = 7
    elif face == 5 or face == 10:
        face = 8
    else:
        face = 9

    fullCode4D = str(face) + morton4D 
    fullCode3D = str(face) + morton3D
    
    return fullCode4D, fullCode3D


def decodeMortonToXY(mortonCode):
    
    # Get face number and morton code
   
    face = mortonCode[0] # First number indicates rhombus face!
    morton = mortonCode[1:len(mortonCode)] # String 
    res = len(morton) #discrete!
    lastDig = int(morton[-1:])
    
    numXValues = totRange # Amount of possible X values
    numYValues = totRange # Amount of possible Y values
    
    # Compute the X, Y, and Z values on rhombus... markers = middle value
    xMarker = numXValues/2 
    yMarker = numYValues/2
    
    for i in range(res-1):    
        if (int(morton[i]) %2 ==0):
            numXValues = numXValues/2 
            xMarker = xMarker - numXValues/2 
        else:
            numXValues = numXValues/2
            xMarker = xMarker + numXValues/2
            
        if (int(morton[i]) <= 1):
            numYValues = numYValues/2
            yMarker = yMarker - numYValues/2
        else:
            numYValues = numYValues/2
            yMarker = yMarker + numYValues/2
             
    # Look at last digit

    if lastDig %2 ==0:
        xMarker = (xMarker - 1) 
    else:
        xMarker = (xMarker)

    if lastDig <=1:
        yMarker = (yMarker -1)
    else:
        yMarker = (yMarker)

    return (xMarker,yMarker,int(face), res)    

def decodeMortonToXYH(mortonCode):
    
    # Get face number and Morton code
    
    face = mortonCode[0] # First number indicates rhombus face!
    morton = mortonCode[1:len(mortonCode)] # String 
    res = len(morton) # discrete!
    lastDig = int(morton[-1:])
    

    
    numXValues = totRange # Amount of possible X values
    numYValues = totRange # Amount of possible Y values
    numZValues = totRange # Amount of possible Z values
    
    # Compute the X, Y, and Z values on rhombus... markers = middle value
    xMarker = numXValues/2 
    yMarker = numYValues/2
    zMarker = numZValues/2
    
    yVals = [0,1,4,5]

    
    for i in range(res-1):    
        if (int(morton[i]) %2 ==0):
            numXValues = numXValues/2 
            xMarker = xMarker - numXValues/2 
        else:
            numXValues = numXValues/2
            xMarker = xMarker + numXValues/2
            
        if (int(morton[i]) in yVals):
            numYValues = numYValues/2
            yMarker = yMarker - numYValues/2
        else:
            numYValues = numYValues/2
            yMarker = yMarker + numYValues/2
            
        if (int(morton[i]) <= 3):
            numZValues = numZValues/2
            zMarker = zMarker - numZValues/2
        else:
            numZValues = numZValues/2
            zMarker = zMarker + numZValues/2
             
    # Look at last digit

    if lastDig %2 ==0:  # 0,2
        xMarker = (xMarker - 1)
    else:
        xMarker = (xMarker) 

    if lastDig in yVals:
        yMarker = (yMarker -1)
    else:
        yMarker = (yMarker)

    if lastDig <= 3:
        zMarker = (zMarker -1)
    else:
        zMarker = (zMarker)

    return (xMarker,yMarker, zMarker, int(face), res)    


def decodeMortonToXYHT(mortonCode):
    
  
    # Get face number and morton code
    
    face = mortonCode[0] # First number indicates rhombus face!

    morton = mortonCode[1:len(mortonCode)] # String 
    res = len(morton) /2

    lastDig = int(morton[-2:])
    
    
    numXValues = totRange # Amount of possible X values
    numYValues = totRange# Amount of possible Y values
    numZValues = totRange # Amount of possible Z values
    numTValues = totRange # Amount of possible T values
    
    # Compute the X, Y, Z, and T values on rhombus... markers = middle value
    xMarker = numXValues/2 
    yMarker = numYValues/2
    zMarker = numZValues/2
    tMarker = numTValues/2
    
    yVals = [0,1,4,5,8,9,12,13]
    zVals = [0,1,2,3,8,9,10,11]
    
    for i in range(res-1):  
        i = i * 2
        if (int(morton[i:i+2]) %2 ==0):
            numXValues = numXValues/2 
            xMarker = xMarker - numXValues/2 
        else:
            numXValues = numXValues/2
            xMarker = xMarker + numXValues/2
            
        if (int(morton[i:i+2]) in yVals):
            numYValues = numYValues/2
            yMarker = yMarker - numYValues/2
        else:
            numYValues = numYValues/2
            yMarker = yMarker + numYValues/2
            
        if (int(morton[i:i+2]) in zVals):
            numZValues = numZValues/2
            zMarker = zMarker - numZValues/2
        else:
            numZValues = numZValues/2
            zMarker = zMarker + numZValues/2
        
        if (int(morton[i:i+2]) <= 7):
            numTValues = numTValues/2
            tMarker = tMarker - numTValues/2
        else:
            numTValues = numTValues/2
            tMarker = tMarker + numTValues/2       
     

    # Look at last digit

    if lastDig %2 ==0:
        xMarker = (xMarker - 1) 
    else:
        xMarker = (xMarker)

    if lastDig in yVals:
        yMarker = (yMarker -1)
    else:
        yMarker = (yMarker)

    if lastDig in zVals:
        zMarker = (zMarker -1)
    else:
        zMarker = (zMarker)

    if lastDig <=7:
        tMarker = (tMarker -1) 
    else:
        tMarker = (tMarker)   
        

    return (xMarker,yMarker, zMarker, tMarker, int(face), res)   

def MortonToLatLong2D(x,y, face, res):
     
    # Scale coordinates to scale of Cartesian system
    
    scaledX = (x/totRange) * (-NEWORIGX *2)
    scaledY = (y/totRange)*(-NEWORIGX*2)


    # Convert coordinates from skewed system to Cartesian system (origin at left)
    
    a = np.array([[1,(-1/sqrt(3))], [1,(1/sqrt(3))]])
    b = np.array([scaledX,scaledY]) 
    x = np.linalg.solve(a, b)
    
    xCoord = x[0]
    yCoord = x[1]

    
    # Get triangle face from rhombus face based on values of y.
    # If y is negative, triangles will be downward oriented
    
    if yCoord >=0:
        if (face == 0):
            face = 1
        elif (face == 1):
            face = 11
        elif (face == 2):
            face = 2
        elif (face == 3):
            face = 12
        elif (face == 4):
            face = 3
        elif (face == 5):
            face = 13
        elif (face == 6):
            face = 4
        elif (face == 7):
            face = 14
        elif (face == 8):
            face = 5
        elif (face == 9):
            face = 15
    else:
        if (face == 0):
            face = 6
        elif (face == 1):
            face = 16
        elif (face == 2):
            face = 7
        elif (face == 3):
            face = 17
        elif (face == 4):
            face = 8
        elif (face == 5):
            face = 18
        elif (face == 6):
            face = 9
        elif (face == 7):
            face = 19
        elif (face == 8):
            face = 10
        elif (face == 9):
            face = 20

    # Translate coordinates to center (origin) of icosahedron triangle, 
    # taking into account triangle orientation 
    
    xOrigin = 0
    yOrigin = 0
    
    if ((face >= 1 and face <=5) or (face>=11 and face <=15)):
        xOrigin = ((-NEWORIGX) - xCoord) * (-1)
        yOrigin = ((-NEWORIGY) - yCoord) * (-1)
    else:
        xOrigin = (-NEWORIGX) - xCoord
        yOrigin = ((-NEWORIGY) + yCoord) * (-1) 
  

    # Equation 17
    
    Azprime = atan2(xOrigin, yOrigin)
    
    # Equation 18
    
    rho = sqrt((pow(xOrigin,2) + pow(yOrigin,2)))

    # Adjust Azprime to fall within 0 to 120 degrees
    
    Azprime_adjust_multiples = 0;
    while Azprime < 0.0:
        Azprime += DEG120
        Azprime_adjust_multiples-=1
    while (Azprime > DEG120):
        Azprime -= DEG120
        Azprime_adjust_multiples+=1


    AzprimeCopy = Azprime

    # Equation 19
    
    AG = (pow(RPRIME,2) * pow(TANLOWERG,2)) / (2 * ((1/(tan(Azprime))) + COTTHETA))


    # Iteration, Azprime (plane) converges to Az (sphere)
    
    for i in range(4):
    
        H = acos((sin(Azprime) * SINUPPERG * COSLOWERG) - (cos(Azprime) * COSUPPERG))

        FAZ = (AG - G - H - Azprime + DEG180)
    
        FPRIMEAZ = (((cos(Azprime) * SINUPPERG * COSLOWERG) + (sin(Azprime) * COSUPPERG)) / (sin(H))) - 1
    
        DeltaAzprime = -FAZ/(FPRIMEAZ)
        
        Azprime = Azprime + DeltaAzprime

    Az = Azprime
 

    # Equations 9-11, 23 to obtain z 
    
    q = atan((TANLOWERG)/(cos(Az) + (sin(Az)*COTTHETA)))

    # eq 10

    dprime = ((RPRIME * TANLOWERG) / (cos(AzprimeCopy) + (sin(AzprimeCopy) * COTTHETA)))

    # eq 11
    
    f = dprime / (2.0 * RPRIME * sin(q / 2.0))

    # eq 23, obtain z
    
    z = 2 * asin((rho)/(2*RPRIME*f))
     
    # Add back 120 degree adjustments to Az

    Az += DEG120 * Azprime_adjust_multiples
     
    # Adjust Az to be clockwise from north (needed for final calculation)
    if (face >=1 and face<=5) or (face>=11 and face<=15):
        if (Az <0):
            Az = (M_PI - (Az * (-1))) + M_PI 
    else:
        if (Az <0):
            Az = M_PI - (Az * (-1))
        else:
            Az = Az + M_PI
    
    z = z * R
    
    # triangle center
    center = icostriangles[face]
    
    lat2, lon2 = vincentyDirect(flattening, R, center.lat * RAD2DEG, center.lon * RAD2DEG, Az * RAD2DEG, z) 
    return lat2, lon2


def MortonToLatLong3D(x,y,h,face, res):
     
    # Convert h/Z to height above/below ellipsoid
    
    height = ((-1) * hrange) + ((h / totRange) * (2 * hrange))


    # Scale coordinates to scale of Cartesian system
    
    scaledX = (x/totRange) * (-NEWORIGX *2)
    scaledY = (y/totRange)*(-NEWORIGX*2)


    # Convert coordinates from skewed system to Cartesian system (origin at left)
    
    a = np.array([[1,(-1/sqrt(3))], [1,(1/sqrt(3))]])
    b = np.array([scaledX,scaledY]) 
    x = np.linalg.solve(a, b)
    
    xCoord = x[0]
    yCoord = x[1]

    
    # Get triangle face from rhombus face based on values of y.
    # If y is negative, triangles will be downward oriented
    
    if yCoord >=0:
        if (face == 0):
            face = 1
        elif (face == 1):
            face = 11
        elif (face == 2):
            face = 2
        elif (face == 3):
            face = 12
        elif (face == 4):
            face = 3
        elif (face == 5):
            face = 13
        elif (face == 6):
            face = 4
        elif (face == 7):
            face = 14
        elif (face == 8):
            face = 5
        elif (face == 9):
            face = 15
    else:
        if (face == 0):
            face = 6
        elif (face == 1):
            face = 16
        elif (face == 2):
            face = 7
        elif (face == 3):
            face = 17
        elif (face == 4):
            face = 8
        elif (face == 5):
            face = 18
        elif (face == 6):
            face = 9
        elif (face == 7):
            face = 19
        elif (face == 8):
            face = 10
        elif (face == 9):
            face = 20

    # Translate coordinates to center (origin) of icosahedron triangle, 
    # taking into account triangle orientation 
    
    xOrigin = 0
    yOrigin = 0
    
    if ((face >= 1 and face <=5) or (face>=11 and face <=15)):
        xOrigin = ((-NEWORIGX) - xCoord) * (-1)
        yOrigin = ((-NEWORIGY) - yCoord) * (-1)
    else:
        xOrigin = (-NEWORIGX) - xCoord
        yOrigin = ((-NEWORIGY) + yCoord) * (-1) 
  

    # Equation 17
    
    Azprime = atan2(xOrigin, yOrigin)
    
    # Equation 18
    
    rho = sqrt((pow(xOrigin,2) + pow(yOrigin,2)))

    # Adjust Azprime to fall within 0 to 120 degrees
    
    Azprime_adjust_multiples = 0;
    while Azprime < 0.0:
        Azprime += DEG120
        Azprime_adjust_multiples-=1
    while (Azprime > DEG120):
        Azprime -= DEG120
        Azprime_adjust_multiples+=1


    AzprimeCopy = Azprime

    #Equation 19
    
    AG = (pow(RPRIME,2) * pow(TANLOWERG,2)) / (2 * ((1/(tan(Azprime))) + COTTHETA))


    # Iteration, Azprime (plane) converges to Az (ellipsoid)
    
    for i in range(4):
    
        H = acos((sin(Azprime) * SINUPPERG * COSLOWERG) - (cos(Azprime) * COSUPPERG))

        FAZ = (AG - G - H - Azprime + DEG180)
    
        FPRIMEAZ = (((cos(Azprime) * SINUPPERG * COSLOWERG) + (sin(Azprime) * COSUPPERG)) / (sin(H))) - 1
    
        DeltaAzprime = -FAZ/(FPRIMEAZ)
        
        Azprime = Azprime + DeltaAzprime

    Az = Azprime
 

    # Equations 9-11, 23 to obtain z 
    
    q = atan((TANLOWERG)/(cos(Az) + (sin(Az)*COTTHETA)))

    # eq 10

    dprime = ((RPRIME * TANLOWERG) / (cos(AzprimeCopy) + (sin(AzprimeCopy) * COTTHETA)))


    # eq 11
    
    f = dprime / (2.0 * RPRIME * sin(q / 2.0))

    
    #eq 23, obtain z
    
    z = 2 * asin((rho)/(2*RPRIME*f))
     
    # Add back 120 degree adjustments to Az

    Az += DEG120 * Azprime_adjust_multiples
     
    # Adjust Az to be clockwise from north (needed for final calculation)
    if (face >=1 and face<=5) or (face>=11 and face<=15):
        if (Az <0):
            Az = (M_PI - (Az * (-1))) + M_PI 
    else:
        if (Az <0):
            Az = M_PI - (Az * (-1))
        else:
            Az = Az + M_PI
    
    z = z * R
    
    # triangle center
    center = icostriangles[face]
    
    lat2, lon2 = vincentyDirect(flattening, R, center.lat * RAD2DEG, center.lon * RAD2DEG, Az * RAD2DEG, z) 
    
    return lat2, lon2, height



def MortonToLatLong4D(x,y,h,t, face, res):
         
    # Convert h to height above sphere/ellipsoid
    
    height = ((-1) * hrange) + ((h /  totRange) * (2 * hrange))


    # Convert t to time in seconds since GPS Time - Jan 6, 1980 12 AM
    # add GPS-UTC offset back
    
    timeSec = round((t/totRange) * trange) + 18
    

    # Scale coordinates to scale of Cartesian system
    
    scaledX = (x/ totRange) * (-NEWORIGX *2)
    scaledY = (y/ totRange)*(-NEWORIGX*2)


    # Convert coordinates from skewed system to Cartesian system 
    # (origin at left of rhombus)
    
    a = np.array([[1,(-1/sqrt(3))], [1,(1/sqrt(3))]])
    b = np.array([scaledX,scaledY]) 
    x = np.linalg.solve(a, b)
    
    xCoord = x[0]
    yCoord = x[1]

    # Get triangle face from rhombus face based on values of y.
    # If y is negative, triangles will be downward oriented

    if yCoord >=0:
        if (face == 0):
            face = 1
        elif (face == 1):
            face = 11
        elif (face == 2):
            face = 2
        elif (face == 3):
            face = 12
        elif (face == 4):
            face = 3
        elif (face == 5):
            face = 13
        elif (face == 6):
            face = 4
        elif (face == 7):
            face = 14
        elif (face == 8):
            face = 5
        elif (face == 9):
            face = 15
    else:
        if (face == 0):
            face = 6
        elif (face == 1):
            face = 16
        elif (face == 2):
            face = 7
        elif (face == 3):
            face = 17
        elif (face == 4):
            face = 8
        elif (face == 5):
            face = 18
        elif (face == 6):
            face = 9
        elif (face == 7):
            face = 19
        elif (face == 8):
            face = 10
        elif (face == 9):
            face = 20

    # Translate coordinates to center (origin) of icosahedron triangle, 
    # taking into account triangle orientation 
    
    xOrigin = 0
    yOrigin = 0
    
    if ((face >= 1 and face <=5) or (face>=11 and face <=15)):
        xOrigin = ((-NEWORIGX) - xCoord) * (-1)
        yOrigin = ((-NEWORIGY) - yCoord) * (-1)
    else:
        xOrigin = (-NEWORIGX) - xCoord
        yOrigin = ((-NEWORIGY) + yCoord) * (-1) 
  

    # Equation 17
    
    Azprime = atan2(xOrigin, yOrigin)
    
    # Equation 18
    
    rho = sqrt((pow(xOrigin,2) + pow(yOrigin,2)))

    # Adjust Azprime to fall within 0 to 120 degrees
    
    Azprime_adjust_multiples = 0;
    while Azprime < 0.0:
        Azprime += DEG120
        Azprime_adjust_multiples-=1
    while (Azprime > DEG120):
        Azprime -= DEG120
        Azprime_adjust_multiples+=1


    AzprimeCopy = Azprime

    #Equation 19
    
    AG = (pow(RPRIME,2) * pow(TANLOWERG,2)) / (2 * ((1/(tan(Azprime))) + COTTHETA))


    # Iteration, Azprime (plane) converges to Az (sphere)
    
    for i in range(4):
    
        H = acos((sin(Azprime) * SINUPPERG * COSLOWERG) - (cos(Azprime) * COSUPPERG))

        FAZ = (AG - G - H - Azprime + DEG180)
    
        FPRIMEAZ = (((cos(Azprime) * SINUPPERG * COSLOWERG) + (sin(Azprime) * COSUPPERG)) / (sin(H))) - 1
    
        DeltaAzprime = -FAZ/(FPRIMEAZ)
        
        Azprime = Azprime + DeltaAzprime

    Az = Azprime
 

    # Equations 9-11, 23 to obtain z 
    
    q = atan((TANLOWERG)/(cos(Az) + (sin(Az)*COTTHETA)))

    # eq 10

    dprime = ((RPRIME * TANLOWERG) / (cos(AzprimeCopy) + (sin(AzprimeCopy) * COTTHETA)))


    # eq 11
    
    f = dprime / (2.0 * RPRIME * sin(q / 2.0))

    
    #eq 23, obtain z
    
    z = 2 * asin((rho)/(2*RPRIME*f))
    
    # Add back 120 degree adjustments to Az

    Az += DEG120 * Azprime_adjust_multiples
    
    # Adjust Az to be clockwise from north (needed for final calculation)
    if (face >=1 and face<=5) or (face>=11 and face<=15):
        if (Az <0):
            Az = (M_PI - (Az * (-1))) + M_PI 
    else:
        if (Az <0):
            Az = M_PI - (Az * (-1))
        else:
            Az = Az + M_PI
    
    z = z * R

   
    # triangle center
    center = icostriangles[face]

    lat2, lon2 = vincentyDirect(flattening, R, center.lat * RAD2DEG, center.lon * RAD2DEG, Az * RAD2DEG, z) 
    return lat2, lon2, height, timeSec


def resolution(beamDiam):
    
    # area, square meters
    area = ((M_PI/4) * pow(beamDiam,2))
    
    # area in square millimeters = importance of point
    
    areamm = area * 1000000
    
    # find closest value in cell areas array
    
    res = (np.abs(cellAreas-area)).argmin()
    return res, areamm


def  vincentyInverse(lat1, lon1, lat2, lon2) : 
    # all values in radians!

    if (lon1 == -M_PI):
        lon1 = M_PI

    L = lon2 - lon1
    tanU1 = (1-flattening) * tan(lat1)
    cosU1 = 1 / sqrt((1 + tanU1*tanU1))
    sinU1 = tanU1 * cosU1
    tanU2 = (1-flattening) * tan(lat2)
    cosU2 = 1 / sqrt((1 + tanU2*tanU2))
    sinU2 = tanU2 * cosU2

    sinS=0
    cosS=0
    Sigma=0
    cosSqAlpha=0
    cos2SigmaM=0
    

    Lambda = L
    LambdaPrime = 0
    iterations = 0
    while (abs(Lambda-LambdaPrime) > 1e-12 and iterations<1000):
        iterations += 1
        sinL = sin(Lambda)
        cosL = cos(Lambda)
        sinSqSigma = (cosU2*sinL) * (cosU2*sinL) + (cosU1*sinU2-sinU1*cosU2*cosL) * (cosU1*sinU2-sinU1*cosU2*cosL)
        if (sinSqSigma == 0):
            break; #co-incident points
        sinS = sqrt(sinSqSigma)
        cosS = sinU1*sinU2 + cosU1*cosU2*cosL
        Sigma = atan2(sinS, cosS)
        sinAlpha = cosU1 * cosU2 * sinL / sinS
        cosSqAlpha = 1 - sinAlpha*sinAlpha
        if (cosSqAlpha != 0):
            cos2SigmaM = (cosS - 2*sinU1*sinU2/cosSqAlpha)
        else:
            cos2SigmaM = 0
        
        C = flattening/16*cosSqAlpha*(4+flattening*(4-3*cosSqAlpha))
        LambdaPrime = Lambda;
        Lambda = L + (1-C) * flattening * sinAlpha * (Sigma + C*sinS*(cos2SigmaM+C*cosS*(-1+2*cos2SigmaM*cos2SigmaM)))
        if (abs(Lambda) > M_PI):
            print('Lambda > π')
    
    if (iterations >= 1000):
        print('Formula failed to converge')

    uSq = cosSqAlpha * (R*R - b*b) / (b*b)
    A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
    B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
    DeltaSigma = B*sinS*(cos2SigmaM+B/4*(cosS*(-1+2*cos2SigmaM*cos2SigmaM)-B/6*cos2SigmaM*(-3+4*sinS*sinS)*(-3+4*cos2SigmaM*cos2SigmaM)))

    s = b*A*(Sigma-DeltaSigma)

    Alpha1 = atan2(cosU2*sinL,  cosU1*sinU2-sinU1*cosU2*cosL)
    Alpha2 = atan2(cosU1*sinL, -sinU1*cosU2+cosU1*sinU2*cosL)

    Alpha1 = (Alpha1 + 2*M_PI) % (2*M_PI); # normalise to 0..360
    Alpha2 = (Alpha2 + 2*M_PI) % (2*M_PI); # normalise to 0..360

    return s/R, Alpha1


def  vincentyDirect(f, a, phi1, lembda1, alpha12, s ) : 
        """ 
        Returns the lat and long of projected point
        given a reference point and a distance and azimuth to project. 
        Returns ( phi2,  lambda2) as a tuple, in decimal degrees. 
        Parameters:

        The code has been originally taken from
        https://isis.astrogeology.usgs.gov/IsisSupport/index.php?topic=408.0 in Javascript,
        and later converted into Python. 
        """ 
        piD4 = math.atan( 1.0 ) 
        two_pi = piD4 * 8.0 
        phi1    = phi1    * piD4 / 45.0 
        lembda1 = lembda1 * piD4 / 45.0 
        alpha12 = alpha12 * piD4 / 45.0 
        if ( alpha12 < 0.0 ) : 
            alpha12 = alpha12 + two_pi 
        if ( alpha12 > two_pi ) : 
            alpha12 = alpha12 - two_pi
        b = a * (1.0 - f) 
        TanU1 = (1-f) * math.tan(phi1) 
        U1 = math.atan( TanU1 ) 
        sigma1 = math.atan2( TanU1, math.cos(alpha12) ) 
        Sinalpha = math.cos(U1) * math.sin(alpha12) 
        cosalpha_sq = 1.0 - Sinalpha * Sinalpha 
        u2 = cosalpha_sq * (a * a - b * b ) / (b * b) 
        A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
            (320 - 175 * u2) ) ) 
        B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) ) 
        # Starting with the approx 
        sigma = (s / (b * A)) 
        last_sigma = 2.0 * sigma + 2.0   # something impossible 
            
        # Iterate the following 3 eqs unitl no sig change in sigma 
        # two_sigma_m , delta_sigma 
        while ( abs( (last_sigma - sigma) / sigma) > 1.0e-9 ):
            two_sigma_m = 2 * sigma1 + sigma 
            delta_sigma = B * math.sin(sigma) * ( math.cos(two_sigma_m) \
                    + (B/4) * (math.cos(sigma) * \
                    (-1 + 2 * math.pow( math.cos(two_sigma_m), 2 ) -  \
                    (B/6) * math.cos(two_sigma_m) * \
                    (-3 + 4 * math.pow(math.sin(sigma), 2 )) *  \
                    (-3 + 4 * math.pow( math.cos (two_sigma_m), 2 )))))
            last_sigma = sigma 
            sigma = (s / (b * A)) + delta_sigma 
        phi2 = math.atan2 ( (math.sin(U1) * math.cos(sigma) +\
            math.cos(U1) * math.sin(sigma) * math.cos(alpha12) ), \
            ((1-f) * math.sqrt( math.pow(Sinalpha, 2) +  \
            pow(math.sin(U1) * math.sin(sigma) - math.cos(U1) * \
            math.cos(sigma) * math.cos(alpha12), 2))))
        lembda = math.atan2( (math.sin(sigma) * math.sin(alpha12 )),\
            (math.cos(U1) * math.cos(sigma) -  \
            math.sin(U1) *  math.sin(sigma) * math.cos(alpha12))) 
        C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq )) 
        omega = lembda - (1-C) * f * Sinalpha *  \
            (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) + \
            C * math.cos(sigma) * (-1 + 2 *\
            math.pow(math.cos(two_sigma_m), 2) ))) 
        lembda2 = lembda1 + omega 
        alpha21 = math.atan2 ( Sinalpha, (-math.sin(U1) * \
            math.sin(sigma) +
            math.cos(U1) * math.cos(sigma) * math.cos(alpha12))) 
        alpha21 = alpha21 + two_pi / 2.0 
        if ( alpha21 < 0.0 ) : 
            alpha21 = alpha21 + two_pi 
        if ( alpha21 > two_pi ) : 
            alpha21 = alpha21 - two_pi 
        phi2 = phi2 * 45.0 / piD4 
        lembda2 = lembda2 * 45.0 / piD4 
        
        return phi2, lembda2

if __name__ == "__main__":

    # DGGS full-resolution number
    maxRes=33  # sub-mllimeter precision
    
    # Specify Z range above/below surface of Earth (in meters)
    hrange = 5000.0
    timeSec = 1140874544
    
    # Specify T range (seconds, GPS Time) Jan 6, 1980 -- Jan 6, 2018
    trange = 1199232018.0
    
    # Total range of values for dimensions
    totRange = pow(2,maxRes)

    # Connect to database
    
    try:
        conn = psycopg2.connect("dbname='dggs' user='postgres' host='localhost' password='serengeti'")
        # Commits every transaction by default 
        conn.autocommit = True
    except:
        print "I am unable to connect to the database"

    cur = conn.cursor()

    # Read ISEA4D stats table, store in list
    # This is a table containing statistics on each resolution in the DGGS

    df = pandas.read_excel(r'C:\Neeraj\DGGS\ISEA4D\ISEAStats.xlsx', sheet_name='Sheet1')
    mat = df.as_matrix()
    
    # 4th column contains cell areas in square meters
    cellAreas = mat[:,4]
     
    p1 = Proj(init='EPSG:28992')
    p2 = Proj(proj='latlong',datum='WGS84')
    
    factor = 2.777777777777778 # Footprint size multiplication factor for Riegl scanner 
    
    # Path to folder containing LAS files, set as appropriate
    path = "D:\\PointCloudData\\Tiled5000_2016"

    # SQL string into which values to insert will be substituted
    SQL2 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
    
    # Variable to store number of point currently being processed
    counter = 0

    # Pre-allocate storage in memory for posting 5 million points at once
    bigList = [0]* (5000000)
   
    start_time1 = time.time()
    
    # Loop through LAS files 
    for fn in os.listdir(path):
        try:
            loc = path + "\\" + fn
            inFile = File(loc, mode = "r")
        except (RuntimeError, TypeError, NameError,Exception):
            continue
        
        x = inFile.x
        y = inFile.y
        z = inFile.z
        timeGPS = inFile.gps_time       
        
        # Distances from scanner in meters and millimeters
        # were stored in Red and Green fields of LAS file 
        disM = inFile.Red
        disMM = inFile.Green
      
        for i in range(len(x)):
            
            # Find closest discrete resolution to point
            
            totDisM = disM[i] + (disMM[i] / 1000.0)
            bd = (totDisM / factor) /1000.0  # beam diameter in meters
            resPoint, contPrec = resolution(bd)
            
            # Project points to lat/long from RD

            lon, lat, height = transform(p1, p2, x[i], y[i], z[i])

            # Find Morton code
            
            code4D, code3D = computeMorton4D(lon, lat, height, timeGPS[i], hrange, trange, maxRes)

            # Create tuple of values 
                
            data = (code3D, code4D, lon,lat,height,resPoint,contPrec)
           
            addStr = SQL2 % data

            bigList[counter] = addStr
            counter += 1

        inFile.close()

        if counter == 5000000:
            
            # Convert big list of insert values into a large string
            bigStr = ''.join(str(v) for v in bigList)
            
            # Create a StringIO object in memory. This is required by COPY_FROM
            output = StringIO(bigStr.decode('utf8'))
            print("Processing Time for 5 million points %s seconds ---" % (time.time() - start_time1))
            
            start_time2 = time.time()
            
            # This is the fastest existing function in Psycopg2 to bulk-load data from Python to
            # PostgreSQL: COPY_FROM
            cur.copy_from(output, 'pointsnew')

            output.close()
            
            print("Post time for 5 million points %s seconds ---" % (time.time() - start_time2))
            
            # Reset counter variable
            counter = 0
            
            # Preallocate space for the next 5 million points 
            bigList = [0]* 5000000

            start_time1 = time.time()

