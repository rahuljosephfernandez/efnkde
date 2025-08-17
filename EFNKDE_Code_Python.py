# "PointArr" is the array carrying the input data (spatial distribution of points). The same used in the code below is only a place holder and does not correspond to the data used in the paper for generating any of the graphs.
# If the data is given in terms of latitudes and longitudes, then a conversion to 2D Cartesian form is needed before using the below code.
# If there are forbidden regions, then the corners of the polygon are to be entered in "BoundaryArray". If there are no forbidden regions, then uncomment the line "BoundaryArray = []" in the code below.
# The function EFNKDE(x,y,BW) returns K(x, y) with bandwidth BW.

import numpy as np
from scipy.spatial import Delaunay
import math 

PointArr = np.array([[4.17,4.29], [4.22,9.19], [-6.94,2.41], [-7.18,5.15], [-2.78,1.56], [-5.31,0.04], [3.37,9.44]])
BoundaryArray = [[-2.0,-2.0],[2.0,-2.0],[2.0,2.0],[-2.0,2.0]]
#BoundaryArray = []

tri = Delaunay(PointArr)
ListOfTriangles = tri.simplices

def SortVertices(a,b):
    return [a,b] if a < b else [b,a]

AllEdgeList=[]
for i in ListOfTriangles:
    AllEdgeList.append(SortVertices(i[0],i[1]))
    AllEdgeList.append(SortVertices(i[1],i[2]))
    AllEdgeList.append(SortVertices(i[0],i[2]))
    
AllEdgeList1=np.unique(AllEdgeList, axis=0)

def FindabcGivenPts(x1,y1,x2,y2):
    return [y2 - y1, x1 - x2, (x2 - x1)*y2 + (y1 - y2)*x2]

def FindIntrsctnPtGivenLines(a1,b1,c1,a2,b2,c2):
    return [(b2*c1-b1*c2)/(a2*b1-a1*b2),(a1*c2-a2*c1)/(a2*b1-a1*b2)]

def FindIfPtInt(PtA,PtB,PtX): 
    xA = PtA[0]
    xB = PtB[0]
    xX = PtX[0]
    yA = PtA[1]
    yB = PtB[1]
    yX = PtX[1]
    if xA == xB: 
        AlphaIntExt = (yX - yA)/(yB - yA) 
    else:
        AlphaIntExt = (xX - xA)/(xB - xA)
    if(AlphaIntExt > 0) and (AlphaIntExt < 1): 
        IntAns = 1
    else:
        IntAns = 0
    return IntAns


def RemoveOverlapWithBoundary(BArr): 
    BArrTmp = BArr
    NB = len(BArr)
    BArrTmp.extend([BArr[0]])
    RemoveCrossTalksArr = []
   
    for iK in range(len(AllEdgeList1)):
        Edgeval = AllEdgeList1[iK]
        a = Edgeval[0]
        b = Edgeval[1]
        Coorda = PointArr[a]
        Coordb = PointArr[b]
        abcEdge = FindabcGivenPts(Coorda[0], Coorda[1], Coordb[0], Coordb[1])
    
        for iB in range(NB):
            Coordc = BArrTmp[iB]
            Coordd = BArrTmp[iB + 1]
            abcBoundary = FindabcGivenPts(Coordc[0], Coordc[1], Coordd[0], Coordd[1])
            IntPt = FindIntrsctnPtGivenLines(abcEdge[0], abcEdge[1], abcEdge[2], abcBoundary[0], abcBoundary[1], abcBoundary[2])
            if (FindIfPtInt(Coorda,Coordb,IntPt)==1) and (FindIfPtInt(Coordc,Coordd,IntPt)==1):
               RemoveCrossTalksArr.append([Edgeval[0],Edgeval[1]])
               break
    
    RemoveCrossTalksArr = np.array(RemoveCrossTalksArr)
    UpdatedAllEdgeList=[]

    for edge in AllEdgeList1:
        AFflag=1
        for forbedge in RemoveCrossTalksArr:
            if (edge[0]==forbedge[0]) and (edge[1]==forbedge[1]):
                AFflag=0
        if AFflag==1:
            UpdatedAllEdgeList.append([edge[0],edge[1]])
    
    return np.array(UpdatedAllEdgeList)

if len(BoundaryArray) > 0:
    AllEdgeList1 = RemoveOverlapWithBoundary(BoundaryArray)


nAllEdgeList1 = len(AllEdgeList1)

def NormOfVec(V):
    return math.sqrt(np.dot(V,V))

def IntInfoGivenEdgePT(x,y,iEdge):
    PtX=(x,y)
    ipta = AllEdgeList1[iEdge,0]
    iptb = AllEdgeList1[iEdge,1]
    Pta = PointArr[ipta]
    Ptb = PointArr[iptb]
    VXma = PtX - Pta
    VXmb = PtX - Ptb
    Vbma = Ptb - Pta
    VXmaDOTVbma = np.dot(VXma,Vbma)
    VXmbDOTVamb = np.dot(VXmb,-Vbma)
    FootOfThePerpendicular = Pta + (VXmaDOTVbma/(np.dot(Vbma,Vbma)))*Vbma
    LengthOfPerp = NormOfVec(PtX - FootOfThePerpendicular)
    if VXmaDOTVbma <= 0.0: 
        yA = NormOfVec(Pta - FootOfThePerpendicular)
    elif VXmbDOTVamb <= 0.0: 
        yA = NormOfVec(Ptb - FootOfThePerpendicular)
    else:
        yA = -NormOfVec(Pta-FootOfThePerpendicular)
    yB = yA + NormOfVec(Vbma)
    return [LengthOfPerp,yA,yB]

def PartialEFNKDE(x,y,BW,iEdgeEF):
    IntInfoval = IntInfoGivenEdgePT(x,y,iEdgeEF)
    sval = IntInfoval[0]
    yAval = IntInfoval[1]
    yBval = IntInfoval[2]
    LEdgeval = yBval - yAval
    Expval = math.exp(-0.5*((sval/BW)**2))
    Denomval = 2*math.sqrt(2*math.pi)*nAllEdgeList1*BW*LEdgeval
    erfsval = math.erf(yBval/(BW*math.sqrt(2)))-math.erf(yAval/(BW*math.sqrt(2)))
    AnsvalEF = Expval * erfsval / Denomval
    return AnsvalEF
    
def EFNKDE(x,y,BW):
    Ans=0
    for iE in range(nAllEdgeList1):
        Ans = Ans + PartialEFNKDE(x,y,BW,iE)
    return Ans
