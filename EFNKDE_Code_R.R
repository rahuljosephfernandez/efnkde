library(deldir)
library(pracma)

sortvertices <- function(a, b) {
  
  if (a < b) {
    
    return(c(a, b))
    
  } else {
    
    return(c(b, a))
    
  }
  
}

getedges <- function(data) {
  
  dat <- deldir(data[,1], data[,2], round = FALSE, digits = 10)
  
  triangles <- deldir::triang.list(dat)
  
  tri <- data.frame()
  
  for (i in 1:length(triangles)) {
    
    tri <- rbind(tri, c(triangles[[i]]$ptNum[1], triangles[[i]]$ptNum[2], triangles[[i]]$ptNum[3]))
    
  }
  
  colnames(tri) <- c("a", "b", "c")
  
  holder <- data.frame()
  
  for (i in 1:nrow(tri)) {
    
    holder <- rbind(holder, c(sortvertices(tri[i,1], tri[i,2])[1], sortvertices(tri[i,1], tri[i,2])[2]))
    
    holder <- rbind(holder, c(sortvertices(tri[i,2], tri[i,3])[1], sortvertices(tri[i,2], tri[i,3])[2]))
    
    holder <- rbind(holder, c(sortvertices(tri[i,1], tri[i,3])[1], sortvertices(tri[i,1], tri[i,3])[2]))
    
  }
  
  colnames(holder) <- c("a", "b")
  
  holder <- unique(holder)
  
  return(holder)
  
}

PointArr <- data.frame(x = c(4.17, 4.22, -6.94, -7.18, -2.78, -5.31, 3.37),
                       y = c(4.29, 9.19, 2.41, 5.15, 1.56, 0.04, 9.44))
BoundaryArray <- data.frame(x = c(-2.0, 2.0, 2.0, -2.0),
                            y = c(-2.0, -2.0, 2.0, 2.0))

AllEdgeList <- getedges(data = PointArr)

FindabcGivenPts <- function(x1, y1, x2, y2) {
  return(c(y2 - y1, x1 - x2, (x2-x1)*y2)+(y1-y2)*x2)
}

FindIntrsctnPtGivenLines <- function(a1,b1,c1,a2,b2,c2) {
  return(c(b2*c1 - b1*c2)/(a2*b1 - a1*b2),
         (a1*c2 - a2*c1)/(a2*b1 - b1*b2))
}

FindIfPtInt <- function(PtA, PtB, PtX) {
  xA <- PtA[1]
  xB <- PtB[1]
  xX <- PtX[1]
  yA <- PtA[2]
  yB <- PtB[2]
  yX <- PtX[2]
  if (xA == xB) {
    AlphaIntExt <- (yX - yA)/(yB - yA)
  } else {
    AlphaIntExt <- (xX - xA)/(xB - xA)
  }
  
  if (AlphaIntExt > 0 & AlphaIntExt < 1) {
    IntAns = 1
  } else {
    IntAns = 0
  }
  return(IntAns)
}

RemoveOverlapWithBoundary <- function(BArr) {
  BArrTmp <- BArr
  NB = length(BArr)
  BArrTmp <- rbind(BArrTmp, BArr[1])
  RemoveCrossTalksArr <- c()
  
  for (iK in 1:length(AllEdgeList)) {
    Edgeval <- AllEdgeList[iK]
    a = Edgeval[1]
    b = Edgeval[2]
    Coorda <- PointArr[a]
    Coordb <- PointArr[b]
    abcEdge <- FindabcGivenPts(Coorda[1], Coorda[2], Coordb[1], Coordb[2])
    
    for (iB in 1:NB) {
      Coordc <- BArrTmp[iB]
      Coordd <- BArrTmp[iB + 1]
      abcBoundary <- FindabcGivenPts(Coordc[1], Coordc[2], Coordd[1], Coordd[2])
      IntPt <- FindIntrsctnPtGivenLines(abcEdge[1], abcEdge[2], abcEdge[3],
                                        abcBoundary[1], abcBoundary[2], abcBoundary[3])
      if (FindIfPtInt(Coorda, Coordb, IntPt) == 1 & FindIfPtInt(Coordc, Coordd, IntPt) == 1) {
        RemoveCrossTalksArr <- c(RemoveCrossTalksArr, Edgeval[1], Edgeval[2])
        break
      }
    }
  }
  UpdatedAllEdgeList <- c()
  for (edge in AllEdgeList) {
    AFflag <- 1
    for (forbedge in RemoveCrossTalksArr){
      if (edge[1] == forbedge[1] & edge[2] == forbedge[2]){
        AFflag <- 0
      }
    }
    if (AFflag == 1){
      UpdatedAllEdgeList <- c(UpdatedAllEdgeList,edge[1]),edge[2])
    }
  }
  return(UpdatedAllEdgeList)
}

if (length(BoundaryArray) > 0) {
  AllEdgeList1 <- RemoveOverlapWithBoundary(BoundaryArray)
}

nAllEdgeList1 <- length(AllEdgeList1)

vecnorm <- function(V) {
  
  norm <- sqrt(sum(V*V))
  
  return(norm)
  
}


intinfogivenedge <- function(x, y, edges, points, iedge) {
  
  ptx <- c(x, y)
  
  ipta <- edges[iedge, 1]
  
  iptb <- edges[iedge, 2]
  
  pta <- c(points[ipta, 1], points[ipta, 2])
  
  ptb <- c(points[iptb, 1], points[iptb, 2])
  
  vxma <- ptx - pta
  
  vxmb <- ptx - ptb
  
  vbma <- ptb - pta
  
  vxma.vbma <- sum(vxma*vbma)
  
  vxmb.vamb <- sum(vxmb * -1*vbma)
  
  footoftheperp <- pta + ((vxma.vbma/(sum(vbma*vbma))) * vbma)
  
  lengthofperp <- vecnorm(V = ptx - footoftheperp)
  
  if (vxma.vbma <= 0) {
    
    ya <- vecnorm(V = pta - footoftheperp)
    
  } else if (vxmb.vamb <= 0) {
    
    ya <- vecnorm(V = ptb - footoftheperp)
    
  } else {
    
    ya = -1 * vecnorm(V = ptb - footoftheperp)
    
  }
  
  yb <- ya + vecnorm(V = vbma)
  
  return(c(lengthofperp, ya, yb))
  
}

partial.efnkde <- function(x, y, bw, points, edges, iedgeef) {
  
  intinfovals <- intinfogivenedge(x = x, y = y, points = points, edges = edges, iedge = iedgeef)
  
  s.val <- intinfovals[1]
  
  ya.val <- intinfovals[2]
  
  yb.val <- intinfovals[3]
  
  l.edgeval <- yb.val - ya.val
  
  exp.val <- exp(-0.5 * ((s.val/bw)^2))
  
  denom.val <- 2 * sqrt(2 * pi) * nrow(edges) * bw * l.edgeval
  
  erfs.val <- erf(yb.val/(bw * sqrt(2))) - erf(ya.val/(bw * sqrt(2)))
  
  ansval.ef <- exp.val * erfs.val / denom.val
  
  return(ansval.ef)
  
}

efnkde <- function(x, y, bw, edges, points) {
  
  ans <- 0
  
  for (i in 1:nrow(edges)) {
    
    ans <- ans + partial.efnkde(x = x, y = y, bw = bw, points = points, edges = edges, iedgeef = i)
    
  }
  
  return(ans)
  
}