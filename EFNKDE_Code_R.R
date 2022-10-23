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

FindIntrsctnPtGivenLines(a1,b1,c1,a2,b2,c2) {
  return(c(b2*c1 - b1*c2)/(a2*b1 - a1*b2),
         (a1*c2 - a2*c1)/(a2*b1 - b1*b2))
}

FindIfPtInt <- function(PtA, PtB, PtX) {
  xA <- ptA[0]
  xB <- 
}

