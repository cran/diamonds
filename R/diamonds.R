#  diamonds.R
#  d white, 9 july 2001, from Splus version of 23 april 1998
#
#  functions for creating graphics for hierarchical
#  diamond partitions, random points and walks, 
#  dual hex lines, center hex lines

#  set up fixed aspect plot space
diamond.plot <- function (cl)
  {
    rx <- range (cl$x[!is.na(cl$y)], na.rm=TRUE)
    ry <- range (cl$y[!is.na(cl$x)], na.rm=TRUE)
    dx <- diff (rx); dy <- diff (ry);
    px <- par ("pin")[1]; py <- par ("pin")[2]
    if (dx/dy > px/py) py <- px * dy/dx else px <- py * dx/dy
    par (pin = c(px, py))
    plot (rx, ry, type="n", axes=FALSE, xlab="", ylab="")
  }

#  vertices of portrait format diamond with center at (x, y)
#  and short axis a
diamond.base <- function (x=0, y=0, a=1)
  {
    list (
    c(x-a/2, y),
    c(x,     y-a*sqrt(3)/2),
    c(x+a/2, y),
    c(x,     y+a*sqrt(3)/2))
  }

#  make x and y vectors of centers of diamonds
diamond.centers <- function (b, d)
  {
    if (d < 2) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      t <- tl[[i]]
      u <- (t[[1]][1] + t[[3]][1]) / 2 
      v <- (t[[2]][2] + t[[4]][2]) / 2
      x <- c(x, u)
      y <- c(y, v) }
    list (x=x, y=y)
  }

#  make x and y vectors of centers of duals of diamonds
diamond.dualcents <- function (b, d)
  {
    if (d < 2) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      t <- tl[[i]]
      u <- t[[1]][1]
      v <- (t[[2]][2] + t[[4]][2]) / 2
      x <- c(x, u)
      y <- c(y, v) }
    list (x=x, y=y)
  }

#  make x and y vectors for dual hexes of diamonds
diamond.dualedges <- function (b, d)
  {
    if (d < 1) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      t <- tl[[i]]
      p1 <- c( (t[[1]][1] + t[[2]][1]) / 2, 
               (t[[1]][2] + t[[2]][2]) / 2)
      p2 <- c( (t[[2]][1] + t[[3]][1]) / 2, 
               (t[[2]][2] + t[[3]][2]) / 2)
      p3 <- c( (t[[3]][1] + t[[4]][1]) / 2, 
               (t[[3]][2] + t[[4]][2]) / 2)
      p4 <- c( (t[[1]][1] + t[[4]][1]) / 2, 
               (t[[1]][2] + t[[4]][2]) / 2)
      p5 <- c(  t[[4]][1],
               (2*t[[1]][2] + 1*t[[4]][2]) / 3)
      p6 <- c(  t[[4]][1],
               (2*t[[1]][2] + 1*t[[2]][2]) / 3)
      x <- c(x,p1[1],p6[1],NA,p2[1],p6[1],NA,p3[1],p5[1],NA,
               p4[1],p5[1],NA,p5[1],p6[1],NA)
      y <- c(y,p1[2],p6[2],NA,p2[2],p6[2],NA,p3[2],p5[2],NA,
               p4[2],p5[2],NA,p5[2],p6[2],NA) }
    list (x=x, y=y)
  }

#  make x and y vectors of edges of diamonds
diamond.edges <- function (b, d)
  {
    if (d < 1) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      w <- unlist (c (tl[[i]], tl[[i]][1]))
      u <- split (w, rep (c("x","y"), length(w)/2))
      x <- c(x, u$x, NA)
      y <- c(y, u$y, NA) }
    list (x=x, y=y)
  }

#  walk on hierarchical diamond structure
diamond.hierwalk <- function (b, d)
  {
    if (d < 2) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      t <- tl[[i]]
      u <- (t[[1]][1] + t[[3]][1]) / 2 
      v <- (t[[2]][2] + t[[4]][2]) / 2
      x <- c(x, u); y <- c(y, v) }
    x <- c(x, NA); y <- c(y, NA);
    list (x=x, y=y)
  }

#  label diamonds at depth d by 4-fold subdivision
diamond.labels <- function (d)
  {
    if (d < 2) return (NULL)
    dt <- list ()
    unlist (labeldeep (dt, "", 2, d))
  }

#  random permutation of labels on hierarchical diamonds
diamond.randlabels <- function (d)
  {
    if (d < 2) return (NULL)
    tl <- hierlabels ("", 1, d-1)
    tl <- randlevel (tl, 2, d)
    unlist (tl)
  }

#  make x and y vectors of random points in diamonds
diamond.randpts <- function (b, d)
  {
    if (d < 2) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      u <- ranptdia (tl[[i]])
      x <- c(x, u[1])
      y <- c(y, u[2]) }
    list (x=x, y=y)
  }

#  random walk on hierarchical diamond structure
diamond.randwalk <- function (b, d)
  {
    if (d < 2) return (NULL)
    tl <- hiercenters (b, d)
    tl <- randlevel (tl, 2, d)
    x <- NULL; y <- NULL;
    for (i in 1:length (tl)) {
      x <- c(x, tl[[i]][1]); y <- c(y, tl[[i]][2]) }
    x <- c(x, NA); y <- c(y, NA);
    list (x=x, y=y)
  }

#  make x and y vectors for short axis of diamonds
#  (completing triangles when combined with edges)
diamond.triedges <- function (b, d)
  {
    if (d < 2) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      t <- tl[[i]]
      x <- c(x,t[[1]][1],t[[3]][1],NA)
      y <- c(y,t[[1]][2],t[[3]][2],NA) }
    list (x=x, y=y)
  }

#  partition diamond at depth d by 4-fold subdivision
diamonds <- function (b, d)
  {
    if (d < 2) return (list (b))
    deep (b, 2, d)
  }

#  make x and y vectors for centered hexes of diamonds in list
hexlines <- function (b, d)
  {
    if (d < 2) return (NULL)
    tl <- diamonds (b, d)
    x <- NULL; y <- NULL;
    if (! is.list (tl[[1]])) tl <- list (tl)
    for (i in 1:length (tl)) {
      t <- tl[[i]]
      p1 <- c(   t[[1]][1], 
              (2*t[[1]][2] + 1*t[[2]][2]) / 3)
      p2 <- c(   t[[2]][1], 
              (1*t[[1]][2] + 2*t[[2]][2]) / 3)
      p3 <- c(   t[[3]][1], 
              (2*t[[1]][2] + 1*t[[2]][2]) / 3)
      p4 <- c(   t[[3]][1], 
              (1*t[[4]][2] + 2*t[[1]][2]) / 3)
      p5 <- c(   t[[4]][1],
              (2*t[[4]][2] + 1*t[[1]][2]) / 3)
      p6 <- c(   t[[1]][1],
              (1*t[[4]][2] + 2*t[[1]][2]) / 3)
      x <- c(x,p1[1],p2[1],p3[1],p4[1],p5[1],p6[1],p1[1],NA)
      y <- c(y,p1[2],p2[2],p3[2],p4[2],p5[2],p6[2],p1[2],NA) }
    list (x=x, y=y)
  }

#  recurse on subdiamonds
deep <- function (b, level, d)
  {
    st <- subdiamonds (b)
    if (level == d) return (st)
    sl <- list ()
    for (i in 1:4) sl <- c (sl, deep (st[[i]], level+1, d))
    sl
  }

#  recurse for labeling
labeldeep <- function (dt, parent, level, d)
  {
    st <- c(paste(parent,"0",sep=""),paste(parent,"1",sep=""),
            paste(parent,"2",sep=""),paste(parent,"3",sep=""))
    if (level == d) return ( list(unlist(c(dt, st))) )
    for (i in 1:4) dt <- labeldeep (dt, st[[i]], level+1, d)
    dt
  }

#  recurse on subcenters
centerdeep <- function (b, level, d)
  {
    if (level == d) return (subcenters (b))
    st <- subdiamonds (b)
    sl <- list ()
    for (i in 1:length(st)) 
      sl <- c (sl, list (centerdeep (st[[i]], level+1, d)))
    sl
  }

#  randomly permute list elements and recurse
#  and flatten list to one level
randlevel <- function (tl, level, d)
  {
    temp <- list ()
    perm <- sample (1:4)
    for (i in 1:4) temp[[i]] <- tl[[perm[i]]]
    for (i in 1:4) tl[[i]] <- temp[[i]]
    if (level == d) return (tl)
    sl <- list ()
    for (i in 1:4) 
      sl <- c(sl, randlevel (tl[[i]], level+1, d))
    sl
  }

#  generate random point in diamond
ranptdia <- function (b)
  {
    x <- runif (1, min=0, max=1)
    y <- runif (1, min=0, max=sqrt(3)/2)
    if ((x < 0.5) && (y > x*sqrt(3))) {
      x <- x + 0.5; y <- y - sqrt(3)/2 }
    if ((x > 0.5) && (y > (1-x)*sqrt(3))) {
      x <- x - 0.5; y <- y - sqrt(3)/2 }
    xmin <- b[[1]][1]; xmax <- b[[3]][1];
    ymin <- b[[2]][2]; ymax <- b[[4]][2];
    x <- x*(xmax - xmin) + xmin
    y <- (y + sqrt(3)/2)*(ymax - ymin)/sqrt(3) + ymin
    c(x, y)
  }

#  make diamond centers at depth d by 4-fold subdivision
hiercenters <- function (b, d)
  {
    if (d < 2) 
      return (list (c((b[[1]][1]+b[[3]][1])/2,(b[[2]][2]+b[[4]][2])/2)))
    centerdeep (b, 2, d)
  }

#  recurse 
hierlabels <- function (parent, level, d)
  {
    st <- list(paste(parent,"0",sep=""),paste(parent,"1",sep=""),
               paste(parent,"2",sep=""),paste(parent,"3",sep=""))
    if (level == d) return (st)
    sl <- list ()
    for (i in 1:length(st)) 
      sl <- c (sl, list (hierlabels (st[[i]], level+1, d)))
    sl
  }

#  partition diamond by 4-fold subdivision
subdiamonds <- function (b)
  {
    p1 <- c( (b[[1]][1] + b[[2]][1]) / 2, 
             (b[[1]][2] + b[[2]][2]) / 2)
    p2 <- c( (b[[2]][1] + b[[3]][1]) / 2, 
             (b[[2]][2] + b[[3]][2]) / 2)
    p3 <- c( (b[[3]][1] + b[[4]][1]) / 2, 
             (b[[3]][2] + b[[4]][2]) / 2)
    p4 <- c( (b[[1]][1] + b[[4]][1]) / 2, 
             (b[[1]][2] + b[[4]][2]) / 2)
    p5 <- c( (b[[1]][1] + b[[3]][1]) / 2, 
             (b[[2]][2] + b[[4]][2]) / 2)
    list (list(b[[1]], p1, p5, p4), list(p1, b[[2]], p2, p5),
          list(p4, p5, p3, b[[4]]), list(p5, p2, b[[3]], p3))
  }

#  centers of partition on diamond
subcenters <- function (b) 
  {
    c1 <- c((b[[1]][1]+b[[2]][1])/2,  b[[1]][2])
    c2 <- c( b[[2]][1],              (b[[1]][2]+b[[2]][2])/2)
    c3 <- c( b[[2]][1],              (b[[1]][2]+b[[4]][2])/2)
    c4 <- c((b[[2]][1]+b[[3]][1])/2,  b[[1]][2])
    list (c1, c2, c3, c4)
  }
