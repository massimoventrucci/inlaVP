
# this function define the main effect model, returning the structure matrix 'R'
# it is used to define the main effects, e.g. f1=m(x1, "rw2")
# NOTE: x values must be increasing (THOUGH they can be unequally spaced...)


# x: covariate
# igmrf.type: the main effect model, e.g. "rw2", "besag"
# R: optional; the structure matrix defined through an inla obj 'Cmatrix'
# g: for 'besag' only; the structure matrix defined through an inla obj 'graph'

m <- function(x, igmrf.type = "rw1", R=NULL, g,
              p.spline=FALSE, n.bsplines=20, deg.bsplines=3){

  ###########  ########### ########### utilities
  define.id <- function(x){
    n <- length(x)
    x0 <- unique(x)
    id.x0 <- rank(x0)#order(x0)
    # find the duplicated and assign the rank to them
    id <- rep(NA, n)
    id.dupl <- duplicated(x)
    if (sum(id.dupl)>0){
      id[!id.dupl] <- rank(x[!id.dupl])# order(x[!id.dupl])
      x.dupl <- x[duplicated(x)]
      dictionary.id <- cbind(x0, id.x0)
      for( j in 1:length(x.dupl)) {
        val <- x.dupl[j]
        id.tmp <- which(val == dictionary.id[,1])
        id[which(id.dupl)[j]] <- as.numeric(dictionary.id[id.tmp,2])#as.numeric(dictionary.id[which(x0==val),2])
      }
    } else {
      id <- id.x0
    }
    id
  }

  build.A <- function(id, n, p) {
    A <- matrix(0, nrow=n, ncol=p)
    for(i in 1:n){
      A[i,id[i]] <- 1
    }
    return(A)
  }

  Qrw2 <- function(x) {
    xx <-  sort(unique(x))
    n <- length(xx)
    if(all(diff(xx)==1))
    {
      Qrw <- INLA:::inla.rw(n, order = 2, scale.model=TRUE)
    }else{
      # structure rw2 for irregularly locations
      Qrw <- INLA:::inla.extract.Q("xx", y ~ f(xx, model="rw2",scale.model=TRUE,
                                               initial=0, diagonal = 0, fixed=T),
                                   data = data.frame(xx, y = rep(0, length(xx))))
    }
    return(Qrw)
  }

  Qrw1 <- function(x) {
    xx <-  sort(unique(x))
    n <- length(xx)
    if(all(diff(xx)==1))
    {
      Qrw <- INLA:::inla.rw(n, order = 1, scale.model=TRUE)
    }else{
      # structure rw1 for irregularly locations
      Qrw <- INLA:::inla.extract.Q("xx", y ~ f(xx, model="rw1",scale.model=TRUE,
                                               initial=0, diagonal = 0, fixed=T),
                                   data = data.frame(xx, y = rep(0, length(xx))))
    }
    return(Qrw)
  }

  bspline <- function(x, xl=min(x), xr=max(x), ndx, bdeg) {
    dx <- (xr - xl) / ndx
    knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
    B <- splines::spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
    B
  }

  ###########  ########### ########### END utilities


  if (p.spline & igmrf.type == "besag") {
    stop("p-spline is not available for igmrf.type = 'besag'")
  } else if (p.spline & igmrf.type != "besag") {
    # P-spline case; only for rw1 and rw2 not besag
    A <- bspline(x=x, ndx=n.bsplines-3, bdeg=deg.bsplines)
    if (igmrf.type == "rw1") {
      if (is.null(R)){
        idx <- 1:ncol(A)
        R <- Qrw1(idx)
        rankdef <- 1
      } else {
        idx <- 1:ncol(A)
        rankdef <- 1
      }
    } else if (igmrf.type == "rw2"){
      if (is.null(R)){
        idx <- 1:ncol(A)
        R <- Qrw2(idx)
        rankdef <- 2
      } else {
        idx <- 1:ncol(A)
        rankdef <- 2
      }
    }
  } else {
    # SPDE RW case in Lindgren and Rue paper (approx equivalent to thin-plate spline for rw2)
    idx <- define.id(x)
    A <- build.A(id=idx, n=length(x), p=length(unique(x)))
    if (igmrf.type == "rw1") {
      if (is.null(R)){
        R <- Qrw1(x)
        rankdef <- 1
      } else {
        rankdef <- 1
      }
    } else if (igmrf.type == "rw2"){
      if (is.null(R)){
        R <- Qrw2(x)
        rankdef <- 2
      } else {
        rankdef <- 2
      }
    } else if (igmrf.type == "besag") {
      if (is.null(R)){
        graph.mat <- INLA::inla.as.sparse(INLA::inla.graph2matrix(g))
        R.tmp <- -graph.mat
        diag(R.tmp) <- apply(graph.mat,1,sum)-1
        R <- INLA::inla.scale.model(R.tmp,
                                    constr = list(
                                      A = matrix(1, 1, nrow(graph.mat)),
                                      e=0))
        rankdef <- 1
      } else {
        rankdef <- 1
      }
    }
  }

  # output list
  return(list(x=sort(unique(x)),
              R=R,
              rankdef=rankdef,
              id = idx,
              A = A))
}

control.interaction <- function(m1, m2, interaction.type=4){
  R1 <- m1$R
  R2 <- m2$R
  n1 <- ncol(R1)
  n2 <- ncol(R2)
  A1 <- m1$A
  A2 <- m2$A
  A.int <- INLA::inla.row.kron(M1=A2, M2=A1)
  Identity <- Matrix:::Diagonal(n1*n2, x = 1)
  Id1 = Matrix:::Diagonal(n1, x = 1)
  Id2 = Matrix:::Diagonal(n2, x = 1)
  Kbeta1 <- kronecker(Matrix:::Matrix(1,nrow=n2),Matrix:::Diagonal(n1))
  Kbeta2 <- kronecker(Matrix:::Diagonal(n2), Matrix:::Matrix(1,nrow=n1))
  f2.inter.const <- kronecker(matrix(1,ncol=n2),diag(n1))
  f1.inter.const <- kronecker(diag(n2),matrix(1,ncol=n1))
  # f2.inter.const <- kronecker(matrix(1,ncol=n2), Matrix:::Diagonal(n1))
  # f1.inter.const <- kronecker(Matrix:::Diagonal(n2),matrix(1,ncol=n1))
  if (interaction.type==4){
    Rkron <- kronecker(R2, R1)
    constr <- list(A=rbind(f1.inter.const,
                           f2.inter.const[-1,]),
                   e=matrix(0, n1+n2-1,1))
  } else if (interaction.type==3){  # spatial trend changes at each time
    Rkron <- kronecker(R2, Matrix:::Diagonal(n1))
    constr <- list(A=rbind(f2.inter.const),
                             e=matrix(0, n1,1))

  } else if (interaction.type==2){  # time trend changes at each space
    Rkron <- kronecker(Matrix:::Diagonal(n2), R1)
    constr <- list(A=rbind(f1.inter.const),
                             e=matrix(0, n2,1))
  } else{
    Rkron <- kronecker(Matrix:::Diagonal(n2), Matrix:::Diagonal(n1))
    constr <- NULL
  }
  return(list(n1=n1, n2=n2, R1=R1, R2=R2, Rkron=Rkron, constr=constr,
              A1=A1, A2=A2, A.int=A.int, n=nrow(A.int),
              Identity=Identity, Id1=Id1, Id2=Id2, Kbeta1=Kbeta1, Kbeta2=Kbeta2))
}
