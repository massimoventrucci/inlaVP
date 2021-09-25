

## VP table and VP plots

# table.type:
# 'mixing' = posterior summaries for the mixing hyperparameters
#            NOTE: gives the estimated contribution of: main vs int; space vs time; time:iid vs time str; space:iid vs space.str)

# 'explained_sd_v1' = we apply the posterior estimates of the mixing hyperpars to the posterior mean of sqrt(1/tau) (i.e. total st dev)
#                     to compute the portion of total st dev attributable to the different model components:
#                     main, int, main:space, main:time, time:iid, time:str, space:iid, space:str
#                     NOTE: this gives the estimated portion of total st dev attributable to:
#                     main, int, main:space, main:time; time:iid, time str; space:iid, space.str
#                     but, row should not be read independently,
#                     e.g. the st dev for 'main:space' is the portion of st dev contributed by 'main' (not the total st dev), that can be attributed to spatial variation

# 'explained_sd_v2' = the posterior mean of the st dev of each model component:
#                     total, main, int, main:space, main:time, time:iid, time:str, space:iid, space:str
#                     NOTE: this gives the estimated portion of total st dev attributable to:
#                     main, int, main:space, main:time; time:iid, time str; space:iid, space.str
#                     but this time each row can be read independently I think...

# 'explained_sd_prop' = like 'explained_sd_v1' but rescaling (dividing for the total st.dev)
#                       to have variance partitioning in the proportion scale (0,1)


vp.table <- function(inla.res,
                     vp.model='striid',
                     table.type = 'mixing',
                     nsim=10000){
  res.hypersample <- INLA:::inla.hyperpar.sample(nsim, inla.res, intern=TRUE)
  if(vp.model=='str'){

    stop('This function only works for model (2); model (1) not done yet')

  } else if(vp.model=='striid') {

    gamma.post <- apply(res.hypersample, 1,
                        FUN=function(x) inlaVP:::theta.to.gamma.striid(x))
    phi.post <- apply(res.hypersample, 1,
                      FUN=function(x) inlaVP:::theta.to.phi.striid(x))
    psi1.post <- apply(res.hypersample, 1,
                       FUN=function(x) inlaVP:::theta.to.psi1.striid(x))
    psi2.post <- apply(res.hypersample, 1,
                       FUN=function(x) inlaVP:::theta.to.psi2.striid(x))
    totsd.post <- apply(res.hypersample, 1,
                        FUN=function(x) 1/sqrt(inlaVP:::theta.to.tau.striid(x)))

    if (table.type=='mixing') {
      res.table <- data.frame(rbind(c(mean(1-gamma.post),
                                      quantile(1-gamma.post, 0.025),
                                      quantile(1-gamma.post, 0.25),
                                      quantile(1-gamma.post, 0.5),
                                      quantile(1-gamma.post, 0.75),
                                      quantile(1-gamma.post, 0.975)),
                                    c(mean(gamma.post),
                                      quantile(gamma.post, 0.025),
                                      quantile(gamma.post, 0.25),
                                      quantile(gamma.post, 0.5),
                                      quantile(gamma.post, 0.75),
                                      quantile(gamma.post, 0.975)),
                                    c(mean(phi.post),
                                      quantile(phi.post, 0.025),
                                      quantile(phi.post, 0.25),
                                      quantile(phi.post, 0.5),
                                      quantile(phi.post, 0.75),
                                      quantile(phi.post, 0.975)),
                                    c(mean(1-phi.post),
                                      quantile(1-phi.post, 0.025),
                                      quantile(1-phi.post, 0.25),
                                      quantile(1-phi.post, 0.5),
                                      quantile(1-phi.post, 0.75),
                                      quantile(1-phi.post, 0.975)),
                                    c(mean(psi1.post),
                                      quantile(psi1.post, 0.025),
                                      quantile(psi1.post, 0.25),
                                      quantile(psi1.post, 0.5),
                                      quantile(psi1.post, 0.75),
                                      quantile(psi1.post, 0.975)),
                                    c(mean(1-psi1.post),
                                      quantile(1-psi1.post, 0.025),
                                      quantile(1-psi1.post, 0.25),
                                      quantile(1-psi1.post, 0.5),
                                      quantile(1-psi1.post, 0.75),
                                      quantile(1-psi1.post, 0.975)),
                                    c(mean(psi2.post),
                                      quantile(psi2.post, 0.025),
                                      quantile(psi2.post, 0.25),
                                      quantile(psi2.post, 0.5),
                                      quantile(psi2.post, 0.75),
                                      quantile(psi2.post, 0.975)),
                                    c(mean(1-psi2.post),
                                      quantile(1-psi2.post, 0.025),
                                      quantile(1-psi2.post, 0.25),
                                      quantile(1-psi2.post, 0.5),
                                      quantile(1-psi2.post, 0.75),
                                      quantile(1-psi2.post, 0.975))))
      rownames(res.table) <- c('main', 'int',
                               'main.space', 'main.time',
                               'main.time.iid', 'main.time.str',
                               'main.space.iid', 'main.space.str')
      colnames(res.table) <- c('mean', 'q2.5', 'q25', 'q50', 'q75', 'q97.5')
    } else if(table.type=='explained_sd_v1'){
      res.table <- data.frame(rbind(c(mean(totsd.post), # tot sd
                                      quantile(totsd.post, 0.025),
                                      quantile(totsd.post, 0.25),
                                      quantile(totsd.post, 0.5),
                                      quantile(totsd.post, 0.75),
                                      quantile(totsd.post, 0.975)),
                                    c(mean(totsd.post*(1-gamma.post)), # explained by main (space+time)
                                      quantile(totsd.post*(1-gamma.post), 0.025),
                                      quantile(totsd.post*(1-gamma.post), 0.25),
                                      quantile(totsd.post*(1-gamma.post), 0.5),
                                      quantile(totsd.post*(1-gamma.post), 0.75),
                                      quantile(totsd.post*(1-gamma.post), 0.975)),
                                    c(mean(totsd.post*gamma.post), # tot sd explained by int
                                      quantile(totsd.post*(gamma.post), 0.025),
                                      quantile(totsd.post*(gamma.post), 0.25),
                                      quantile(totsd.post*(gamma.post), 0.5),
                                      quantile(totsd.post*(gamma.post), 0.75),
                                      quantile(totsd.post*(gamma.post), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*phi.post)), # explained by space main effect
                                      quantile(totsd.post*((1-gamma.post)*phi.post), 0.025),
                                      quantile(totsd.post*((1-gamma.post)*phi.post), 0.25),
                                      quantile(totsd.post*((1-gamma.post)*phi.post), 0.5),
                                      quantile(totsd.post*((1-gamma.post)*phi.post), 0.75),
                                      quantile(totsd.post*((1-gamma.post)*phi.post), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*(1-phi.post))), # explained by time main effect
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)), 0.025),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)), 0.25),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)), 0.5),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)), 0.75),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*(1-phi.post)*psi1.post)), # explained by time iid
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*psi1.post), 0.025),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*psi1.post), 0.25),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*psi1.post), 0.5),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*psi1.post), 0.75),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*psi1.post), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*(1-phi.post)*(1-psi1.post))), # explained by time str
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.025),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.25),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.5),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.75),
                                      quantile(totsd.post*((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*phi.post*psi2.post)), # explained by space iid
                                      quantile(totsd.post*((1-gamma.post)*phi.post*psi2.post), 0.025),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*psi2.post), 0.25),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*psi2.post), 0.5),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*psi2.post), 0.75),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*psi2.post), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*phi.post*(1-psi2.post))), # explained by space str
                                      quantile(totsd.post*((1-gamma.post)*phi.post*(1-psi2.post)), 0.025),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*(1-psi2.post)), 0.25),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*(1-psi2.post)), 0.5),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*(1-psi2.post)), 0.75),
                                      quantile(totsd.post*((1-gamma.post)*phi.post*(1-psi2.post)), 0.975))))
      rownames(res.table) <- c('SD', 'SD_main', 'SD_int',
                               'SD_main.space', 'SD_main.time',
                               'SD_main.time.iid', 'SD_main.time.str',
                               'SD_main.space.iid', 'SD_main.space.str')
      colnames(res.table) <- c('mean', 'q2.5', 'q25', 'q50', 'q75', 'q97.5')

    } else if(table.type=='explained_sd_v2'){
      res.table <- data.frame(rbind(c(mean(totsd.post), # tot sd
                                      quantile(totsd.post, 0.025),
                                      quantile(totsd.post, 0.25),
                                      quantile(totsd.post, 0.5),
                                      quantile(totsd.post, 0.75),
                                      quantile(totsd.post, 0.975)),
                                    c(mean(totsd.post*sqrt(1-gamma.post)), # explained by main (space+time)
                                      quantile(totsd.post*sqrt(1-gamma.post), 0.025),
                                      quantile(totsd.post*sqrt(1-gamma.post), 0.25),
                                      quantile(totsd.post*sqrt(1-gamma.post), 0.5),
                                      quantile(totsd.post*sqrt(1-gamma.post), 0.75),
                                      quantile(totsd.post*sqrt(1-gamma.post), 0.975)),
                                    c(mean(totsd.post*gamma.post), # tot sd explained by int
                                      quantile(totsd.post*sqrt(gamma.post), 0.025),
                                      quantile(totsd.post*sqrt(gamma.post), 0.25),
                                      quantile(totsd.post*sqrt(gamma.post), 0.5),
                                      quantile(totsd.post*sqrt(gamma.post), 0.75),
                                      quantile(totsd.post*sqrt(gamma.post), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*phi.post)), # explained by space main effect
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post), 0.025),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post), 0.25),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post), 0.5),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post), 0.75),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post), 0.975)),
                                    c(mean(totsd.post*((1-gamma.post)*(1-phi.post))), # explained by time main effect
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)), 0.025),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)), 0.25),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)), 0.5),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)), 0.75),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)), 0.975)),
                                    c(mean(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*psi1.post)), # explained by time iid
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*psi1.post), 0.025),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*psi1.post), 0.25),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*psi1.post), 0.5),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*psi1.post), 0.75),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*psi1.post), 0.975)),
                                    c(mean(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*(1-psi1.post))), # explained by time str
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.025),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.25),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.5),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.75),
                                      quantile(totsd.post*sqrt((1-gamma.post)*(1-phi.post)*(1-psi1.post)), 0.975)),
                                    c(mean(totsd.post*sqrt((1-gamma.post)*phi.post*psi2.post)), # explained by space iid
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*psi2.post), 0.025),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*psi2.post), 0.25),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*psi2.post), 0.5),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*psi2.post), 0.75),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*psi2.post), 0.975)),
                                    c(mean(totsd.post*sqrt((1-gamma.post)*phi.post*(1-psi2.post))), # explained by space str
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*(1-psi2.post)), 0.025),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*(1-psi2.post)), 0.25),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*(1-psi2.post)), 0.5),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*(1-psi2.post)), 0.75),
                                      quantile(totsd.post*sqrt((1-gamma.post)*phi.post*(1-psi2.post)), 0.975))))
      rownames(res.table) <- c('SD', 'SD_main', 'SD_int',
                               'SD_main.space', 'SD_main.time',
                               'SD_main.time.iid', 'SD_main.time.str',
                               'SD_main.space.iid', 'SD_main.space.str')
      colnames(res.table) <- c('mean', 'q2.5', 'q25', 'q50', 'q75', 'q97.5')

    }
  }
  return(res.table)
}




vp.plotpost <- function(inla.res,
                         vp.model='striid',
                         nsim=10000,
                         ...){
  res.hypersample <- INLA:::inla.hyperpar.sample(nsim, inla.res, intern=TRUE)
  if(vp.model=='str'){
    stop('This function only works for model (2); model (1) not done yet')
  } else if(vp.model=='striid') {
    gamma.post <- apply(res.hypersample, 1,
                        FUN=function(x) inlaVP:::theta.to.gamma.striid(x))
    phi.post <- apply(res.hypersample, 1,
                      FUN=function(x) inlaVP:::theta.to.phi.striid(x))
    psi1.post <- apply(res.hypersample, 1,
                       FUN=function(x) inlaVP:::theta.to.psi1.striid(x))
    psi2.post <- apply(res.hypersample, 1,
                       FUN=function(x) inlaVP:::theta.to.psi2.striid(x))
    totsd.post <- apply(res.hypersample, 1,
                        FUN=function(x) 1/sqrt(inlaVP:::theta.to.tau.striid(x)))
  }
  hist(gamma.post, freq = F, breaks = 30, col=2, ...)
  hist(phi.post, freq = F, breaks = 30, add=TRUE, col=3)
  hist(psi1.post, freq = F, breaks = 30, col=4, add=TRUE)
  hist(psi2.post, freq = F, breaks = 30, col=5, add=TRUE)
  legend('topright', legend=c('gamma', 'phi', 'psi1', 'psi2'),
         col=2:5, lty=1, lwd=3, bty='n')
}



vp.postsample <- function(inla.res,
                          vp.model='striid',
                          nsim=10000,
                          ...){
  res.hypersample <- INLA:::inla.hyperpar.sample(nsim, inla.res, intern=TRUE)
  if(vp.model=='str'){
    stop('This function only works for model (2); model (1) not done yet')
  } else if(vp.model=='striid') {
    gamma.post <- apply(res.hypersample, 1,
                        FUN=function(x) inlaVP:::theta.to.gamma.striid(x))
    phi.post <- apply(res.hypersample, 1,
                      FUN=function(x) inlaVP:::theta.to.phi.striid(x))
    psi1.post <- apply(res.hypersample, 1,
                       FUN=function(x) inlaVP:::theta.to.psi1.striid(x))
    psi2.post <- apply(res.hypersample, 1,
                       FUN=function(x) inlaVP:::theta.to.psi2.striid(x))
    totsd.post <- apply(res.hypersample, 1,
                        FUN=function(x) 1/sqrt(inlaVP:::theta.to.tau.striid(x)))
  }
  return(data.frame(gamma=gamma.post,
                    phi=phi.post,
                    psi1=psi1.post,
                    psi2=psi2.post,
                    totsd=totsd.post))
}


vp.plot <- function(inla.res,
                    vp.model='striid',
                    table.type='mixing',
                    cex.sources=0.95,
                    title.plot='',
                    nsim=10000, ...){
  if(vp.model=='str'){
    stop('This function only works for model (2); model (1) not done yet')
  } else if(vp.model=='striid') {

    dfrx <- vp.table(inla.res=inla.res, table.type=table.type)

    if(table.type=='mixing'){
      # 'mixing'
      label.sources <- c('SOURCE',
                         'main', 'int',
                         'main:space','main:time',
                         'time:iid','time:str',
                         'space:iid','space:str')
      par(mar=c(5,5,5,2)+0.1)
      topleft.vert <- 13
      id.sources <- c(1,2,4,5,7,8,10,11, 13) # seq(2,16,2)
      plot(c(0,1),
           c(0,topleft.vert),
           type='n',
           xlab='', ylab='',
           axes=F, ...)
      axis(1, at=seq(0,1,by=0.1), labels=seq(0,1,by=0.1), lwd.ticks=0.5)
      axis(2, at=id.sources, labels = rev(label.sources),
           tick=F,
           lwd.ticks=0.5,
           cex.axis=cex.sources,
           #hadj= 1, padj = 1,
           las=2)
      axis(3, at=seq(0,1,by=0.1), labels=seq(0,1,by=0.1), lwd.ticks=0.5)
      dfrx.tmp <- dfrx[8:1,]
      for(i in 1:8){
        segments(x0=dfrx.tmp[i,'q2.5'], y0=id.sources[i],
                 x1=dfrx.tmp[i,'q97.5'], y1=id.sources[i], lwd=1)
        segments(x0=dfrx.tmp[i,'q25'], y0=id.sources[i],
                 x1=dfrx.tmp[i,'q75'], y1=id.sources[i], lwd=2.5)
        points(x=dfrx.tmp[i,'q50'], y=id.sources[i], pch=16, cex=1)
      }
      segments(x0=0, y0=-.5, x1=0, y1=topleft.vert+.5)
      segments(x0=1, y0=-.5, x1=1, y1=topleft.vert+.5)
      title(main=title.plot,
            #adj=0,
            line=3)

    } else if(table.type=='explained_sd_v1' | table.type=='explained_sd_v2'){
      # 'explained_sd'
      label.sources <- c('SOURCE',
                         'total sd', 'main', 'int',
                         'main:space','main:time',
                         'time:iid','time:str',
                         'space:iid','space:str')
      par(mar=c(5,5.5,5,2)+0.1)
      topleft.vert <- 15
      right.bound <- max(dfrx)+0.01
      id.sources <- c(1,2,4,5,7,8,10,11,12, 15) # seq(2,16,2)
      plot(c(0,right.bound),
           c(0,topleft.vert),
           type='n',
           xlab='', ylab='',
           axes=F,
           ...)
      axis(1, at=round(seq(0,right.bound, by = right.bound/5),2),
           labels=round(seq(0,right.bound, by = right.bound/5),2),
           lwd.ticks=0.5)
      axis(2, at=id.sources, labels = rev(label.sources),
           tick=F,
           lwd.ticks=0.5,
           cex.axis=cex.sources,
           #hadj= 1, padj = 1,
           las=2)
      axis(3, at=round(seq(0,right.bound, by = right.bound/5),2),
           labels=round(seq(0,right.bound, by = right.bound/5),2),
           lwd.ticks=0.5)
      dfrx.tmp <- dfrx[9:1,]
      for(i in 1:9){
        segments(x0=dfrx.tmp[i,'q2.5'], y0=id.sources[i],
                 x1=dfrx.tmp[i,'q97.5'], y1=id.sources[i], lwd=1)
        segments(x0=dfrx.tmp[i,'q25'], y0=id.sources[i],
                 x1=dfrx.tmp[i,'q75'], y1=id.sources[i], lwd=2.5)
        points(x=dfrx.tmp[i,'q50'], y=id.sources[i], pch=16, cex=1)
      }
      segments(x0=0, y0=-.5, x1=0, y1=topleft.vert+.5)
      segments(x0=right.bound, y0=-.5, x1=right.bound, y1=topleft.vert+.5)
      title(main=title.plot,
            #adj=0,
            line=3)
    }
  }
}




# vp.plot2 <- function(inla.res,
#                      vp.model='striid',
#                      table.type='explained_sd_v1',
#                      cex.sources=0.95,
#                      eps.right=0.01,  #control right boundary of the plot
#                      right.bound.sd=NULL,
#                      title.plot='',
#                      nsim=10000, ...){
#   if(vp.model=='str'){
#     stop('This function only works for model (2); model (1) not done yet')
#   } else if(vp.model=='striid') {
#
#     dfrx <- vp.table(inla.res=inla.res, table.type=table.type)
#
#     # 'explained_sd'
#     label.sources <- c('SOURCE',
#                        'tot st_dev:', 'main', 'int',
#                        'main:', 'space','time',
#                        'time:', 'iid','str',
#                        'space:', 'iid','str')
#     par(mar=c(5,5,5,2)+0.1)
#     topleft.vert <- 15#18
#     bottomleft.vert <- 0
#     right.bound <- NULL
#     right.bound <- ifelse(is.null(right.bound.sd),
#                           max(dfrx)+eps.right, right.bound.sd)
#     id.sources <- c(1,2,3,  5,6,7,  9,10,11,   13,14,15, 16) - 1#, 18)
#     id.breaks <- c(4,8,12)-1
#     plot(c(bottomleft.vert,right.bound),
#          c(bottomleft.vert,topleft.vert),
#          type='n',
#          xlab='', ylab='',
#          axes=F,
#          ...)
#     axis(1, at=round(seq(0,right.bound, by = right.bound/5),2),
#          labels=round(seq(0,right.bound, by = right.bound/5),2),
#          lwd.ticks=0.5)
#     axis(2, at=id.sources, labels = rev(label.sources),
#          tick=F,
#          lwd.ticks=0.5,
#          cex.axis=cex.sources,
#          #hadj= 1, padj = 1,
#          las=2)
#     axis(3, at=round(seq(0,right.bound, by = right.bound/5),2),
#          labels=round(seq(0,right.bound, by = right.bound/5),2),
#          lwd.ticks=0.5)
#     dfrx.tmp <- rbind(dfrx,
#                       dfrx[c('SD_main','SD_main.space','SD_main.time'),])
#     dfrx.tmp <- dfrx.tmp[c('SD_main.space.str',
#                            'SD_main.space.iid',
#                            'SD_main.space1',
#                            'SD_main.time.str',
#                            'SD_main.time.iid',
#                            'SD_main.time1',
#                            'SD_main.time',
#                            'SD_main.space',
#                            'SD_main1',
#                            'SD_int',
#                            'SD_main',
#                            'SD'),]
#     for(i in 1:12){
#       segments(x0=dfrx.tmp[i,'q2.5'], y0=id.sources[i],
#                x1=dfrx.tmp[i,'q97.5'], y1=id.sources[i], lwd=1)
#       segments(x0=dfrx.tmp[i,'q25'], y0=id.sources[i],
#                x1=dfrx.tmp[i,'q75'], y1=id.sources[i], lwd=2.5)
#       points(x=dfrx.tmp[i,'q50'], y=id.sources[i], pch=16, cex=1)
#     }
#     segments(x0=0, y0=-.5, x1=0, y1=topleft.vert+.5)
#     segments(x0=right.bound, y0=-.5, x1=right.bound, y1=topleft.vert+.5)
#     ## line breaks
#     segments(x0=c(0,0,0), y0=id.breaks,
#              x1=rep(right.bound,3), y1=id.breaks, lty=3, col='black')
#     title(main=title.plot,
#           #adj=0,
#           line=3)
#   }
# }
#
#
