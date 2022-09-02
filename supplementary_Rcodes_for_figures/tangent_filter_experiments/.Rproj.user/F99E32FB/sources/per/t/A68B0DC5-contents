foo <- function(gg_data, legendPosition = "none"){
  ggplot(data = gg_data, aes(x = Time, y = Estimate, col = Start)) + 
    geom_line(alpha = 0.5) + scale_color_viridis_d() +
    geom_hline(yintercept = trueTheta, col = "red") + 
    scale_x_continuous(expand = c(0, 100)) + scale_y_continuous(expand = c(0, 0)) +
    theme_bw() + 
    theme(legend.position = legendPosition,
          axis.title = element_text(size = title_size),
          axis.text = element_text(size = axis_size))
}
smoothEst <- function(vecEst, n0 = 500){
  n <- length(vecEst)
  if(n0 >= n)
    n0 <- max(1, floor(n / 10))
  output <- vecEst
  tailEst <- vecEst[(n0 + 1) : n]
  output[(n0 + 1) : n] <- cumsum(tailEst) * 1 / ((n0 + 1) : n - n0)
  output
}

foo_legend <- function(name = "oneObs_oneStart_severalGrads"){
  thetaEst <- read.table(paste0("simulation_results/est_", name, ".txt"),
                         h = F, sep = ";")
  thetaEst <- thetaEst[, order(thetaEst[1, ])]
  theta_gg <- do.call(rbind.data.frame,
                      lapply(1:ncol(thetaEst), function(i){
                        data.frame(Estimate = thetaEst[, i], Time = times,
                                   Start = factor(rep(i, nrow(thetaEst)),
                                                  levels = 1:ncol(thetaEst)))
                      }))
  thetaMoy <- apply(thetaEst, 2, smoothEst)
  theta_gg_smooth <- do.call(rbind.data.frame,
                             lapply(1:ncol(thetaMoy), function(i){
                               data.frame(Estimate = thetaMoy[, i], Time = times,
                                          Start = factor(rep(i, nrow(thetaEst)),
                                                         levels = 1:ncol(thetaMoy)))
                             }))
  foo_tmp <- function(gg_data){
    ggplot(data = gg_data, aes(x = Time, y = Estimate, col = Start)) + 
      geom_line(alpha = 0.5) + scale_color_viridis_d(labels=paste(seq(0.5, 1, length = 6))) +
      geom_hline(yintercept = trueTheta, col = "red") + 
      scale_x_continuous(expand = c(0, 100)) + scale_y_continuous(expand = c(0, 0))+
      guides(colour = guide_legend(title = expression(kappa))) +
      theme_bw() + 
      theme(legend.position = c(0.75, 0.75), legend.title = element_text(size = title_size),
            legend.text = element_text(size = axis_size),
            axis.title = element_text(size = title_size), legend.key = element_blank(),
            legend.background = element_blank(),
            axis.text = element_text(size = axis_size))
  }
  
  savePdf(name)
  print(foo_tmp(theta_gg))
  dev.off()
  embedPdf(name)
  name_smooth <- paste(name, "smooth", sep = "_")
  savePdf(name_smooth)
  print(foo_tmp(theta_gg_smooth))
  dev.off()
  embedPdf(name_smooth)
}
