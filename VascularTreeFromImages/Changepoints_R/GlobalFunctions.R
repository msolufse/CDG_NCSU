LargeExcel = function(fname){
  #getting info about all the excel sheets
  sheets = readxl::excel_sheets(fname)
  tibble = lapply(sheets,function(x) readxl:: read_excel(fname,sheet=x))
  data_frame = lapply(tibble,as.data.frame)
  
  #assigning names to data frames
  names(data_frame) = sheets
  print(data_frame)
}


absmin = function(x){
  x = abs(x)
  min(x[x>0], na.rm = TRUE)
}


closest = function(x,y){
  which(abs(x-y)==min(abs(x-y)))
}


FindChangepoints = function(data, q, plot,name){suppressWarnings({
  for (i in 1:length(data)){
    r = data[[i]]$Radius
    l = data[[i]]$Length
    mod = length(r)-2
    rmod = r[1:mod]
    lmod = l[1:mod]
    
    # detect number of change points
    if (length(l)>10){
      cpt = cpt.meanvar(r,penalty = "BIC", method = 'BinSeg', Q=q)
      
      #build regression model
      fit_lm = lm(r~l)
      t = try(segmented(fit_lm, seg.Z = ~l, npsi = length(cpt@cpts)),silent = T)
      if("try-error" %in% class(t)){
        fit_lm = lm(rmod~lmod)
        t2 = try(segmented(fit_lm, seg.Z = ~lmod, npsi = length(cpt@cpts)))
        if("try-error" %in% class(t2)){
          mod2 = mod-4
          rmod2 = r[1:mod2]
          lmod2 = l[1:mod2]
          fit_lm = lm(rmod2~lmod2)
          t3 = try(segmented(fit_lm, seg.Z = ~lmod2, npsi = length(cpt@cpts)))
          if("try-error" %in% class(t3)){
            fit_segmented = fit_lm
            slopes = coef(fit_segmented)
          } else{
            fit_segmented = segmented(fit_lm, seg.Z = ~lmod2, npsi = length(cpt@cpts))
            f = try(slope(fit_segmented),silent = T)
            if("try-error" %in% class(f)){
              slopes = coef(fit_segmented)
            }else{
              slopes = slope(fit_segmented)
              slopes = slopes[["lmod2"]][1:(length(cpt@cpts)+1)]}
            
          }
        } else{
          fit_segmented = segmented(fit_lm,seg.Z=~lmod,npsi = length(cpt@cpts))
          f = try(slope(fit_segmented),silent = T)
          if("try-error" %in% class(f)){
            slopes = coef(fit_segmented)
          }else{
            slopes = slope(fit_segmented)
            slopes = slopes[["lmod"]][1:(length(cpt@cpts)+1)]}
        }
      } else {
        fit_lm = lm(r~l)
        fit_segmented = segmented(fit_lm, seg.Z = ~l, npsi = length(cpt@cpts))
        f = try(slope(fit_segmented),silent = T)
        if("try-error" %in% class(f)){
          slopes = coef(fit_segmented)
        }else{
          slopes = slope(fit_segmented)
          slopes = slopes[["l"]][1:(length(cpt@cpts)+1)]}
      }
      
      
      #Add location of changepoints to new column "p"
      data[[i]]$p  = c(fit_segmented$psi[,2],rep(NA,nrow(data[[i]])-length(fit_segmented$psi[,2]))) 
      data[[i]]$c  = c(slopes,rep(NA,nrow(data[[i]])-length(slopes)))
      data[[i]]$cpx = c(rep(NA,nrow(data[[i]])))
      data[[i]]$cpy = c(rep(NA,nrow(data[[i]])))
      data[[i]]$cpz = c(rep(NA,nrow(data[[i]])))
      
      for (j in 1:length(fit_segmented$psi[,2])){
        p   = fit_segmented$psi[,2]
        if (length(p)!= 0){
          cpl = closest(l,p[j])
          data[[i]]$cpx[j] = data[[i]]$x[cpl]
          data[[i]]$cpy[j] = data[[i]]$y[cpl]
          data[[i]]$cpz[j] = data[[i]]$z[cpl]
        }
        else{
          data[[i]]$cpx[1] = "NO CHPT"
          data[[i]]$cpy[1] = "NO CHPT"
          data[[i]]$cpz[1] = "NO CHPT"
        }
      }}
    else{
      data[[i]]$p   = c(rep(NA,nrow(data[[i]])))
      data[[i]]$c   = c(rep(NA,nrow(data[[i]])))
      data[[i]]$cpx = c(rep(NA,nrow(data[[i]])))
      data[[i]]$cpy = c(rep(NA,nrow(data[[i]])))
      data[[i]]$cpz = c(rep(NA,nrow(data[[i]])))
      data[[i]]$cpx[1] = "NO CHPT"
      data[[i]]$cpy[1] = "NO CHPT"
      data[[i]]$cpz[1] = "NO CHPT"
    }
    
    
    #plot changepoints
    if (plot == TRUE){
      if (length(p)!= 0){
        plot(fit_segmented,col = 'blue',lwd = 3, main = names(data[i]), ylab = "Radius (cm)", xlab = "Length (cm)")
        points(l,r,col='gray',lwd = 1)
        lines.segmented(fit_segmented, col = 'red')
        points.segmented(fit_segmented, col = 'red')
        legend(x = 'topright',legend = c('Fit','Data','Changepoint'), col = c('blue','gray','red'), lty = c(1,NA,NA), pch = c(NA,1,1),lwd = 2)
      }
      else {
        plot(l,r)
      }} else{
        
      }
    
    
  }})
  write.xlsx(data,name,sheetname = names(data))
}