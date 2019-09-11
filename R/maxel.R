#(broom)

conf_interval <- function(ci_x=NULL, ci_y=NULL) {
  # Computes maximum effect limit for specified confidence interval of x and optionally y
  # ci_x [1y numeric]: vector of lower and upper confidence limtis of x
  # ci_y [1y numeric]: vector of lower and upper confidence limtis of y
  # 
  # ci_x = c(.5, 2); ci_y = c(6, 9)
  # maxel.conf_interval(ci_x, ci_y)
  
  if (! is.numeric(ci_x)) stop("ci_x=NULL: confidence interval for x must be defined")
  
  #Check to make sure x is ordered correctly
  if (ci_x[2] < ci_x[1]) ci_x = rev(ci_x)
      
  if (is.numeric(ci_y)) {
    #Check to make sure x is ordered correctly
    if (ci_y[2] < ci_y[1]) ci_y = rev(ci_y)
    
    # Maximum difference between opposite bounds of interval
    maxel = max(abs(ci_x[2] - ci_y[1]),abs(ci_x[1] - ci_y[2]))

  } else {
    # CI limit furthest from zero, reported as a magnitude
    maxel = max(abs(ci_x))

  }
  return(maxel)
}


normal_unpaired <- function(x=NULL,y=NULL, var.equal=FALSE, 
                                  conf.level=0.96, alternative = "two.sided") {
  # Computes maximum effect limit for specified confidence interval of x and optionally y
  # ci [numeric vector]: samples from x
  # ci_y [numeric vector]: samples from y
  # 
  # x = c(4, 5, 3, 2, 4, 6, 4, 3, 4); y=c(8, 9, 7, 8, 10, 6, 5, 8, 7, 9, 8)
  # maxel.normal_unpaired(x,y)
  
  # Verify that numeric vector present for x1
  if (! is.numeric(x)) stop("ci=NULL: samples for x must be defined")
  
  if (is.numeric(ci_y)) {
    #Check to make sure x is ordered correctly
    if (ci_y[2] < ci_y[1]) ci_y = rev(ci_y)
    
    # Run ttest to acquire confidence interval
    tt = t.test(x,y,paired=FALSE,var.equal=var.equal, conf.level=conf.level,
                alternative=alternative)
    
    # Maximum difference between opposite bounds of interval
    maxel = maxel.conf_interval(tt$conf.int)
    
  } else {
    
    # Run ttest to acquire confidence interval
    tt = t.test(x,var.equal=var.equal, conf.level=conf.level,
                alternative=alternative)
    
    # CI limit furthest from zero, reported as a magnitude
    maxel=  maxel.conf_interval(tt$conf.int) 
  }
  return(maxel)
}
