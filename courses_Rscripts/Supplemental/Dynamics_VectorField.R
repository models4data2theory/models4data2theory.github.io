### FUNCTION DEFINITION

#! Outer product with vector function ###
# expand.outer joins functionality of "expand.grid" and "outer", so that
# a vector function can be applied, values of which are then stored in collumns.
# Combinations of "x" and "y" are expanded into vectors of length n*m, where
# n = length(x), m = length(y) and this ordering is preserved also in 'values'.
# Matrix of values (list item 'values') and expanded variables (list items 'x','y') are returned.
# The expand.outer function is called by the plotVectorField function that follows
expand.outer <- function(x, y, vecfun) {
  xy.pairs <- expand.grid(x = x,
                          y = y,
                          KEEP.OUT.ATTRS = FALSE)
  x.exp <- xy.pairs$x
  y.exp <- xy.pairs$y
  list(
    values = matrix(vecfun(x.exp, y.exp), nrow = 2, byrow = TRUE),
    x = x.exp,
    y = y.exp
  )
}

# vector field plot function
# grid.points can be defined for both axes at once or separately

plotVectorField <-
  function(vecfun,
           xlim,
           ylim,
           arrow.length = 0.08,
           grid.points = c(15, 15)) {
    gp <- if (length(grid.points) > 1)
      grid.points
    else
      rep(grid.points, 2)
    maxlength <- c(diff(xlim), diff(ylim)) / (gp - 1) * 0.9
    
    #prepare data
    x0 <- seq(xlim[1], xlim[2], length = gp[1])
    y0 <- seq(ylim[1], ylim[2], length = gp[2])
    xy.data <- expand.outer(x0, y0, vecfun)
    x0 <- xy.data$x
    y0 <- xy.data$y
    
    dx <- xy.data$values[1, ]
    dy <- xy.data$values[2, ]
    #scale
    k <- min(maxlength / c(max(abs(dx)), max(abs(dy))))
    x1 <- x0 + k * dx
    y1 <- y0 + k * dy
    #plot
    par(pty = 's', xaxs = 'i', yaxs = 'i')
    plot.default(
      range(x0, x1),
      range(y0, y1),
      xlab = "",
      ylab = "",
      type = "n",
      frame.plot = T
    )
    suppressWarnings(
    arrows(x0,
           y0,
           x1,
           y1,
           length = arrow.length,
           angle = 20,
           code = 2)
    )
  }

##########################################################################
##########################################################################
##########################################################################
### EXAMPLE

# For the 2-variable model
# dxdt=x*(1-x-0.6*y)
# dydt=y*(1-y-0.5*x)
# plotted for an x-range and y-range of 0 to 2

# plotVectorField(function(x,y) {c(x*(1-x-0.6*y),y*(1-y-0.5*x))}, xlim=c(0,2),ylim=c(0,2),arrow.length=0.03)
#
# # add isoclines
# iso.x<-function(x){-0.5*x + 1}
# iso.y<-function(x){-(x - 1)/0.6}
# curve(iso.x,add=T)
# curve(iso.y,add=T,lty=2)
#
# legend('topright',legend=c('x','y'),lty=c(1,2),bg='white',cex=0.5)
# title(xlab='x',ylab='y')
