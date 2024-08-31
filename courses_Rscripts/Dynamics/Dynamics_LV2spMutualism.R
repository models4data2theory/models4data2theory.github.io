# Script for running the 2 species mutualism model and
# generating a plot that can be manipulated.

# Version History:
# September 2023 (Sophia Simon)
# * original version for ECL 200A

### STANDARD VERSION OF MODEL ###

# Simulate Lotka-Volterra mutualism without stochasticity
lvMutualism = function(
	plotStyle="phase",
	N1_0=100,
	r1=0.1,
	K1=150,
	alpha12=0.5,
	N2_0=99,
	r2=0.1,
	K2=200,
	alpha21=0.5,
	endTime=100,
	printCode=FALSE
)
{
	require(deSolve, quietly=TRUE)
	
	if(printCode)
	{
		cat(sprintf(
			'lvMutualism(plotStyle=\"%s\", N1_0=%d, r1=%.1f, K1=%d, alpha12=%.1f, N2_0=%d, r2=%.1f, K2=%d, alpha21=%.1f, endTime=%d)\n',
			plotStyle, N1_0, r1, K1, alpha12, N2_0, r2, K2, alpha21, endTime
		))
	}
	
	dNdtFunc = function(t, states, param) 
	{
		with(as.list(c(states, param)), 
		{
			dN1dt = r1 * N1 * ((K1-N1+alpha12*N2) / K1)
			dN2dt = r2 * N2 * ((K2-N2+alpha21*N1) / K2)
			list(c(dN1dt, dN2dt))
		})
	}
	
	# Parameter, initial state, and times for deSolve
	params = c(r1=r1, K1=K1, alpha12=alpha12, r2=r2, K2=K2, alpha21=alpha21) 
	states = c(N1=N1_0, N2=N2_0) 
	times = seq(0, endTime, by = 0.01)
	
	# Redirect standard error output to /dev/null
	# sink("/dev/null")

	# Execute ODE
	out = as.data.frame(ode(y = states, times = times, func = dNdtFunc, parms = params, atol=1e-3, rtol=1e-3))
	
	# Stop redirecting standard error output
	# sink()
	
	mainTitle = sprintf(
		'N1_0=%d, r1=%.1f, K1=%d, alpha12=%.1f,\nN2_0=%d, r2=%.1f, K2=%d, alpha21=%.1f',
		N1_0, r1, K1, alpha12, N2_0, r2, K2, alpha21
	)
	
	# Plot results
	if(plotStyle == "time")
	{
		lvmutTimePlot(
			mainTitle=mainTitle,
			times=times, N1=out$N1, N2=out$N2
		)
	}
	else if(plotStyle == "phase")
	{
		lvmutPhasePlot(
			mainTitle=mainTitle,
			N1=out$N1, N2=out$N2,
			N1_0=N1_0, r1=r1, K1=K1, alpha12=alpha12,
			N2_0=N2_0, r2=r2, K2=K2, alpha21=alpha21
		)
	}
	
	return(invisible(out))
}
		
# Interactive controller
lvMutualismInteractive = function(
	N1_0=NULL,
	r1=NULL,
	K1=NULL,
	alpha12=NULL,
	N2_0=NULL,
	r2=NULL,
	K2=NULL,
	alpha21=NULL,
	endTime=NULL
)
{
	require(manipulate, quietly=TRUE)
	
	manipExpr = expression(invisible(
		lvMutualism(
			plotStyle=plotStyle,
			N1_0=N1_0,
			r1=r1,
			K1=K1,
			alpha12=alpha12,
			N2_0=N2_0,
			r2=r2,
			K2=K2,
			alpha21=alpha21,
			endTime=endTime,
			printCode=TRUE
		)))
	
	controls = list()
	controls$`_expr` = manipExpr
	controls$plotStyle = picker(
		"Phase Space" = "phase",
		"Abundances Over Time" = "time",
		label = "Plot Style"
	)
	if(is.null(N1_0))
		controls$N1_0 = slider(min = 1, max = 200, initial = 100, step = 1, ticks = FALSE)
	if(is.null(r1))
		controls$r1 = slider(min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE)
	if(is.null(K1))
		controls$K1 = slider(min = -10, max = 100, initial = 50, step = 1, ticks = FALSE)
	if(is.null(alpha12))
		controls$alpha12 = slider(
			min = 0.0, max = 1.5, initial = 0.2, step = 0.1, ticks = FALSE
		)
	if(is.null(N2_0))
		controls$N2_0 = slider(min = 1, max = 200, initial = 100, step = 1, ticks = FALSE)
	if(is.null(r2))
		controls$r2 = slider(min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE)
	if(is.null(K2))
		controls$K2 = slider(min = -10, max = 100, initial = 50, step = 1, ticks = FALSE)
	if(is.null(alpha21))
		controls$alpha21 = slider(
			min = 0.0, max = 1.5, initial = 0.5, step = 0.1, ticks = FALSE
		)
	if(is.null(endTime))
		controls$endTime = slider(
			min = 0, max = 400, initial = 100, step = 1, ticks = FALSE,
			label = "End Time"
		)
	
	do.call(manipulate, controls)
}


### PLOTTING ###

# Plot time course
lvmutTimePlot = function(
	mainTitle=NULL,
	times, N1, N2, K1, K2
)
{
  if (is.na(max(N1)) | is.na(max(N2))) {
    plot(times, N1[1:length(times)], type = "l", col="red",
         lwd=2, xlab="", ylab="",
         ylim=c(0, 100)
    )
    lines(times, N2[1:length(times)], col="blue", lwd=2,)
  } else {
    plot(times, N1, type = "l", col="red",
         lwd=2, xlab="", ylab="",
         ylim=c(0, 1.2*max(max(N1), max(N2)))
    )
    lines(times, N2, col="blue", lwd=2,)
  }
	
	
	legend('topright', c('N1', 'N2'), col=c("red", "blue"), lty=c(1,1))
	title(main=mainTitle, xlab="time", ylab="abundance", cex.main=0.8)
}

# Plot LV vector field

lvmutPhasePlot = function(
    mainTitle=NULL,
    N1, N2,
    N1_0, r1=0.1, K1=150, alpha12=0.5, N2_0, r2=0.1, K2=200, alpha21=0.5
)
{
  options(warn=-1) # Prevents zero-length arrow warnings from being printed
  
  plot.new()
  plot.window(
    xlim = c(1.3 * min(0, -K1/alpha12), 1.2 * max(N2_0, K2, K1/alpha12, (alpha21*K1+K2)/(1-alpha12*alpha21))),
    ylim = c(1.3 * min(0, -K2/alpha21), 1.2 * max(N1_0, K1, K2/alpha21, (alpha21*K2+K1)/(1-alpha12*alpha21)))
  )
  
  title(main=mainTitle, xlab=expression(N[2]), ylab=expression(N[1]), cex.main=0.8)
  
  axis(side=1, pos=0)
  axis(side=2, pos=0)
  
  plotaLine = function(x1, y1, x2, y2, col) {
    # Calculate the slope (m) and intercept (b) of the line
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    # Define a range for x values (e.g., from -100 to 100 for illustration)
    x_range = seq(-1000, 1000, by = 1)
    # Calculate the corresponding y values for the line
    y_range = m * x_range + b
    # Plot the line
    lines(x_range, y_range, col = col)
  }
  
  # Plot N2 nullcline: dN2/dt = 0
  # intercept on N2 axis: N1 = 0, N2 = K2
  # intercept on N1 axis: N1 = -K2/alpha21, N2 = 0
  #lines(x = c(0, K2), y = c(-K2/alpha21,0), col='red')
  plotaLine(K2,0,0,-K2/alpha21,'red')
  text(x = K2 + 0.15, y= -0.15 * K2, labels = expression(K[2]), col = 'red')
  text(x =0.15 * K2, y= -K2 / alpha21 -0.075*K2 / alpha21 , labels = expression(-K[2] / alpha[21]), col = 'red')
  
  # Plot N1 nullcline: dN1/dt = 0
  # intercept on N1 axis: N1 = K1, N2 = 0
  # intercept on N2 axis: N1 = 0, N2 = -K1/alpha12
  #lines(x = c(0, -K1/alpha12), y = c(K1,0), col='blue')
  plotaLine(0,K1,-K1/alpha12,0,'blue')
  text(x = -K1 / alpha12 - 0.12*K1 / alpha12, y=0.15 * K1, labels = expression(-K[1] / alpha[12]), col = 'blue')
  text(x = -0.075*K1, y=K1 + 0.075 * K1, labels = expression(K[1]), col = 'blue')
  
 
  dN1dtFunc = function(N)
  {
    r1 * N[1] * (1 - (N[1] - alpha12 * N[2])/K1)
  }
  
  dN2dtFunc = function(N)
  {
    r2 * N[2] * (1 - (N[2] - alpha21 * N[1])/K2)
  }
  
  # Calculate the maximum values of N1 and N2 based on parameters
  max_N1 = max(N1_0, K1, K2/alpha21, (alpha21*K2+K1)/(1-alpha12*alpha21))
  max_N2 = max(N2_0, K2, K1/alpha21, (alpha21*K1+K2)/(1-alpha12*alpha21))
  
  points = expand.grid(
    seq(0, 1.1 * max(N1_0, K1, K2/alpha21, (alpha21*K2+K1)/(1-alpha12*alpha21)), by=max_N1/10),
    seq(0, 1.1 * max(N2_0, K2, K1/alpha21, (alpha21*K1+K2)/(1-alpha12*alpha21)), by=max_N2/10)
    )
  
  dN1dt = apply(X=points, MARGIN=1, FUN=dN1dtFunc)
  dN2dt = apply(X=points, MARGIN=1, FUN=dN2dtFunc)
  
#  magnitude = sqrt(mean(dN1dt*dN1dt + dN2dt*dN2dt)) / 20
  
#  startN1 = points[,1] - dN1dt / magnitude
#  endN1 = points[,1] + dN1dt / magnitude
#  startN2 = points[,2] - dN2dt / magnitude
# endN2 = points[,2] + dN2dt / magnitude
  
#  arrows(y0=startN1, x0=startN2, y1=endN1, x1=endN2, length=0.05)

  points(y=N1[1], x=N2[1], pch=4,col='magenta')
  lines(y=N1, x=N2, type='l',lwd=2,col='magenta')
  points(y=N1[length(N1)], x=N2[length(N2)], pch=24,col='magenta')
  
  legend('topleft', c('N1 nullcline', 'N2 nullcline'),
         col=c("blue", "red"), lty=c(1,1))
}

# Plot time course for competition
lvTimePlot = function(
    mainTitle=NULL,
    times, N1, N2
)
{
  plot(times, N1, type = "l", col="red",
       lwd=2, xlab="", ylab="",
       ylim=c(0, 1.2*max(max(N1), max(N2)))
  )
  lines(times, N2, col="blue", lwd=2,)
  
  legend('topright', c('N1', 'N2'), col=c("red", "blue"), lty=c(1,1))
  title(main=mainTitle, xlab="time", ylab="abundance", cex.main=0.8)
}

# Plot LV vector field for competition
lvPhasePlot = function(
    mainTitle=NULL,
    N1, N2,
    N1_0, r1=0.1, K1=150, alpha12=0.5, N2_0, r2=0.1, K2=200, alpha21=0.5
)
{
  options(warn=-1) # Prevents zero-length arrow warnings from being printed
  
  plot.new()
  plot.window(
    xlim=c(0, 1.1 * max(N1_0, K1, K2/alpha21)),
    ylim = c(0, 1.1 * max(N2_0, K2, K1/alpha12))
  )
  
  title(main=mainTitle, xlab=expression(N[1]), ylab=expression(N[2]), cex.main=0.8)
  
  axis(side=1, pos=0)
  axis(side=2, pos=0)
  
  # Plot N1 nullcline: dN1/dt = 0
  # intercept on N1 axis: N1 = K1, N2 = 0
  # intercept on N2 axis: N1 = 0, N2 = K1/alpha12
  lines(x = c(K1, 0), y = c(0, K1/alpha12), col='red')
  mtext(text=expression(K[1]), at = K1, side = 1, line=-2, col='red')
  mtext(text=expression(K[1] / alpha[12]), at = K1/alpha12, side = 2, line=-2, col='red')
  
  # Plot N2 nullcline: dN2/dt = 0
  # intercept on N2 axis: N2 = K2, N1 = 0
  # intercept on N1 axis: N2 = 0, N1 = K2/alpha21
  lines(x = c(0, K2/alpha21), y = c(K2, 0), col='blue')
  mtext(text=expression(K[2] / alpha[21]), at = K2/alpha21, side = 1, line=-2, col='blue')
  mtext(text=expression(K[2]), at = K2, side = 2, line=-2, col='blue')
  
  dN1dtFunc = function(N)
  {
    r1 * N[1] * (1 - (N[1] + alpha12 * N[2])/K1)
  }
  
  dN2dtFunc = function(N)
  {
    r2 * N[2] * (1 - (N[2] + alpha21 * N[1])/K2)
  }
  
  points = expand.grid(
    seq(0, 1.1 * max(N1_0, K1, K2/alpha21, -K1), by=2),
    seq(0, 1.1 * max(N2_0, K2, K1/alpha12, -K2), by=2)
  )
  
  dN1dt = apply(X=points, MARGIN=1, FUN=dN1dtFunc)
  dN2dt = apply(X=points, MARGIN=1, FUN=dN2dtFunc)
  
#  magnitude = sqrt(mean(dN1dt*dN1dt + dN2dt*dN2dt)) / 20
  
#  startN1 = points[,1] - dN1dt / magnitude
#  endN1 = points[,1] + dN1dt / magnitude
#  startN2 = points[,2] - dN2dt / magnitude
#  endN2 = points[,2] + dN2dt / magnitude
  
#  arrows(y0=startN1, x0=startN2, y1=endN1, x1=endN2, length=0.05)
  
  points(x=N1[1], y=N2[1], pch=4,col='magenta')
  lines(x=N1, y=N2, type='l',lwd=2,col='magenta')
  points(x=N1[length(N1)], y=N2[length(N2)], pch=24,col='magenta')
  
  legend('topright', c('N1 nullcline', 'N2 nullcline'),
         col=c("red", "blue"), lty=c(1,1))
}

