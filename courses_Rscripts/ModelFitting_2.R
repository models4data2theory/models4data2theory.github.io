# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Your challenge is to fit the binomial model to 
# our data on fertilization success rates
# (i.e. counts of total and fertilized ovules)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import the GoogleSheet directly
dat <- read.csv('https://docs.google.com/spreadsheets/d/e/2PACX-1vTyyqthYKtedUdgRCOE37ec-oA4TzY6Mq8glR9bWr8ORhGQjWZlkeIuM5AgdGa8-zHE9pJma8C3n4_n/pub?gid=0&single=true&output=csv')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If that doesn't work, 
# then paste the above link into a browser 
# to save the GoogleSheet csv file to your desktop.

# Import downloaded csv file using either prompt...
# dat <- read.csv(file.choose())

# ... or by specifying its location.
# (un-comment following 2 lines to use)
# dat.loc <- 'MyDesktopLocation'
# dat <- read.csv(paste0(dat.loc, 'Delphinium_fertilization_success - Sheet1.csv'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove rows with no data
dat <- subset(dat, !is.na(dat$Total.count) 
                 & !is.na(dat$Fertilized.count))

# Extract vector of Total.counts (draws)
N <- dat$Total.count

# Extract vector of Fertilized.counts (successes)
k <- dat$Fertilized.count

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To do:

# 1) Define the negative log-likelihood function for binomial



# 2) Provide initial guess for binomial 'probability'



# 3) Fit the model



# 4) Overlay the parameter estimate on 
#    histogram of *Proportion* Fertilized




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


