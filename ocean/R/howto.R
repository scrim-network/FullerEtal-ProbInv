# written by Robert W. Fuller on 20101003

source('ipcc.R')

outfiles <- F

# show hindcast
ipccPlotOceanTemp(outfiles=outfiles)

# show predictions (not drift corrected)
ipccPlotPredict(outfiles=outfiles)

# show control runs
ipccPlotControl(outfiles=outfiles)

# show drift corrected predictions
ipccPlotDriftCorrect(outfiles=outfiles)

# generate ocean_temp.txt
source('nathan.R')

# write ocean_temp.txt
writeNathan()

# read ocean_temp.txt
data <- readNathan()

# print the data from ocean_temp.txt
#print(data)
