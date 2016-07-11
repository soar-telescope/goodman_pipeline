#IRAF script to delete keywords known to contain non-ascii characters
#simon torres
#2016-06-29
#SOAR Observatory


string impat
impat = "*.fits"

print ("Removing keyworks known for containing NON-ASCII characters\n")
hedit (images=impat, fields="PARAM0,PARAM61,PARAM62,PARAM63", verify-, delete+)

print ("Bad keywords removed\n")

