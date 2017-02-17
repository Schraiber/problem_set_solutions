def compute_min_max(l):
	#need to convert the entries in the list from strings to ints
	intList = []
	for item in l:
		intList.append(int(item))
	#use built-in functions to compute min and max
	theMin = min(intList)
	theMax = max(intList)
	#return two things
	return theMin, theMax

myFile = open("tab_separated.txt")
myOutput = open("problem_8_output.txt","w")
for line in myFile:
	strippedLine = line.strip()
	splitLine = strippedLine.split()
	myMin, myMax = compute_min_max(splitLine)
	myOutput.write("The min of the line is %d and the max is %d\n"%(myMin,myMax))
	
	
