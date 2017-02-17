theList = []
while True:
	print "Input a number. Put -1 if you want to stop."
	theInput = int(raw_input())
	if theInput == -1:
		break
	theList.append(theInput)
sumOfList = sum(theList)
newList = [theList[0]*theList[1], theList[2]*theList[3]]
print "The sum of all your numbers is", sumOfList
print "The new list is", newList
