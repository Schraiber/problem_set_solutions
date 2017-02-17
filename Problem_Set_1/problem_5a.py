theList = []
for i in range(4):
	print "Input a number"
	theList.append(int(raw_input()))
sumOfList = sum(theList)
newList = [theList[0]*theList[1], theList[2]*theList[3]]
print "The sum of all your numbers is", sumOfList
print "The new list is", newList
