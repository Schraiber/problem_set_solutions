def product_of_first_and_last(theList):
	product = theList[0]*theList[-1] #alternatively, theList[0]*theList[len(theList)-1]
	return product

userList = []
for i in range(5):
	print "Input a number"
	theNumber = int(raw_input())
	userList.append(theNumber)

theProduct = product_of_first_and_last(userList)
print "The product of the first and last entries is", theProduct
