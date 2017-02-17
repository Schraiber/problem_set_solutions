while True:
	print "Input a number. Input -1 if you want to stop."
	theNumber = int(raw_input())
	if theNumber == -1:
		break
	if theNumber % 2 == 0:
		print "The number is even. Half the number is", theNumber/2
	else:
		print "The number is odd. Twice the number is", theNumber*2
