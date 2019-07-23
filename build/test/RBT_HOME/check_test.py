import sys

TOLERANCE = 0.01

try:
	ref=[item2.split() for item2 in [item.rstrip() for item in open(sys.argv[1]).readlines()]]
	test=[item2.split() for item2 in [item.rstrip() for item in open(sys.argv[2]).readlines()]]
except IndexError:
	print("Need to specify two file names.")
	sys.exit(2)
except FileNotFoundError:
	print("Specified file(s) not found.")
	sys.exit(1)

refcoords=ref[4:48]
testcoords=test[4:48]

refscore=float(ref[123][0])
testscore=float(test[123][0])

checkcoords=[]
for i in range(0,len(refcoords)):
	for j in [0,1,2]:
		checkcoords.append(abs(float(refcoords[i][j])-float(testcoords[i][j])))
error = False
for check in checkcoords:
	if check > TOLERANCE:
		error = True
if abs(refscore - testscore) > TOLERANCE:
	error = True

if error:
	print("The test failed, please check the compilation is OK and no errors were raised.")
else:
	print("The test succeeded! The results agree with the reference ones.\nHave fun using rDock!!")
