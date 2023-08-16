#Sage code to accompany the paper "Kapranov degrees" by Joshua Brakensiek, Christopher Eur, Matt Larson, and Shiyue Li
#Author: Matt Larson
#
#
#
#


#fivelist = load('fivelist.sobj')
#sixlist = load('sixlist.sobj')
#above files must be loaded


#tests if a collection of indices satisfies the cerberus condition
#returns 0 if it does not satisfy, 1 if it satisfies
#number don't have to be 1 through n
def cerberus(u):
	dim = len(u)
	for s in subsets(range(dim)):
		if (len(s) == 0):
			continue
		total = u[s[0]]
		for i in s:
			total = total.union(u[i])
		if (len(total) - 3 < len(s)):
			return 0
	return 1


#takes a partition P, Q of n and a pair (set, index)
#P and Q should be sets
#returns 0 if it goes to P, 1 if it goes to Q, and the corresponding set on M_{0, P \cup 0}
def divy(P, Q, ind):
	if (P.issuperset({ind[1]})):
		restriction = ind[0].intersection(P)
		if (len(restriction) == 1):
			newset = ind[0].intersection(Q)
			return [1, [newset.union({0}), 0]]
		elif (len(restriction) == len(ind[0])):
			return [0, [restriction, ind[1]]]
		return [0, [restriction.union({0}), ind[1]]]
	restriction = ind[0].intersection(Q)
	if (len(restriction) == len(ind[0])):
		return [1, [restriction, ind[1]]]
	elif(len(restriction) > 1):
		return [1, [restriction.union({0}), ind[1]]]
	return [0, [ind[0].intersection(P).union({0}), 0]]

#takes a partition P, Q and a bunch of indices 
#outputs two sets of indices, one on M_{0, P \cup 0} and the other on M_{0, Q \cup 0}
def restrictToBoundary(P, Q, ind):
	Pindex = []
	Qindex = []
	for a in ind:
		result = divy(P, Q, a)
		if (result[0] == 0):
			Pindex.append(result[1])
		else:
			Qindex.append(result[1])
	return [Pindex, Qindex]



#takes a collection of n-3 subsets of [n] and i_j in each of the sets
#returns the Kapranov degree
#ex: kapranovDegree([[{1,2,3,4}, 1], [{1,2,3,4,5}, 2], [{3,4,5,6}, 6]])
#If check=False, then it will not check if the input is valid (consists of n - 3 subsets of [n], i_j lies in the jth set)
#if knownCerberus = True, then it will not check if the input satisfies Cerberus
def kapranovDegree(setsystem, check=True, knownCerberus=False):
	if (check == True):
		if (not isValidKapranov(setsystem)):
			print("Not valid")
			return False
	if (knownCerberus == False):
		a = [s[0] for s in setsystem]
		if (not cerberus(a)):
			return 0
	n = len(setsystem) + 3
	if (n == 3 or n==4):
		return 1
	if (n==5):
		for S in fivelist:
			if (S[0] == setsystem):
				return S[1]
	if (n==6):
		for S in sixlist:
			if (S[0] == setsystem):
				return S[1]

	#if all psi classes are the same and the system satisfies Cerberus, then degree is 1
	allpsi = set([s[1] for s in setsystem])
	if (len(allpsi) == 1):
		return 1
	#removes an element from the first set, running the recursion
	list1 = list(setsystem[0][0])
	removed = list1[0]
	if (list1[0] == setsystem[0][1]):
		removed = list1[1]
	newsetsystem = setsystem.copy()
	deletedsetsystem = setsystem.copy()
	deletedsetsystem.remove(setsystem[0])
	newsetsystem[0][0].remove(removed)
	total = set(range(1, n+1))
	rest = total.difference(setsystem[0][0])
	rest = rest.difference({removed})
	partitions = subsets(rest)
	contribrest = 0
	for part in partitions:
		part = set(part)
		P = part.union({setsystem[0][1], removed})
		Q = total.difference(P)
		restrictions = restrictToBoundary(P, Q, deletedsetsystem)
		if (len(restrictions[0]) != len(P) - 2):
			continue
		if (len(restrictions[1]) != len(Q) - 2):
			continue
		contribrest += kapranovDegree(standardForm(restrictions[0]), check=False)*kapranovDegree(standardForm(restrictions[1]), check=False)
	return contribrest + kapranovDegree(newsetsystem, check=False)


#returns 1 if the set system is a valid Kapranov degree
def isValidKapranov(setsystem):
	groundset = set()
	for S in setsystem:
		if (not set([S[1]]).issubset(S[0])):
			return 0
		groundset = groundset.union(S[0])
	if (len(setsystem) != len(groundset) - 3):
		return 0
	return 1

#takes a set system on P \cup 0
#returns a corresponding set system on {1, ... |P| + 1}
def standardForm(setsystem):
	groundset = set()
	reindexing = []
	for S in setsystem:
		newelements = S[0].difference(groundset)
		groundset = groundset.union(S[0])
		for j in list(newelements):
			reindexing.append((j, len(reindexing) + 1))
	reindexdict = dict(reindexing)
	newsets = []
	for S in setsystem:
		newset = set([reindexdict[j] for j in S[0]])
		newsets.append([newset, reindexdict[S[1]]])
	return newsets

#takes a collection of n-3 subsets of size 4, each contained in [n]
#computes the cross ratio degree
#ex: crossRatioDeg([{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16},{1,5,9,13},{1,6,10,14},{1,7,11,15},{1,8,12,16},{2,5,10,15},{3,5,12,14},{4,6,12,15},{2,6,11,16},{3,7,9,13}])
#if check = False, will not check if the input is valid
def crossRatioDeg(setsystem, check=True):
	#converts problem into a Kapranov degree
	if (not isValidCross(setsystem)):
		print("Not valid")
		return False
	kapsetsystem = []
	for S in setsystem:
		listS = list(S)
		kapsetsystem.append([S, listS[0]])
	return kapranovDegree(kapsetsystem)

#returns 1 if the set system is a valid cross ratio degree
def isValidCross(setsystem):
	groundset = set()
	for S in setsystem:
		if (len(S) != 4):
			return 0
		groundset = groundset.union(S)
	if (len(setsystem) != len(groundset) - 3):
		return 0
	return 1

