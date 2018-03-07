from heapq import heappush,heappop

def rearrange(Q,e,givenNodes):
	newQ = []
	for element in Q:
		towardNext, fedge = element
		heappush(newQ,(len(e[fedge][0]-givenNodes),fedge))
	return newQ

def bfs_hyperpath(x,S,e,fstar):
	reachable = set()
	givenNodes = set()
	visited = set()
	Q = []
	#initialization
	reachable.add(x)
	givenNodes.add(x)
	for v in S:
		givenNodes.add(v)
	if x in fstar:
		for edge in fstar[x]:
			visited.add(edge)
			heappush(Q,(len(e[edge][0]-givenNodes),edge))

	while len(Q) > 0:
		towardNext, fedge = heappop(Q)
		#source condition
		if towardNext == 0:			
			t = e[fedge][1]
			if t not in reachable:
				reachable.add(t)
				givenNodes.add(t)
				Q = rearrange(Q,e,givenNodes)
				if t in fstar:
					for edge in fstar[t]:
						if edge not in visited:
							visited.add(fedge)
							heappush(Q,(len(e[edge][0]-givenNodes),edge))

	'''
	visitedN = set()
	visitedE = set()
	Q = deque()

	for v in X:
		visitedN.add(v)
		Q.append(v)
	for v in S:
		visitedN.add(v)

	while len(Q) > 0:
		v = Q.popleft()
		if v in fstar:
			for fIndex in fstar[v] - visitedE:
				t = e[fIndex][1]
				if t in Q:
					continue;
				towardNext = True
				for s in e[fIndex][0]:
					if s not in visitedN:
						towardNext = False
						break;
				if towardNext:
					visitedE.add(fIndex)
					Q.append(t)
					visitedN.add(t)
	return visitedN
	'''
	return reachable