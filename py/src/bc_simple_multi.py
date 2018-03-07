import time
from common import read_hypergraph,write_bc,adjSimpleEdge
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import deque

def bfs_simple(n,s,e):
	#data structure
	dist = {}
	pred = {}
	sigma = {}

	Q = deque()
	S = deque()

	#initialization
	dist[s] = 0;
	sigma[s] = 1;
	pred[s] = set();
	Q.append(s)
	while len(Q) > 0:
		u = Q.popleft()
		S.append(u)
		if u in e:
			for v in e[u]:
				if v not in dist :
					dist[v] = dist[u] + 1
					Q.append(v)
				if dist[v] == dist[u] + 1:
					if v not in pred:
						pred[v] = set()
					pred[v].add(u)
					if v not in sigma:
						sigma[v] = 0
					sigma[v] = sigma[v] + sigma[u]
	return pred,sigma,S

def bc_simple(n,e):

	startTime = time.time()
	newE= adjSimpleEdge(n,e)
	pool = ProcessPoolExecutor(max_workers=3)
	BC = [0.0]*n

	futures = [pool.submit(bc_simple_run, i,n,newE) for i in range(n)]
	j = 0
	for x in as_completed(futures):
		if(j%500==0):
			tempTime = time.time()
			#print(j,tempTime-startTime)
		j +=1
		value = x.result()
		for key in value:
			BC[key] += value[key]

	endTime = time.time()
	print('Total time:',endTime-startTime)
	return BC

def bc_simple_run(node,n,e):
	BC = {}
	#print(node,'start')
	pred,sigma,S = bfs_simple(n,node,e)
	delta = {}

	while len(S)>0:
		u = S.pop()	
		if u not in delta:
			delta[u] = 0
		for v in pred[u]:
			if v not in delta:
				delta[v] = 0
			delta[v] = delta[v] + float(sigma[v])/sigma[u]*(1+delta[u])
		if u is not node:
			if delta[u] > 0.0:
				BC[u]=delta[u]
	#print(node,'end')
	#print(len(BC))
	return BC

def bc_simple_main():
	dataFile = "test.txt"
	#coda
	#dataFile = "coda.txt"
	n,e = read_hypergraph("data/"+dataFile)	
	BC = bc_simple(n,e)
	write_bc(BC,"bc/bcr_"+dataFile)
	

if __name__ == '__main__':
    bc_simple_main()