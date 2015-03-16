### utilities have been normalized to fall in the interval [0,1]
### to obtain the  expeted utility with unscaled utilities apply the transformation 20*MEU-2H, where MEU is the value obtained with this formulation
#
## HORIZON
import sys
if len(sys.argv) < 2:
    print "Usage:", sys.argv[0], "H, where H is a positive integer specifying the planning horizon (number of steps)."
    exit(0)
H = int(sys.argv[1])
A = 21 # position sA
B = 25 # position sB
states = [1, 2, 3, 4, 5, 9, 14, 13, 12, 17, 16, 21, 22, 19, 20, 24, 25] # possibility space of S
actions = ['r','n','e','s','w'] # possibility space of D
start = 1 # initial state
VARS =  ['S',         'N','E','SO','W','A','B','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15','X16']
sizeC = [ len(states),  2,  2,   2,  2,  2,  2,   5,   5,   5,   5,   5,   5,   5,   5,   5,    5,    5,    5,    5,    5,    5,    5]
print "/**********************"
print " * MAZE LIMID"
print " * HORIZON:", H
print " * START:", start
# probabilistic finite state machine specifying the state transition distribution
# (s,s',action) -> prob
#          1 2 3 4 5       6 7 8 9 10
grid = [[0,0,0,0,0,0,0],[0,1,1,1,1,1,0],[0,0,0,0,1,0,0],[0,0,1,1,1,0,0],[0,1,1,0,1,1,0],[0,1,1,0,1,1,0],[0,0,0,0,0,0,0]]
for y in range(7):
    print " *",
    for x in range(7):
        if grid[6-y][x]:
            print "%2d" % (5*(6-y-1)+x),
        else:
            print "  ",
    print
T = {}
for y in range(1,6):        
    for x in range(1,6):
        s = 5*(y-1)+x
        if grid[y][x]:
                #states.append(s)
                T[(s,s,'r')] = 1.0 # rest
                if grid[y+1][x]:
                    # move north
                    T[(s,5*y+x,'n')] = 0.9
                    T[(s,s,'n')] = 0.1
                else:
                    T[(s,s,'n')] = 1.0
                if grid[y][x+1]:
                    # move east
                    T[(s,s+1,'e')] = 0.9
                    T[(s,s,'e')] = 0.1
                else:
                    T[(s,s,'e')] = 1.0
                if grid[y-1][x]:
                    # move south
                    T[(s,5*(y-2)+x,'s')] = 0.9
                    T[(s,s,'s')] = 0.1
                else:
                    T[(s,s,'s')] = 1.0
                if grid[y][x-1]:
                    # move west
                    T[(s,5*(y-1)+x-1,'w')] = 0.9
                    T[(s,s,'w')] = 0.1
                else:
                    T[(s,s,'w')] = 1.0
                    
print " **********************/"
print "LIMID" # header
## NUM VARS (CHANCE, DECISION, VALUE)
print len(VARS)*H, 16*H, 2*H
vid = {} # dictionary mapping variables into their indexes
## CARDINALITIES
# CHANCE VARS
for t in range(H):
    for i in range(len(VARS)):
        vid[VARS[i]+'_'+str(t+1)]=len(vid)
        print sizeC[i],
    print
# DEC VARS
for t in range(H):
    for i in range(1,17):
        vid['D'+str(i)+'_'+str(t+1)]=len(vid)
        print "5",
print

## PARENT SETS
# CHANCE VARS
for t in range(H):
    h = '_'+str(t+1)
    hp = '_'+str(t)
    if t == 0:
        print 0   # S
    else:
        print 2, '\n\t', vid['S'+hp], vid['X16'+hp] # S' | S X16
    print 1, '\n\t', vid['S'+h] # N | S
    print 1, '\n\t', vid['S'+h] # E | S
    print 1, '\n\t', vid['S'+h] # SO | S
    print 1, '\n\t', vid['S'+h] # W | S
    if t == 0:
        print 0 # A
        print 0 # B
    else:
        print 2, '\n\t', vid['S'+hp], vid['A'+hp] # A' | S A
        print 2, '\n\t', vid['S'+hp], vid['B'+hp] # B' | S B
    print 5, '\n\t', vid['D1'+h], vid['N'+h], vid['E'+h], vid['SO'+h], vid['W'+h] # X1 | D1 N E SO W
    for i in range(2,17):
        print 6, '\n\t', vid['X'+str(i-1)+h], vid['D'+str(i)+h], vid['N'+h], vid['E'+h], vid['SO'+h], vid['W'+h] # Xi | Xi-1 Di N E SO W
# DEC VARS
for t in range(H):
    h = '_'+str(t+1)
    for i in range(16):
        print 0 # Di+1
# UTIL VARS
for t in range(H):
    h = '_'+str(t+1)
    print 3, '\n\t', vid['S'+h], vid['A'+h], vid['B'+h] # U(S,A,B)
    print 1, '\n\t', vid['X16'+h] # U(X16) 

    
## CPTS
for t in range(H):
    if t==0:
        # P(S)
        print len(states)
        print '\t',
        for s in states:
            if s == start:
                print "1.0",
            else:
                print "0.0",
        print
    else:
        # P(S'|S,X16)
        print len(states)*len(states)*5
        for d in actions:
            for s1 in states:
                print '\t',
                for s2 in states:
                    if (s1,s2,d) in T:
                        print T[(s1,s2,d)],
                    else:
                        print 0.0,
                print
    # P(N|S)
    print 2*len(states)
    for s in states:
        x,y = ((s-1) % 5)+1, int((s-1)/5)+1
        if grid[y+1][x]==0:
            print "\t0.99 0.01"
        else:
            print "\t0.1  0.9"
    # P(E|S) 
    print 2*len(states)
    for s in states:
        x,y = ((s-1) % 5)+1, int((s-1)/5)+1
        if grid[y][x+1]==0:
            print "\t0.99 0.01"
        else:
            print "\t0.1  0.9"
    # P(SO|S)
    print 2*len(states)
    for s in states:
        x,y = ((s-1) % 5)+1, int((s-1)/5)+1
        if grid[y-1][x]==0:
            print "\t0.99 0.01"
        else:
            print "\t0.1  0.9"
    # P(W|S)
    print 2*len(states)
    for s in states:
        x,y = ((s-1) % 5)+1, int((s-1)/5)+1
        if grid[y][x-1]==0:
            print "\t0.99 0.01"
        else:
            print "\t0.1  0.9"
    if t==0:
        # P(A)
        print "2\n\t0.0 1.0"
    else:
        # P(A'|S,A)
        print 2*len(states)*2
        for a in [ 'y', 'n' ]:
            for s in states:
                if a=='y' or s==A:
                    print "\t1.0 0.0"
                else:
                    print "\t0.0 1.0"
    if t==0:
        # P(B)
        print "2\n\t0.0 1.0"
    else:
        # P(B'|S,B)
        print 2*len(states)*2
        for b in [ 'y', 'n' ]:
            for s in states:
                if b=='y' or s==B:
                    print "\t1.0 0.0"
                else:
                    print "\t0.0 1.0"
    # P(X1|D,N,E,SO,W)
    print 5*5*16
    for w in range(2):
        for so in range(2):
            for e in range(2):
                for n in range(2):
                    for d1 in range(5):
                        print '\t',
                        for x1 in range(5):
                            if (n,e,so,w)==(0,0,0,0) and (x1==d1):
                                print 1.0,
                            elif (n,e,so,w)==(0,0,0,0) and (x1!=d1):
                                print 0.0,
                            else:
                                print 0.2,
                        print
    for p in range(1,16):
        print 5*5*5*16
        # P(Xi+2|Xi+1,D,N,E,SO,W)
        for w in range(2):
            for so in range(2):
                for e in range(2):
                    for n in range(2):
                        for d2 in range(5):
                            for x1 in range(5):
                                print '\t',
                                for x2 in range(5):
                                    if ((n+2*e+4*so+8*w!=p) and (x2==x1)) or ((n+2*e+4*so+8*w==p) and (x2==d2)):
                                        print 1.0,
                                    else:
                                        print 0.0,
                                print

## UTILS
for t in range(H):
    # U(S,A,B)
    print len(states)*2*2
    print '\t',    
    for b in ['y','n']:
        for a in ['y','n']:
            for s in states:
                if a=='n' and b=='n' and s==B: # visited B first
                    print 0.25,
                elif a=='n' and b=='n' and s==A: # visited A first
                    print 0.5,
                elif a=='y' and b=='n' and s==B: # visited A then B
                    print 1.0,
                else:
                    print 0.05,
    print
    # U(X16)    
    print "5\n\t0.05  0.0  0.0  0.0  0.0"
