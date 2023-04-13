G = AdditiveAbelianGroup([16,16,16,64,64,64])#The Group description here, only abelian p groups
B = G.gens()#The basis is e_1,e_2,...,e_N
Q = list(B)
p=2 #the prime p
e = [int(math.log(q.order(),p)) for q in Q]
QQ = [1,2,4,2,8,16]#the q coefficients
K = sum([QQ[i]*Q[i] for i in range(len(QQ))])
M = [2,1,4,8,16,2]#the m's
M2 = list(M)

I_M = set(range(len(Q)))
I_Q = set(range(len(Q)))

def shuffle():
    M1=list(M)
    Q1 = list(QQ)
    a=1
    p_ord=0
    while true:
        if all(c==0 for c in M1):
            break
        l=k=-1
        M1 = [c//a for c in M1]
        Q1 = [c//a for c in Q1]
        #print("M1:",M1)
        #print("Q1:",Q1)
        for i in I_M:
            if M1[i]%p!=0 and i>k:
                k=i
        #print("k:",k)
        for j in I_Q:
            if Q1[j]%p!=0 and j>l:
                l=j
        #print("l:",l)
        if Q[k].order() > p_ord and l!=k and l!=-1 and k!=-1:
            Q[k],Q[l]=Q[l],Q[k]
            QQ[k],QQ[l]=QQ[l],QQ[k]
        p_ord = Q[k].order() 
        a=p
        #print("p_ord: ",p_ord)
        #print(QQ)
    K = sum([QQ[i]*Q[i] for i in range(len(QQ))])
    return Q
def Reduce(K):
    while True:
        if all(M[i]%p==0 for i in I_M):
            for i in I_M:
                M[i]=M[i]//p
            for i in I_Q:
                QQ[i]= QQ[i]//p
            continue
        else:
            break
    
    K= sum([QQ[i]*Q[i] for i in I_Q])
    return K

def check(K,M):
    if K.order()== 1 or all(a == 0 for a in M) :
        if all(a == 0 for a in M) and K.order()==1 :
            print(Q)
            return 1
        else:
            print("No solution exists")
            return 1
    L = sum([M[i]*Q[i] for i in range(len(Q))])
    if K.order() == L.order() and all(rho(p**e,QQ)==rho(p**e,M) for e in range(int(math.log(K.order(),p)))):
        
        return 0
    else:
        print("No solution exists")
        return 1
    
def rho(l,S):
    r = 0
    S = [(l*S[i])%Q[i].order() for i in range(len(S))]
    while True:
        if all(a%p==0 for a in S):
            r= r+1
            S = [a/p for a in S]
        else:
            break
    return r

def solve(K):
    while true:
        K = Reduce(K)
        k = 0
        #print("I_M: ", I_M)
        for i in I_M:
            if M[i]%p != 0 and i>k:
                k = i
                
        #print(M)
        if K.order() == p**e[k]:
            P[k] = inv(M[k],k)*(K-sum([M[i]*Q[i] for i in I_M])+M[k]*Q[k])
            return P
        else:
            
            I_M1= set()
            I_Q1= set()
            #print("I_M: ", I_M, "I_M1: ", I_M1)
            #print("I_Q: ", I_Q, "I_Q1: ", I_Q1)
            for i in I_M:
                if (p**e[k]*M[i]*Q[i])== G.zero():
                    I_M1.add(i)
            #print("k: ", k,"I_M: ", I_M, "I_M1: ", I_M1)
            for i in I_Q:
                if (p**e[k]*QQ[i]*Q[i])== G.zero():
                    I_Q1.add(i)
            #print("I_Q: ", I_Q, "I_Q1: ", I_Q1)
            K1 = sum([QQ[i]*Q[i] for i in I_Q1])
            P[k]= inv(M[k],k)*(K1-sum([M[i]*Q[i] for i in I_M1])+M[k]*Q[k])
            I_M.difference_update(I_M1)
            #print("k: ", k,"I_M: ", I_M, "I_M1: ", I_M1)
            
            I_Q.difference_update(I_Q1)
            K = K-K1
            #print(K,"I_M:",I_M)
            
                        

def inv(m,j):
    q1 = m
    q2 = p**e[j]
    g,x,y = gcdExtended(q1,q2)
    return x%p**e[j]
    
def gcdExtended(a, b):
    # Base Case
    if a == 0 :
        return b,0,1
             
    gcd,x1,y1 = gcdExtended(b%a, a)
     
    # Update x and y using results of recursive
    # call
    x = y1 - (b//a) * x1
    y = x1
     
    return gcd,x,y
print("G: ",G)
print("K: ",K)
print("M: ",M)
#print("after shuffling K: ",K)
#print("after shuffling QQ: ",QQ)
#print("after shuffling Q: ",Q)

Result=check(K,M)
if Result != 1:
    P = list(shuffle())
    print("solution: ",solve(K))
    print(sum([M2[i]*P[i] for i in range(len(QQ))]))
else:
    quit()
