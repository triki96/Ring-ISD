reset()

# compute partial gaussian elimination on ring-matrix
def PGE_ring(H):    
    n = H.ncols()
    r = H.nrows()
        
    ok = 1
    i = 0
    while (i<r)&(ok):
        #swap rows so that you have a pivot
        j = r-1-i
        while((not (H[j,n-1-i]).is_unit())&(j>0)):
            j = j-1
        if ((not (H[j,n-1-i]).is_unit())):
            ok = 0
        else: #swap rows
            tmp = H[r-1-i,:]
            H[r-1-i,:] = H[j,:]
            H[j,:] = tmp
            #scale row so that you have the pivot
            scale_coeff = H[r-1-i,n-1-i]^-1
            H[r-1-i,:] = scale_coeff*H[r-1-i,:]
            for v in range(r):
                if v!= i:
                    scale_coeff = H[r-1-v,n-1-i]
                    H[r-1-v,:] = H[r-1-v,:] - scale_coeff*H[r-1-i,:]
        i+=1;
    return ok, H;

# compute the GV bound for random linear codes
def gv_bound(n,k,p):
    w = 1
    while ((binomial(n,w)*(p-1)^w)<p^(n-k)):
        w +=1    
    return w-1

def swapCols(i,j,H):
    temp = H[:, i]
    H[:, i] = H[:, j]
    H[:, j] = temp
    return H

#Set code parameters
q = 4
p = 2
k = 5
n = 20

d = gv_bound(n,k,p) #compute minimum distance
w = floor(d/2) #error vector weight

print("[q, n, k, d] = ",[q, n, k, d])
print("Doing SDP with w = ",w)

#Generate code
Zq = Integers(q)
A = random_matrix(Zq,n-k,k)
H = block_matrix(Zq,1,2,[identity_matrix(Zq,n-k),A])

#Sample random error vector
e = matrix(Zq,1,n)
supp_e = Combinations(n,w).random_element() #random support
for i in supp_e:
    val = Zq.random_element()
    while val == 0:
        val = Zq.random_element()
    e[0,i] = val

#compute syndrome
s = e*H.transpose()


################# PROJECT THE CODE #####################

Fp = GF(p)
sub_H = H.change_ring(Fp)
sub_s = s.change_ring(Fp)

#We test Prange
solution_found = 0
while solution_found == 0:
    P = Permutations(n).random_element().to_matrix().change_ring(Fp) #sample permutation
    new_sub_H = sub_H*P #apply permutation
    if rank(new_sub_H[:,0:n-k])<(n-k):
        continue
    #Do Gaussian elimination
    U = new_sub_H[:,0:n-k]^-1
    new_sub_s = sub_s*U.transpose()
    #Count Hamming weight
    weight_candidate = n-k-new_sub_s.list().count(0)
    if weight_candidate <= w:
        #invert permutation
        perm_solution = matrix(Fp,1,n)
        perm_solution[0,0:n-k] = new_sub_s
        e_prime = perm_solution*P^-1
        solution_found = 1


################ CHECK PARTIAL SOLUTION ######################
        
print("Is it working? ",e.change_ring(Fp) == e_prime)
supp_e = vector(e).support()
    
print("e:  ", e)
print("e': ",e_prime)

for i in supp_e:
    print("Pos ",i,": ",e[0,i], "--> = ",e_prime[0,i])

################ FIND THE  ORIGINAL SOLUTION ###################
################ using these new informations ###################
    
# We use the recovered information about e_prime
e_prime_supp = [item[1] for item in e_prime.support()]


#We test Prange
solution_found = 0
while solution_found == 0:
	P = Permutations(n).random_element().to_matrix().change_ring(Zq) #sample permutation
	new_H = H*P #apply permutation
	for i in range(len(e_prime_supp)):
		new_H = swapCols(i,e_prime_supp[i],new_H) # metto a sx le colonne per cui il supp della sol è non nullo
	if PGE_ring(new_H[:,0:n-k])[0]==0:
		continue
	#Do Gaussian elimination
	U = new_H[:,0:n-k]^-1
	new_s = s*U.transpose()
	#Count Hamming weight
	weight_candidate = n-k-new_s.list().count(0)
	if weight_candidate <= w:
		# reverse columns
		perm_solution = matrix(Zq,1,n)
		perm_solution[0,0:n-k] = new_s
		for i in range(len(e_prime_supp)):
			perm_solution = swapCols(i,e_prime_supp[i],perm_solution) # rigiro le entrate
		#invert permutation
		e_prime = perm_solution*P^-1
		solution_found = 1

print("Find sol: ", e_prime)
