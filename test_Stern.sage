reset();
load('list_sorting.sage');
#load('HJ_W.sage');
#load('attack_utils.sage')
load('stern_utils.sage');

#######################################################


##################################

#Use a seed for the simulation
seed_val = randrange(2^128);
set_random_seed(seed_val);

#Set code parameters and attack parameters

num_codes = 10000000;

n = 100; #code length
k = 50;
q = 2; #finite field size

Fq = GF(q);

Fq_star = [x for x in Fq.list()[1:]];

H = random_matrix(Fq,n-k,n)

#parameters
w = 13;
p = 2;
ell = 6;

t = 10; #num of ISD calls


#Sample H
#filename = "./Results/values_"+str(n)+"_"+str(w)+"_"+str(p)+"_"+str(ell);

#with open(filename, "w") as f:
#    f.write(" ");
#f.close();


avg_hat_N = 0;

for num_attempt in range(t):

    print(num_attempt);

    #Apply permutation and PGE
    X = stern_isd(H, Fq, Fq_star, n, k, p, ell, w);
    avg_hat_N+=len(X);


succ_pr = binomial(floor((k+ell)/2), p)*binomial(ceil((k+ell)/2), p)*binomial(n-k-ell,w-2*p)/binomial(n,w);

emp_Nw = avg_hat_N/(t*succ_pr);
th_Nw = binomial(n,w)*(q-1)^(w-1)*q^(k-n);
print("Emp. = ",emp_Nw*1.,", Th. = ", th_Nw*1.);


#        with open(filename, "a") as f:
#            f.write("\nNum Codes = "+str(id_code));
#            f.write("\nEmp Num Coll = "+str(average_num_coll/(id_code + 1.))+", Th Num Coll = "+str(th_num_collision));
#            f.write("\nEmp Num Sol = "+str(average_num_found/(id_code + 1.))+", Th Num Sol = "+str(th_num_sol));
#            f.write("\n--------------------------------");
#        f.close();
