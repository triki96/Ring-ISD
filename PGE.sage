#reset();

#r = 10;
#n = 20;
#ell = 6;
#q = 3;

#Fq = GF(q);
#H = random_matrix(Fq,r,n);

def PGE(n,r,ell,Fq,H):

#copy_H = copy(H);

    ok = 1;
    i = 0;
    while (i<(r-ell))&(ok):

        #swap rows so that you have a pivot
        j = r-1-i;
        while(H[j,n-1-i]==0)&(j>=0):
            j = j-1;

        #if no valid element has been found, report failure
        if (H[j,n-1-i]==0):
            ok = 0;

        else: #swap rows
            tmp = H[r-1-i,:];
            H[r-1-i,:] = H[j,:];
            H[j,:] = tmp;

            #scale row so that you have the pivot
            scale_coeff = H[r-1-i,n-1-i]^-1;
            H[r-1-i,:] = scale_coeff*H[r-1-i,:];

            #create zeros
            for v in range(r):
                if v!= i:
                    scale_coeff = H[r-1-v,n-1-i];
                    H[r-1-v,:] = H[r-1-v,:] - scale_coeff*H[r-1-i,:];

#            print(H)
#            print("--------------");
        i+=1;

    return ok, H;
