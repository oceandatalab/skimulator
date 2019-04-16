import numpy as np
import matplotlib.pyplot as plt

#/********************************************************
#* bplsine_init_uniform *
#* *
#* initialize an allocated bspline struct with *
#* uniform knots between a and b *
#* *
#* - bspline : a pointer to an allocated bspline *
#* - a, b : interval of the bspline *
#* */

class spline_circular:
    def __init__(self,k,gp2):
 
        self.k = k
        self.g = gp2 - 1 # number of interior knots 
        self.Nknots = self.g + 2 + 2 * (self.k+1) # total number of knots 
        self.knots = np.zeros([self.Nknots])
 
        self.Nvals = self.g + self.k + 1
 
        self.Ncoefs = self.Nvals
        self.coefs = np.zeros([self.Ncoefs])
  
        self.dp = np.zeros([self.k + 1])
        self.dm = np.zeros([self.k + 1])

    def init_uniform(self,a,b):

        self.a = a
        self.b = b

        k = self.k
        g = self.g
        Nknots = self.Nknots

        self.knots = (np.arange(Nknots)-k)*(b-a)/(g+1)+a

        #for i in range(k):
        #    self.knots[i] = a;
        #
        #for i in np.arange(k,k + g+1,1):
        #    self.knots[i] = a + (i - k) / (g + 1.) * (b - a)
        #
        #for i in np.arange(k + g + 1,Nknots,1): 
        #    self.knots[i] = b
    
    def find_interval_uniform(self,x):
        
        k = self.k
        g = self.g
        a = self.a
        b = self.b

        j = np.floor(((x - a) / (b - a) * (g + 1.)) + k).astype('int')

        if (x < a).any():
            j[np.where(x < a)] = -1
        if (x > b).any():
            j[np.where(x > b)] = -1
                
        if (j > g + k).any():
                j[np.where(j > g + k)] = g+k

        return(j)

    def bspline_value (self, x):
        k = self.k
        g = self.g
        a = self.a
        b = self.b
        dp = self.dp
        dm = self.dm
        n=len(x)

    # Find interval 
        i=self.find_interval_uniform (x)
        if (i < 0).any():
            return([0,0,0,0],-1,-2)

    # interval of vals affected 
        istart = i - k
        iend = i
            
        vals = np.zeros([k+1,n])
        dp = np.zeros([self.k + 1,n])
        dm = np.zeros([self.k + 1,n])

        vals[0,:] = 1.
        for s in np.arange(1,k+1,1):
            dp[s - 1,:] = self.knots[i + s] - x
            dm[s - 1,:] = x - self.knots[i + 1 - s]
            prev = np.zeros([n])
            old = vals[0 ,:] 
            for r in range(1,s+1,1):
                M = old[:] / (dp[r - 1,:] + dm[s + 1 - r - 1,:])
                old = 1*vals[r ,:] 
                vals[r - 1 ,:] = prev[:] + dp[r - 1,:] * M[:]
                vals[r ,:] = dm[s + 1 - r - 1,:] * M[:]
                prev = 1*vals[r ,:]

        return (vals,istart%(g+1),iend)
     
    def init_val(self,yy):
        self.coefs=yy
        
    def get_val(self,xx):

        vals,istart,iend=self.bspline_value(xx)
        sy=0*xx
        for i in range(4):
            sy+=vals[i,:]*self.coefs[istart[:]+i]
        return(sy)

