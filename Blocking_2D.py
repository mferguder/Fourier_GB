def divisors(n):
    temp=sorted(list(reduce(list.__add__, ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0))))
    divs=[] #list
    for i in range(len(temp)-1):
        if temp[i]!=temp[i+1]: divs.append(temp[i])
    divs.append( temp[len(temp)-1] )
    return divs
#----------------------------------------------------------------------
def BlockFunction(x,a,b):
    result=np.zeros(len(x))
    for i in range(len(x)): result[i]=np.real(a*cmath.sqrt(1-b*b-2*b*(1-np.power(b+0j,np.power(2,x[i])))/np.power(2,x[i]))/(1-b))
    return result
#----------------------------------------------------------------------
def Blocking_2D(xqt):
    data=xqt.transpose()
    g=data.shape[0]
    l=len(data[0,0])
    divs=divisors(l)
    blomax=int(0.5+np.log2(l))   
    dxq2=np.zeros((g,g))
    bloxydys=np.zeros((g,g,3,len(divs)-2))
    for ix in range(g):
        for iy in range(g):
            hi2t=data[ix,iy]
            blox=np.log2(divs)
            blox=blox[0:-2]
            bloy=np.zeros(len(divs)-2)
            dbloy=np.zeros(len(divs)-2)        
            for j in range(len(divs)-2):
                sz=divs[j]
                temp=np.zeros(int(np.floor(l/float(sz))));
                for k in range(int(np.floor(l/float(sz)))):
                    temp[k]=np.mean( hi2t[ k*sz:(k+1)*sz ] )
                bloy[j]=np.std(temp)/np.sqrt(len(temp)-1)
                dbloy[j]=bloy[j]/np.sqrt(2*len(temp)-2)
            popt, pcov = curve_fit(BlockFunction, blox, bloy, p0=[bloy[0],0.2], sigma=dbloy)
            dxq2[ix,iy]=np.sqrt((1+popt[1])/(1-popt[1]))*bloy[0]
            #dhq2[i]=np.std(hq2t[i])
            bloxydys[ix,iy]=np.array([blox,bloy,dbloy])
    return dxq2,bloxydys
#----------------------------------------------------------------------
def Get_ZdZ_n_rand(f,df,pows):
    #df/=np.sqrt(time-t0+1)
    Z = np.zeros(X.shape)
    dZ = np.zeros(X.shape)
    for n in range(g):
        for m in range(g):
            qx=(n%g)-g*(n>g/2-1)
            qy=(m%g)-g*(m>g/2-1)
            if min(qx,qy)<-gg or max(qx,qy)>gg: continue
            q=np.sqrt(qx**2+qy**2)*(2*np.pi/Lavg)
            Z[qx+gg,qy+gg]=np.random.normal( (f[n,m])*(q**pows), (df[n,m])*(q**pows),1 )
            dZ[qx+gg,qy+gg]=(df[n,m])*(q**pows)
    Z[gg,gg]=np.nan
    return Z,dZ

