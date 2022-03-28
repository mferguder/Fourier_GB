
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
def Blocking(xq2t,qq):
    l=len(xq2t[0])
    divs=divisors(l)
    blomax=int(0.5+np.log2(l))
    dxq2=np.zeros(len(qq))
    bloxydys=np.zeros((len(qq),3,len(divs)-2))
    for i in range(len(qq)):
        xi2t=xq2t[i]
        blox=np.log2(divs)
        blox=blox[0:-2]
        bloy=np.zeros(len(divs)-2)
        dbloy=np.zeros(len(divs)-2)
        for j in range(len(divs)-2):
            sz=divs[j]
            temp=np.zeros(int(np.floor(l/float(sz))));
            for k in range(int(np.floor(l/float(sz)))):
                temp[k]=np.mean( xi2t[ k*sz:(k+1)*sz ] )
            bloy[j]=np.std(temp)/np.sqrt(len(temp)-1)
            dbloy[j]=bloy[j]/np.sqrt(2*len(temp)-2)
        popt, pcov = curve_fit(BlockFunction, blox, bloy, p0=[bloy[0],0.2], sigma=dbloy)
        dxq2[i]=np.sqrt((1+popt[1])/(1-popt[1]))*bloy[0]
        #dxq2[i]=np.std(xq2t[i])
        bloxydys[i]=np.array([blox,bloy,dbloy])
    return dxq2#,bloxydys
#----------------------------------------------------------------------
def Spectra(qq,hq2t,nllq2t,nLq2t):
    z1=np.mean(hq2t,axis=0)*np.power(qq,4)
    z2=np.mean(nllq2t,axis=0)*np.power(qq,2)
    z3=np.mean(nLq2t,axis=0)
    dz1=Blocking(hq2t.T,qq)*np.power(qq,4)
    dz2=Blocking(nllq2t.T,qq)*np.power(qq,2)
    dz3=Blocking(nLq2t.T,qq)
    return z1,z2,z3,dz1,dz2,dz3
