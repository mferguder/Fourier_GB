import scipy.special as sc
maxlim=32
pows=4
Nboot=200
badness=0.9999999
qmax=5.5
def Bootstrp_2D(qx,qy,zx,dzx,xq_func,_xq_func,px,qmax):
#     lim=Get_lim(maxlim,zx,dzx,xq_func,_xq_func,1-1e-5)
#     Z,dZ=Get_ZdZ_n_rand(zx,dzx)
    xdatar=np.vstack((qx, qy))
#     xdatar,Zr,dZr=Get_fit_data(xdata,Z,dZ,lim)
    
    a,b = curve_fit(_xq_func,   xdatar, zx, px, sigma=dzx)
    chi_s=np.power(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx),2)
    c=sc.gammainc((len(qxr)-len(a))/2.,sum(chi_s)/4.)   
#     c=stats.chi2.cdf(sum(chi_s) , len(qx)-len(a)-1)
    
    qzdzx=np.concatenate((np.array([qx]),np.array([qy]),np.array([zx]),np.array([dzx])))
    params_xq=np.zeros((Nboot,len(a)))
    dparams_xq=np.zeros((Nboot,len(a),len(a)))
#     if lim==4: return params_xq,qzdzx

    while c>badness or qmax**2<qx[-1]**2+qy[-1]**2:
        if c<1: print( len(qx),c,np.mean(abs(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx))))
        qx=qx[:-1]
        qy=qy[:-1]
        zx=zx[:-1]
        dzx=dzx[:-1]
        xdatar=xdatar[:,:-1]
        a,b = curve_fit(_xq_func,   xdatar, zx, px, sigma=dzx)
        chi_s=np.power(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx),2)
        c=sc.gammainc((len(qx)-len(a))/2.,sum(chi_s)/4.)
#         c=stats.chi2.cdf(sum(chi_s) , len(qx)-len(a)-1)
    print( len(qx),c,np.mean(abs(np.divide(np.subtract(zx,xq_func(qx,qy,*a)),dzx))))
    qzdzx=np.concatenate((np.array([qx]),np.array([qy]),np.array([np.sqrt(qx**2+qy**2)]),np.array([zx]),np.array([dzx])))
    
    i=0
    while i<Nboot:
        zxk=np.random.normal(zx,dzx,len(qx))
        try:
            params_xq[i], dparams_xq[i] = curve_fit(_xq_func,   xdatar, zx, px, sigma=dzx)
            i+=1
        except RuntimeError: print("Error - curve_fit failed")
    params_xq =np.transpose(params_xq)
#     print( params_xq,qzdzx)
    return params_xq,qzdzx

guess_prms = [(1./23,1./25)]
p0 = [p for prms in guess_prms for p in prms]

params_hq,qzdz1=Bootstrp_2D(qxr,qyr,z1,dz1,new_hq,_new_hq,p0,qmax)
