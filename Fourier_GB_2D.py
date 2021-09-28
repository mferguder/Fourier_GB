
def Fourier_GB(head,tail,surf,Lt,inorm):
    for Wavenumbers in range(1):
        Q=np.zeros((g,g))
        Qx=np.zeros((g,g))
        Qy=np.zeros((g,g))
        waves=np.array([[0,1],[1,0]])
        ignore=np.ones((g,g))
        ignore[0,1]=0
        ignore[1,0]=0
        for n in range(g):
            for m in range(g):
                if n+m==0: continue
                qx=(n%g)-g*(n>g/2-1)
                qy=(m%g)-g*(m>g/2-1)
                isqx=int(qx==-g/2)*2-1
                isqy=int(qy==-g/2)*2-1
                Qx[n,m]=qx
                Qy[n,m]=qy
                Q[n,m]=qx**2+qy**2
                if [qx,qy] in waves.tolist() or [isqx*qx,isqy*qy] in waves.tolist(): continue
                else:
                    waves=np.append(waves, np.array([[qx,qy]]), axis=0)
                    ignore[n,m]=0
        ignorer=np.ravel(ignore) # ignores conjugate amplitudes
        Qxr=np.ravel(Qx)
        Qyr=np.ravel(Qy)
        Q2r=np.ravel(Q)
        Qxr=Qxr[np.logical_not(ignorer)]
        Qyr=Qyr[np.logical_not(ignorer)]
        Q2r=Q2r[np.logical_not(ignorer)]
        Q2=np.unique(Q2r)
        q=Q2**.5*(2*np.pi/np.mean(Lt))

        Q2r_sorted=np.sort(Q2r)
        sorter=np.zeros((len(Q2r),len(Q2)))
        j=0
        for i in range(len(Q2r)):
            if Q2r_sorted[i]!=Q2r_sorted[i-1]: j+=1
            sorter[i,j-1]+=1
        sorter/=np.sum(sorter, axis=(0))[None,:]
    
    L=Lt[:,None]
    ln=np.sum((tail-head)**2,axis=2)**.5       # lipid length              time * N
    D_z=abs(tail[:,:,2]-head[:,:,2])           # lipid z                   time * N
    if inorm==0: norms=ln                      # normalization factors     time * N
    elif inorm==1: norms=D_z
    ns=(tail-head)/norms[:,:,None]             # lipid directors           time * N * 3
    n_agg=np.zeros((time,3,2,g,g))             # grid directors            time * 3 * 2 * g * g
    h_agg=np.zeros((time,2,g,g))               # grid surface              time * 2 * g * g

    z_c=np.mean(tail[:,:,2],axis=1)            # midplane                  time
    alphas=((head[:,:,2].T>z_c).T)             # leaflet index             time * N
    gxs=(tail[:,:,0]%L*g//L)                   # x grid index              time * N
    gys=(tail[:,:,1]%L*g//L)                   # y grid index              time * N
    gxs_h=(surf[:,:,0]%L*g//L)                 # x grid index              time * N
    gys_h=(surf[:,:,1]%L*g//L)                 # y grid index              time * N
                                               # tilted or flipping lipids time * N  
    ignore=(D_z/ln<.038)+(((head[:,:,2].T<z_c).T)!=(head[:,:,2]<tail[:,:,2])) 
    ignore_h=(D_z*0==1)
    
  # aggs: 1D list of 5D indexes (it, iN, alpha, gx, gy)
    for it in range(time):
        if it%400==0: print(it)
        aggs=(np.arange(N)*2*g*g+alphas[it]*g*g+gxs[it]*g+gys[it]).astype(int)
        Nagxy=np.zeros((N*2*g*g))             # grid tensor               N * 2 * g * g
        Nagxy[aggs]=np.ravel(np.logical_not(ignore[it]))
        Nagxy=Nagxy.reshape(N,2,g,g)
        k_n= np.sum(Nagxy, axis=(0))                 # degeneracies           2 * g * g
        k1 = k_n[None,:,:,:]
        k3 = np.repeat(k1, 3, axis=0)
        Nagxy/=k1

        aggs_h=(np.arange(N)*2*g*g+alphas[it]*g*g+gxs_h[it]*g+gys_h[it]).astype(int)
        Nagxy_h=np.zeros((N*2*g*g))
        Nagxy_h[aggs_h]=np.ravel(np.logical_not(ignore_h[it])) # h does not ignore tilted lipids
        Nagxy_h=Nagxy_h.reshape(N,2,g,g)
        k_h= np.sum(Nagxy_h, axis=(0))                 # degeneracies           time * 2 * g * g
        k1_h = k_h[None,:,:,:]
        k3_h = np.repeat(k1_h, 3, axis=0)
        Nagxy_h/=k1_h
    
#     ret = np.einsum('ijk,ijlmn->iklmn', ns, Nagxy)   actually slows down
        n_agg[it]=np.tensordot(ns[it], Nagxy, axes=(0,0))
        h_agg[it]=np.tensordot(surf[it,:,2], Nagxy_h, axes=(0,0))

        n_agg_neigh=np.array([np.roll(n_agg[it], 1, axis=-1),np.roll(n_agg[it], -1, axis=-1),
                              np.roll(n_agg[it], 1, axis=-2),np.roll(n_agg[it], -1, axis=-2)])
        h_agg_neigh=np.array([np.roll(h_agg[it], 1, axis=-1),np.roll(h_agg[it], -1, axis=-1),
                              np.roll(h_agg[it], 1, axis=-2),np.roll(h_agg[it], -1, axis=-2)])
        n_agg[it][k3==0]= np.nanmean(n_agg_neigh,axis=0)[k3==0]
        h_agg[it][k_h==0]= np.nanmean(h_agg_neigh,axis=0)[k_h==0]
    
    #normalization
    if inorm==1: n_agg[:,:2]/=abs(n_agg[:,None,2])
    elif inorm==0: n_agg[:,:2]/=(np.sum(n_agg**2,axis=1)**.5)[:,None,:]
    
    nqxyzs=np.fft.fft2(n_agg)*Lt[:,None,None,None,None]/g/g/2
    hqs=np.fft.fft2(h_agg)*Lt[:,None,None,None]/g/g/2
    nqx=np.reshape(nqxyzs[:,0,0]-nqxyzs[:,0,1],(len(n_agg),g*g))  #  divide by 2 ?
    nqy=np.reshape(nqxyzs[:,1,0]-nqxyzs[:,1,1],(len(n_agg),g*g))
    hqt=np.reshape(hqs[:,0]+hqs[:,1],(len(h_agg),g*g))
    nqx=nqx[:,np.logical_not(ignorer)]                            #     time * q_not_unique
    nqy=nqy[:,np.logical_not(ignorer)]
    hqt=hqt[:,np.logical_not(ignorer)]
    nllq2t=np.abs((Qxr*nqx+Qyr*nqy)**2/Q2r)                       #     time * q_not_unique
    nLq2t =np.abs((Qyr*nqx-Qxr*nqy)**2/Q2r)
    hq2t  =np.abs(hqt)**2
    nllq2t=np.matmul(nllq2t[:,Q2r.argsort()],sorter)              #     time * q_unique
    nLq2t=np.matmul(nLq2t[:,Q2r.argsort()],sorter)
    hq2t=np.matmul(hq2t[:,Q2r.argsort()],sorter)
    return hq2t,nllq2t,nLq2t,q # h_agg,n_agg
