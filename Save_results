def UseParams_hq(Nboot,params_hq):
    j=[] #list
    a=np.divide(params_hq[1],np.multiply(np.multiply(params_hq[0],params_hq[2]),4))
    Rm1s_hq=np.zeros(Nboot)
    Rm2s_hq=np.zeros(Nboot)
    for i in range(Nboot):
        if np.imag(cmath.sqrt(1+4*a[i]))==0: Rm1s_hq[i],Rm2s_hq[i]=np.roots([a[i],-1,-1])
        else: j.append(i)
    for ii in range(len(j)):
        i=j[len(j)-ii-1]
        params_hq=[list(params_hq[0][:i])+list(params_hq[0][i+1:]),list(params_hq[1][:i])+list(params_hq[1][i+1:]),list(params_hq[2][:i])+list(params_hq[2][i+1:])]
        Rm1s_hq=list(Rm1s_hq[:i])+list(Rm1s_hq[i+1:])
        Rm2s_hq=list(Rm2s_hq[:i])+list(Rm2s_hq[i+1:])
        #if ii==len(j)-1: print("No real solution for Rm_hq", len(j), "times")
    params_hq=np.array(params_hq)
    K_hq=np.mean(1/params_hq[0])
    dK_hq=np.std(1/params_hq[0])
    Kte_hq=np.mean(1/params_hq[1])
    dKte_hq=np.std(1/params_hq[1])
    Qc2_hq=np.mean(1/params_hq[2])
    dQc2_hq=np.std(1/params_hq[2])
    Rm1_hq=np.mean(Rm1s_hq)
    dRm1_hq=np.std(Rm1s_hq)
    Rm2_hq=np.mean(Rm2s_hq)
    dRm2_hq=np.std(Rm2s_hq)
    Kd1_hq=np.mean(np.divide(Rm1s_hq,params_hq[0]))
    Kd2_hq=np.mean(np.divide(Rm2s_hq,params_hq[0]))
    dKd1_hq=np.std(np.divide(Rm1s_hq,params_hq[0]))
    dKd2_hq=np.std(np.divide(Rm2s_hq,params_hq[0]))
    Kt1_hq=np.mean(np.divide(np.add(Rm1s_hq,1),params_hq[1]))
    dKt1_hq=np.std(np.divide(np.add(Rm1s_hq,1),params_hq[1]))
    Kt2_hq=np.mean(np.divide(np.add(Rm2s_hq,1),params_hq[1]))
    dKt2_hq=np.std(np.divide(np.add(Rm2s_hq,1),params_hq[1]))
    return K_hq,dK_hq,Kte_hq,dKte_hq,   Qc2_hq,dQc2_hq,Rm1_hq,dRm1_hq,   Rm2_hq,dRm2_hq,Kd1_hq,dKd1_hq,   Kd2_hq,dKd2_hq,Kt1_hq,dKt1_hq,   Kt2_hq,dKt2_hq,params_hq
#----------------------------------------------------------------------
def UseParams_nqll(params_nqll):
    K_nqll=np.mean(1/params_nqll[0])
    dK_nqll=np.std(1/params_nqll[0])
    Qc2_nqll=np.mean(1/params_nqll[1])
    dQc2_nqll=np.std(1/params_nqll[1]) 
    return K_nqll,dK_nqll,Qc2_nqll,dQc2_nqll
#----------------------------------------------------------------------
def UseParams_nqL(params_nqL):
    Kt_nqL=np.mean(1/params_nqL[0])
    dKt_nqL=np.std(1/params_nqL[0])
    Kw_nqL=np.mean(params_nqL[1])
    dKw_nqL=np.std(params_nqL[1])  
    return Kt_nqL,dKt_nqL,Kw_nqL,dKw_nqL
#----------------------------------------------------------------------
def UseParams_nllL(params_nllL):
    K_nllL=np.mean(1/params_nllL[0])
    dK_nllL=np.std(1/params_nllL[0])
    Kt_nllL=np.mean(1/params_nllL[1])
    dKt_nllL=np.std(1/params_nllL[1])
    Rm_nllL=np.mean(params_nllL[2])
    dRm_nllL=np.std(params_nllL[2])
    Kw_nllL=np.mean(params_nllL[3])
    dKw_nllL=np.std(params_nllL[3])
    Kd_nllL=np.mean(np.multiply(1/params_nllL[0],params_nllL[2]))
    dKd_nllL=np.std(np.multiply(1/params_nllL[0],params_nllL[2]))
    Qc2_nllL=np.mean( np.multiply( np.multiply( params_nllL[0],   1/params_nllL[1]),
                                   np.multiply( 1/params_nllL[2], 1/params_nllL[2]) ) )*4
    dQc2_nllL=np.std( np.multiply( np.multiply( params_nllL[0],   1/params_nllL[1]),
                                   np.multiply( 1/params_nllL[2], 1/params_nllL[2]) ) )*4
    return K_nllL,dK_nllL,   Kt_nllL,dKt_nllL,Rm_nllL,dRm_nllL,   Kw_nllL,dKw_nllL,Kd_nllL,dKd_nllL
#----------------------------------------------------------------------
def UseParams_hll(Nboot,params_hll):
    j=[] #list
    a=np.divide(params_hll[1],np.multiply(np.multiply(params_hll[0],params_hll[2]),4))
    Rm1s_hll=np.zeros(Nboot)
    Rm2s_hll=np.zeros(Nboot)
    for i in range(Nboot):
        if np.imag(cmath.sqrt(1+4*a[i]))==0: Rm1s_hll[i],Rm2s_hll[i]=np.roots([a[i],-1,-1])
        else: j.append(i)
    for ii in range(len(j)):
        i=j[len(j)-ii-1]
        params_hll=[list(params_hll[0][:i])+list(params_hll[0][i+1:]),list(params_hll[1][:i])+list(params_hll[1][i+1:]),list(params_hll[2][:i])+list(params_hll[2][i+1:])]
        Rm1s_hll=list(Rm1s_hll[:i])+list(Rm1s_hll[i+1:])
        Rm2s_hll=list(Rm2s_hll[:i])+list(Rm2s_hll[i+1:])
        if ii==len(j)-1: print("No real solution for Rm_hll", len(j), "times")
    params_hll=np.array(params_hll)
    K_hll=np.mean(1/params_hll[0])
    dK_hll=np.std(1/params_hll[0])
    Kte_hll=np.mean(1/params_hll[1])
    dKte_hll=np.std(1/params_hll[1])
    Qc2_hll=np.mean(1/params_hll[2])
    dQc2_hll=np.std(1/params_hll[2])
    Rm1_hll=np.mean(Rm1s_hll)
    dRm1_hll=np.std(Rm1s_hll)
    Rm2_hll=np.mean(Rm2s_hll)
    dRm2_hll=np.std(Rm2s_hll)    
    Kd1_hll=np.mean(np.divide(Rm1s_hll,params_hll[0]))
    dKd1_hll=np.std(np.divide(Rm1s_hll,params_hll[0]))
    Kd2_hll=np.mean(np.divide(Rm2s_hll,params_hll[0]))
    dKd2_hll=np.std(np.divide(Rm2s_hll,params_hll[0]))
    Kt1_hll=np.mean(np.divide(np.add(Rm1s_hll,1),params_hll[1]))
    dKt1_hll=np.std(np.divide(np.add(Rm1s_hll,1),params_hll[1]))
    Kt2_hll=np.mean(np.divide(np.add(Rm2s_hll,1),params_hll[1]))
    dKt2_hll=np.std(np.divide(np.add(Rm2s_hll,1),params_hll[1]))
    return K_hll,dK_hll,Kte_hll,dKte_hll,   Qc2_hll,dQc2_hll,Rm1_hll,dRm1_hll,   Rm2_hll,dRm2_hll,Kd1_hll,dKd1_hll,   Kd2_hll,dKd2_hll,Kt1_hll,dKt1_hll,   Kt2_hll,dKt2_hll
#----------------------------------------------------------------------
def UseParams_hllL(params_hllL):
    K_hllL=np.mean(1/params_hllL[0])
    dK_hllL=np.std(1/params_hllL[0])
    Kt_hllL=np.mean(1/params_hllL[1])
    dKt_hllL=np.std(1/params_hllL[1])
    Rm_hllL=np.mean(params_hllL[2])
    dRm_hllL=np.std(params_hllL[2])
    Kw_hllL=np.mean(params_hllL[3])
    dKw_hllL=np.std(params_hllL[3])
    return K_hllL,dK_hllL,Kt_hllL,dKt_hllL,Rm_hllL,dRm_hllL,Kw_hllL,dKw_hllL
#----------------------------------------------------------------------
def OtherInfo(Nboot,values,qzdzi):
    return np.array([ np.concatenate((np.array([Nboot]),
                                      np.array([len(qzdzi[0])]),
                                      np.array([len(values)]),
                                      np.concatenate((qzdzi), axis=0),
                                      np.array(values),
                                      np.zeros(Nboot-len(qzdzi[0])*3-3-len(values))
                                     ))])
#----------------------------------------------------------------------

results=np.concatenate((params_hq,   OtherInfo(Nboot,UseParams_hq(Nboot,params_hq),  qzdz1),
              params_nqll, OtherInfo(Nboot,UseParams_nqll(params_nqll),    qzdz2),
              params_nqL,  OtherInfo(Nboot,UseParams_nqL(params_nqL),      qzdz3),
              params_nllL, OtherInfo(Nboot,UseParams_nllL(params_nllL),    qzdz4),
              params_hll,  OtherInfo(Nboot,UseParams_hll(Nboot,params_hll),qzdz5)))
# np.save('surfaces_'+lipid_name,results)
