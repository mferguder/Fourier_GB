def Get_lim(xmaxlim,xqtavg,dxqtavg,func,func2,cmax):
    xlim=xmaxlim
    Z,dZ=Get_ZdZ_n(xqtavg,dxqtavg,pows)
    xdatar,Zr,dZr=Get_fit_data(xdata,Z,dZ,xlim)
    xdatar2=np.sqrt(xdatar[0]**2+xdatar[1]**2)
    popt, pcov = curve_fit(func2, xdatar, Zr, p0,sigma=dZr)
    c=stats.chi2.cdf(sum(np.power(np.divide(np.subtract(Zr,func(xdatar[0],xdatar[1],*popt)),dZr),2)) , len(xdatar[0])-len(popt))
    while c>cmax:
        xlim-=1
        xdatar,Zr,dZr=Get_fit_data(xdata,Z,dZ,xlim)
        xdatar2=np.sqrt(xdatar[0]**2+xdatar[1]**2)
        popt, pcov = curve_fit(func2, xdatar, Zr, p0,sigma=dZr)
        c=stats.chi2.cdf(sum(np.power(np.divide(np.subtract(Zr,func(xdatar[0],xdatar[1],*popt)),dZr),2)) , len(xdatar[0])-len(popt))
    print c
    return xlim
def Get_ignore(xdata,lim):
    ignore=[]
    for i in range(len(xdata[0])):
        qx=xdata[0,i]/(2*np.pi/Lavg)
        qy=xdata[1,i]/(2*np.pi/Lavg)
        if qx**2+qy**2>lim or qx<-qy or (qy<1 and qx==-qy): ignore=np.append(ignore,i)
    return ignore
def Get_fit_data(xdata,Z,dZ,lim):
    ignore=Get_ignore(xdata,lim)
    xdatar=np.delete(xdata    , ignore, 1)
    Zr=np.delete(Z.ravel(), ignore, 0)
    dZr=np.delete(dZ.ravel(), ignore, 0)
    return xdatar,Zr,dZr


gg=4
xmin, xmax, nx = -gg,gg,2*gg+1 #-12, 11, 24
ymin, ymax, ny = -gg,gg,2*gg+1 #-12, 11, 24
x, y = np.linspace(xmin, xmax, nx)*(2*np.pi/Lavg), np.linspace(ymin, ymax, ny)*(2*np.pi/Lavg)
X, Y = np.meshgrid(x, y)
xdata = np.vstack((X.ravel(), Y.ravel()))

Nb=20
for i in range(Nb):
    Z,dZ=Get_ZdZ_n_rand(hqtavg,dhqtavg,pows)
    xdatar,Zr,dZr=Get_fit_data(xdata,Z,dZ,lim)
    [popts[0,i],popts[1,i]], pcov = curve_fit(_new_hq, xdatar, Zr, p0,sigma=dZr)
