from numpy import *
from rhoratefcns import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
import os
from construction import LoadData

GENERATE = False

MsunV = 4.83
Gconst = 6.67259e-8
realMsun = 1.989e33
Rsun = 6.9599e10
pc = 3.1e18
km = 10**5
yr = 365*24*3600

rtest = arange(-7,7,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

def getrate(model,partial = False):
    Menc,psi,Jc2,g,G,f = 0,1,2,3,4,5
    if partial == False:
        if model.generate == False:
            seton = {Menc:"OFF",psi:"OFF",Jc2:"OFF",g:"OFF",G:"OFF",f:"OFF"}
        if model.generate == True:
            seton = {Menc:"ON",psi:"ON",Jc2:"ON",g:"ON",G:"ON",f:"ON"}
    if partial != False:
        seton = partial    
    try:        
        os.system('echo $DISPLAY > tempdisplay')
        displays = LoadData('tempdisplay')
        os.system('rm -f tempdisplay')
        display = displays[0]
        if display == '':
            plottinglist = {Menc:False,psi:False,Jc2:False,g:False,G:False,f:False}
        if display != '':
            plottinglist = {Menc:['r','M'],psi:['r',r'$\psi$'],Jc2:['E',r'$J_c^2$'],g:['E','g'],G:['E','G'],f:['E','f']}
        
        exps = {Menc:[3-model.g,0],psi:[-1,-1],Jc2:[-1,-1],g:[model.b-0.5,model.g-0.5],G:[model.b-4,model.g-4],f:[model.b-1.5,model.g-1.5]}

        model.statfile.write('GALAXY: {0}\n'.format(model.name))

        model.statfile.write('Menc:\n')
        rarray,rchange,rstart = rgrid([model],4,-6,0.03)
        Mencgood = compute([model],funcMenc,rtest,[4,-6,0.03],rgrid,exp[Menc],plottinglist[Menc],seton[Menc])

        if Mencgood == 0:
            model.statfile.write('Failed to evaluate Menc')
            model.statfile.close()
            model.pdfdump.close()
            return 0,0,0,0,0,0
    
        elif Mencgood != 0:
            
            model.statfile.write('\npsi:\n')
            pprereqs = [model,'Model',Mencgood,'Menc']
            psigood = compute(pprereqs,funcpsi,rtest,[4.3,-6,0.03],rgrid,exp[psi],plottinglist[psi],seton[psi])
    
            if psigood == 0:
                model.statfile.write('Failed to evaluate psi')
                model.statfile.close()
                model.pdfdump.close()
                return Mencgood,0,0,0,0,0

            elif psigood != 0:
                Jprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
                
                model.statfile.write('\nJc2:\n')
                Jc2good = compute(Jprereqs,funcJc2,rtest,[3,-4,0.01],Egrid,exps[Jc2],plottinglist[Jc2],seton[Jc2])

                if Jc2good == 0:
                    model.statfile.write('Failed to evaluate Jc2')
                    model.statfile.close()
                    model.pdfdump.close()
                    return Mencgood,psigood,0,0,0,0
                
                elif Jc2good != 0:

                    lgprereqs = [model,'Model',psigood,"psi"]
                    
                    model.statfile.write('\ng:\n')
                    ggood = compute(lgprereqs,funclg,rtest,[3,-3,0.1],Egrid,exps[g],plottinglist[g],seton[g])
                
                    if ggood == 0:
                        model.statfile.write('Failed to evaluate g')
                        model.statfile.close()
                        model.pdfdump.close()
                        return Mencgood,psigood,Jc2good,0,0,0
                
                    elif ggood != 0:
                        bGprereqs = [model,'Model',psigood, "psi",ggood,"g"]
                        
                        Gtest = arange(-5,5,0.01)
                        Gtest = append(Gtest,40)
                        Gtest = insert(Gtest,0,-40)
                        Gtest = 10**Gtest
                        
                        model.statfile.write('\nG:\n')
                        Ggood = compute(bGprereqs,funcbG,Gtest,[3,-3,0.1],Egrid,exps[G],plottinglist[G],seton[G])
    
                        if model.memo == True:
                            model.p1bG = {}
                            model.p2bG = {}
                            model.p3bG = {}
                    
                        if Ggood == 0:
                            model.statfile.write('Failed to evaluate G')
                            model.statfile.close()
                            model.pdfdump.close()
                            return Mencgood,psigood,Jc2good,ggood,0,0
                        
                        elif Ggood != 0:
                            
                            fprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
                            
                            ftest = arange(-3,5,0.01)
                            ftest = append(ftest,40)
                            ftest = insert(ftest,0,-40)
                            ftest = 10**ftest
                            
                            model.statfile.write('\nf:\n')
                            fgood = compute(fprereqs,funcf,ftest,[5,-3,0.03],Egrid,exp[f],plottinglist[f],seton[f])
                            if fgood == 0:
                                model.statfile.write('Failed to evaluate f')
                                model.statfile.close()
                                model.pdfdump.close()
                                print('\a')
                                return Mencgood,psigood,Jc2good,ggood,Ggood,0
                        
                            elif fgood != 0:
                                model.statfile.close()
                                model.pdfdump.close()
                                print('\a')
                                return Mencgood,psigood,Jc2good,ggood,Ggood,fgood
    except KeyboardInterrupt:
        model.statfile.write('\n\nFunction creation cancelled')
        model.statfile.close()
        model.pdfdump.close()
        raise
        
if __name__ == '__main__':

    alpha = 1.0
    beta = 4.0
    gamma = 1.5
    r0pc = 1.0
    rho0 = 1e5
    MBH_Msun = 1e3
    name = 'testform'

    from rhomodels import NukerModelRho
    model = NukerModelRho(name,alpha,beta,gamma,r0pc,rho0,MBH_Msun,GENERATE,memo = False)
    Mencgood,psigood,Jc2good,ggood,Ggood,fgood = getrateplot(model)




