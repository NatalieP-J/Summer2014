from numpy import *
from rhoratefcns import *
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
import os
from construction import LoadData
from construction import pklread

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

utest = arange(-7,0,0.01)
utest = insert(utest,0,-40)
utest = 10**utest

def displaycheck():
    os.system('echo $DISPLAY > tempdisplay')
    displays = LoadData('tempdisplay')
    os.system('rm -f tempdisplay')
    display = displays[0]
    if display == '':
        return False
    elif display != '':
        return True

def existcheck(directory,dcheck):
    seton = {}
    plottinglist = {}
    gvals = {}
    prereqs = [[],['Menc'],['Menc','psi'],['psi'],['psi','g'],['Menc','psi'],['Jc2','G','f']]
    strnames = ['Menc','psi','Jc2','g','G','f','dgdlnrp']
    pvals = [['r','M'],['r',r'$\psi$'],['E',r'$J_c^2$'],['E','g'],['E','G'],['E','f'],[r'$u^2$',r'$\frac{dg}{dlnr_p}$']]
    for i in range(len(strnames)):
        prechecks = array([])
        for j in range(len(prereqs[i])):
            prechecks = append(prechecks,gvals[prereqs[i][j]])
        prepass = len(prechecks) - len(prechecks[where(prechecks == True)])
        try:
            vals = pklread('{0}/{1}.pkl'.format(directory,strnames[i]))
            gcheck = goodcheck(vals[:,1])
            gvals[strnames[i]] = gcheck
            if gcheck == True and prepass == 0:
                seton[strnames[i]] = 'OFF'
                if dcheck == True:
                    plottinglist[strnames[i]] = pvals[i]
                if dcheck == False:
                    plottinglist[strnames[i]] = False
            if gcheck != True or prepass != 0:
                seton[strnames[i]] = 'FAIL'
                plottinglist[strnames[i]] = False
        except IOError:
            gvals[strnames[i]] = False
            if prepass == 0:
                seton[strnames[i]] = 'ON'
                if dcheck == True:
                    plottinglist[strnames[i]] = pvals[i]
                if dcheck == False:
                    plottinglist[strnames[i]] = False
            if prepass != 0:
                seton[strnames[i]] = 'FAIL'
                plottinglist[strnames[i]] = False
    return seton,plottinglist
            
def getrate(model,partial = False):
    dcheck = displaycheck()
    Menc,psi,Jc2,g,G,f,rate = 0,1,2,3,4,5,6
    if partial == False:
        if model.generate == False:
            seton,plottinglist = existcheck(model.directory,dcheck)
        if model.generate == True:
            seton = {'Menc':"ON",'psi':"ON",'Jc2':"ON",'g':"ON",'G':"ON",'f':"ON",'dgdlnrp':"ON"}
            if dcheck == False:
                plottinglist = {'Menc':False,'psi':False,'Jc2':False,'g':False,'G':False,'f':False,'dgdlnrp':False}
            if dcheck == True:
                plottinglist = {'Menc':['r','M'],'psi':['r',r'$\psi$'],'Jc2':['E',r'$J_c^2$'],'g':['E','g'],'G':['E','G'],'f':['E','f'],'dgdlnrp':[r'$u^2$',r'$\frac{dg}{dlnr_p}$']}
    elif partial != False:
        seton = partial 
        plottinglist = {'Menc':False,'psi':False,'Jc2':False,'g':False,'G':False,'f':False,'dgdlnrp':False}   
    try:                
        exps = {'Menc':[3-model.g,0],'psi':[-1,-1],'Jc2':[-1,-1],'g':[model.b-0.5,model.g-0.5],'G':[model.b-4,model.g-4],'f':[model.b-1.5,model.g-1.5],'dgdlnrp':[2,0]}

        sh = {'Menc':[4,-6,0.03],'psi':[4.3,-6,0.03],'Jc2':[3,-4,0.01],'g':[3,-3,0.1],'G':[3,-3,0.1],'f':[5,-3,0.03],'rate':[0,-4,0.04]}

        model.statfile.write('GALAXY: {0}\n'.format(model.name))

        model.statfile.write('Menc:\n')
        up,down,step = sh[Menc]
        rarray,rchange,rstart = rgrid([model],up,down,step)
        Mencgood = compute([model],funcMenc,rtest,sh['Menc'],rgrid,exps['Menc'],plottinglist['Menc'],seton['Menc'])

        if Mencgood == 0:
            model.statfile.write('Failed to evaluate Menc')
            model.statfile.close()
            model.pdfdump.close()
            return 0,0,0,0,0,0,0
    
        elif Mencgood != 0:
            
            model.statfile.write('\npsi:\n')
            pprereqs = [model,'Model',Mencgood,'Menc']
            psigood = compute(pprereqs,funcpsi,rtest,sh['psi'],rgrid,exps['psi'],plottinglist['psi'],seton['psi'])
    
            if psigood == 0:
                model.statfile.write('Failed to evaluate psi')
                model.statfile.close()
                model.pdfdump.close()
                return Mencgood,0,0,0,0,0,0

            elif psigood != 0:
                Jprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
                
                model.statfile.write('\nJc2:\n')
                Jc2good = compute(Jprereqs,funcJc2,rtest,sh['Jc2'],Egrid,exps['Jc2'],plottinglist['Jc2'],seton['Jc2'])

                if Jc2good == 0:
                    model.statfile.write('Failed to evaluate Jc2')
                    model.statfile.close()
                    model.pdfdump.close()
                    return Mencgood,psigood,0,0,0,0,0
                
                elif Jc2good != 0:

                    lgprereqs = [model,'Model',psigood,"psi"]
                    
                    model.statfile.write('\ng:\n')
                    ggood = compute(lgprereqs,funclg,rtest,sh['g'],Egrid,exps['g'],plottinglist['g'],seton['g'])
                
                    if ggood == 0:
                        model.statfile.write('Failed to evaluate g')
                        model.statfile.close()
                        model.pdfdump.close()
                        return Mencgood,psigood,Jc2good,0,0,0,0
                
                    elif ggood != 0:
                        bGprereqs = [model,'Model',psigood, "psi",ggood,"g"]
                        
                        Gtest = arange(-5,5,0.01)
                        Gtest = append(Gtest,40)
                        Gtest = insert(Gtest,0,-40)
                        Gtest = 10**Gtest
                        
                        model.statfile.write('\nG:\n')
                        Ggood = compute(bGprereqs,funcbG,Gtest,sh['G'],Egrid,exps['G'],plottinglist['G'],seton['G'])
    
                        if model.memo == True:
                            model.p1bG = {}
                            model.p2bG = {}
                            model.p3bG = {}
                    
                        if Ggood == 0:
                            model.statfile.write('Failed to evaluate G')
                            model.statfile.close()
                            model.pdfdump.close()
                            return Mencgood,psigood,Jc2good,ggood,0,0,0
                        
                        elif Ggood != 0:
                            
                            fprereqs = [model,'Model',Mencgood,"Menc",psigood,"psi"]
                            
                            ftest = arange(-3,5,0.01)
                            ftest = append(ftest,40)
                            ftest = insert(ftest,0,-40)
                            ftest = 10**ftest
                            
                            model.statfile.write('\nf:\n')
                            fgood = compute(fprereqs,funcf,ftest,sh['f'],Egrid,exps['f'],plottinglist['f'],seton['f'])
                            if fgood == 0:

                                model.statfile.write('Failed to evaluate f')
                                model.statfile.close()
                                model.pdfdump.close()

                                return Mencgood,psigood,Jc2good,ggood,Ggood,0,0
                        
                            elif fgood != 0:

                                rprereqs = [model,'Model',Jc2good,'Jc2',Ggood,'Ggood',fgood,'fgood']
                                model.statfile.write('\nrate:\n')
                                rategood = compute(rprereqs,funcdgdlnrp,utest,sh['dgdlnrp'],stdgrid,exps['dgdlnrp'],plottinglist['dgdlnrp'],seton['dgdlnrp'])
                                
                                if rategood == 0:
                                    model.statfile.close()
                                    model.pdfdump.close()
                                    print('\a')
                                    return Mencgood,psigood,Jc2good,ggood,Ggood,fgood,0

                                elif rategood != 0:
                                    model.statfile.close()
                                    model.pdfdump.close()
                                    if plottinglist == {Menc:['r','M'],psi:['r',r'$\psi$'],Jc2:['E',r'$J_c^2$'],g:['E','g'],G:['E','G'],f:['E','f'],rate:[r'$u^2$',r'$\frac{dg}{dlnr_p}$']}:
                                        os.system('mv {0}/{1}_master.pdf {0}/{1}_complete.pdf'.format(model.directory,model.name))
                                    print('\a')
                                    return Mencgood,psigood,Jc2good,ggood,Ggood,fgood,rategood
                
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





