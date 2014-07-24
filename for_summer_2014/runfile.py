
rtest = arange(-7,7,0.01)
rtest = append(rtest,40)
rtest = insert(rtest,0,-40)
rtest = 10**rtest

Mencgood = compute([],["Menc",Menc],funcMenc,rtest,[4,-6,0.03],rgrid,[3-model.g,0],[[2,0,3-model.b,4*pi*model.rho(rchange)*(rchange**3)],['r','M'],False])

psigood = compute([],["psi",psi],funcpsi,rtest,[4.3,-6,0.03],rgrid,[-1,-1],[False,['r','$\psi$'],False])

Jprereqs = [Mencgood,"Menc",psigood,"psi"]

Jc2good = compute(Jprereqs,["Jc2",Jc2],funcJc2,rtest,[3,-4,0.01],Egrid,[-1,-1],[False,['E','Jc2'],False])

lgprereqs = [psigood,"psi"]

ggood = compute(lgprereqs,["g",g],funclg,rtest,[3,-3,0.1],Egrid,[model.b-0.5,model.g-0.5],[False,['E','g'],False])

bGprereqs = [psigood, "psi",ggood,"g"]

Gtest = arange(-5,5,0.01)
Gtest = append(Gtest,40)
Gtest = insert(Gtest,0,-40)
Gtest = 10**Gtest

Ggood = compute(bGprereqs,["G",G],funcbG,Gtest,[3,-3,0.1],Egrid,[model.b-4,model.g-4],[False,['E','G'],False])

psibG_memo = {}
part2bG_memo = {}
part3bG_memo = {}

fprereqs = [Mencgood,"Menc",psigood,"psi"]

ftest = arange(-3,5,0.01)
ftest = append(ftest,40)
ftest = insert(ftest,0,-40)
ftest = 10**ftest

fgood = compute(fprereqs,["f",f],funcf,ftest,[5,-3,0.03],Egrid,[model.b-1.5,model.g-1.5],[False,['E','f'],False])
