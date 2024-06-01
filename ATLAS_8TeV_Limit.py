import ROOT
import os,argparse
from ATLAS_table import *
from tabulate import tabulate


pjoin = os.path.join
##----------------------------

parser =argparse.ArgumentParser()
parser.add_argument('--filename'   ,type=str  ,required=True               ,help='histogram file location'              )
parser.add_argument('--mass'       ,type=int  ,required=True               ,help='mass value of the signal'             )
parser.add_argument('--printTable' ,type=bool ,required=False,default=False,help='print yield in nice table'            )
parser.add_argument('--SR'         ,type=str  ,required=True ,default=1    ,help='SR between 1 to 20'                   )
parser.add_argument('--sfactor'    ,type=int  ,required=False,default=1    ,help='multiply signal yield for limit tools')
parser.add_argument('--mode'       ,type=str  ,required=False,default='nom',help='mode on signal (nominal, up, down)'   )

args = parser.parse_args()
#---------------------------------
#####################################################
# VLL XSEC and ATLAS DATA
#####################################################
## M100, M150, M200, M250

VLL_8TeV_NLO_Xsec = {
"M100":0.577584,#pb
"M150":0.128,
"M200":0.0414,
"M250":0.016,    
}

ATLAS_DATA_LUMI = 20.3*1000 #in pb

#----------------------------------------------------
## function to get no of events from the histogram
def get_nevents(h,binlo,binhi,mode=0):
    binWidth = h.GetBinWidth(1)
    ibin_begin = 0
    nevents = 0
    while(h.GetBinCenter(ibin_begin)< binlo):
        ibin_begin+=1
    
    ibin_end = ibin_begin
    while(h.GetBinCenter(ibin_end)<binhi):
        ibin_end+=1
    
    for i in range(ibin_begin,ibin_end):
        #print(i,ibin_begin,ibin_end)
        if(mode==0):nevents += h.GetBinContent(i)
        elif(mode==1):nevents += (h.GetBinContent(i)+h.GetBinError(i))
        elif(mode==-1):nevents += (h.GetBinContent(i)-h.GetBinError(i))
        
        #print(h.GetNbinsX())
        if(i>h.GetNbinsX()+1):break  ##put high no as an end value leads to high ibin_end> total bin of an
                                     ##histogram. So, just terminate if ibin_end>total bins of hist
                                     
    return ((nevents*ATLAS_DATA_LUMI*VLL_8TeV_NLO_Xsec[f'M{args.mass}'])/200000)    

#print signal table
def printSignalTable(info):        
    table=[]
    for k,v in info.items():
        table.append([k]+v)
    
    print(tabulate(table,headers=['SRName','Sel']+Category,tablefmt='grid',numalign='left',floatfmt='.2f'))

#print bkg,dbkg and obs table
def printPaperTable(info):
    table=[]
    for i in range(1,21):
        cat_yield=[]
        cat_yield.append(f'SR{i}')
        for cat in Category:
            for item in info:
                a=cat+'_'+'SR'+str(i)
                b=item[0]+'_'+item[1]
                if a==b :
                    #print(a,b)
                    #cat_yield.append((f'{item[2]}+{item[3][0]}-{item[3][1]}',item[4]))
                    cat_yield.append((f'{item[2]}',item[4]))
                    
        ##
        #print(cat_yield)
        table.append(cat_yield)
        
    print(tabulate(table,headers=['SRName']+Category,tablefmt='grid',numalign='left',floatfmt='.2f'))                  


##----------------------------------------------------

print(args.filename)
#open the histogram rootfile
rf = ROOT.TFile.Open(args.filename,'READ')

## SR category and variables
Category = ['3L_onZ','3L_offZ',' 3L_noOOSF','2L1T_onZ','2L1T_offZ','2L1T_noOOSF']
Variable = ['HTlep','pTmin','MET_HTb150','MET_HTa150','Meff','Meff_METa100']

## make a list of histogram names that are needed
histNames = []
for cat in Category:
    for var in Variable:
        histNames.append(f'{cat}_{var}')

##make a dict containing all the dictionary
rootHist={}
for hname in histNames:
    h = rf.Get(hname)
    rootHist[hname] = rf.Get(hname)
    #print(hname)
    #print(f"Hist Integral: {hname} = {h.Integral()}")

#------------
## Get the signal yield
def getSignalTable(SignalRegion,mode_):
    signal={}
    for SR in SignalRegion:
        NAME,DESC,BINLOW,BINHI,VAR = SR
        SignalYield=[]
        SignalYield.append(DESC)
    
        for cat in Category:                
            signal[NAME]=[]
            nSig = get_nevents(rootHist[cat+'_'+VAR],BINLOW,BINHI,mode_)
            if(nSig<0.01):nSig=0.01
            SignalYield.append(nSig)

        ##    
        signal[NAME]=SignalYield

    return signal   
   
## create a dictionary of s,b,db and obs
def get_yield(signalTable,paperTable,nReg):
    yield_table={}
    for i in range(1,nReg+1):
        yield_table['SR'+str(i)]={}
        s=[]
        b=[]
        db=[]
        obs=[]
        for cat in Category:
            for item in paperTable:
                a_=cat+'_'+'SR'+str(i)
                b_=item[0]+'_'+item[1]
                if a_==b_:
                    s=signalTable['SR'+str(i)][1:]
                    b.append(item[2])
                    db.append(item[3][0])
                    obs.append(item[4])
                    
        yield_table['SR'+str(i)]['s']=s
        yield_table['SR'+str(i)]['b']=b
        yield_table['SR'+str(i)]['db']=db
        yield_table['SR'+str(i)]['obs']=obs
     
    #
    return yield_table
       
##---------------------------------------------------------
#               Statistical Combination    
##----------------------------------------------------------
from ROOT import RooFit, RooStats, RooWorkspace, RooArgSet, RooRealVar
from array import array

def combLim(region,mass):
    
    # Define the signal, the expected background, observation, and background efficiency    
    signalFactor = args.sfactor
    
    s        = array('d',[i*signalFactor for i in region['s']])
    mainMeas = array('d',region['obs'])   
    bkgMeas  = array('d',region['b'])    
    dbMeas   = array('d',region['db'])
    
   
    
    print("s="  ,s)
    print("b="  ,mainMeas)
    print("db=" ,bkgMeas)
    print("obs=",dbMeas)
    
    for i in range(6):
        if(bkgMeas[i]!=0):
            dbMeas[i] = dbMeas[i] / bkgMeas[i]
        else:
            dbMeas[i]=dbMeas[i]
    f = ROOT.RooStats.NumberCountingPdfFactory()
    wspace = RooWorkspace()
    f.AddModel(s, 6, wspace, "TopLevelPdf", "masterSignal")
    f.AddData(mainMeas, bkgMeas, dbMeas, 6, wspace, "ObservedNumberCountingData")

    # Uncomment to see structure of workspace
    # wspace.Print()

    # Define the null hypothesis for the calculator
    mu = wspace.var("masterSignal")
    poi = RooArgSet(mu)
    nullParams = RooArgSet("nullParams")
    nullParams.addClone(mu)
    nullParams.setRealValue("masterSignal", 0)

    # Create a calculator for doing the hypothesis test
    plc = RooStats.ProfileLikelihoodCalculator(wspace.data("ObservedNumberCountingData"), wspace.pdf("TopLevelPdf"), poi, 0.05, nullParams)

    # Get confidence interval
    paramsOfInterest = nullParams
    plc.SetParameters(paramsOfInterest)
    lrint = plc.GetInterval()
    lrint.SetConfidenceLevel(0.95)

    #Get upper and lower limits
    print("lower limit on master signal =", lrint.LowerLimit(mu))
    print("upper limit on master signal =", lrint.UpperLimit(mu))
    print("signal factor                =", signalFactor)
    results ={}
    
    results[f'{mass}']=[lrint.UpperLimit(mu)*signalFactor,lrint.UpperLimit(mu)*signalFactor*VLL_8TeV_NLO_Xsec[f'M{args.mass}'],VLL_8TeV_NLO_Xsec[f'M{args.mass}']]

    ## results[]=[r-value,observed xsec, theory xsec]

    return results


##------run comblim---

#getTable
ATLAS_8TeV_paperTable = HTlep_table + pTmin_table + MET_HTb150_table + MET_HTa150_table + Meff_table + Meff_METa100_table

limit_result={}

uncmode_=[]
uncmode_.append(args.mode)

for uncmode in uncmode_:
    
    if(uncmode=='nom'):mode_=0;
    if(uncmode=='up'     ):mode_=1;
    if(uncmode=='down'   ):mode_=-1;
    
    signalTable = getSignalTable(SignalRegion,mode_)
    yield_table = get_yield(signalTable,ATLAS_8TeV_paperTable,20)
    
    
    #for i in yield_table.keys():
    #    print()
    #    print("calculating limit for {i}.....")
    #    print()
    #    limit_result[i+'_'+uncmode] = combLim(yield_table[i],args.mass)

    ##do it individually for converging the limit setting tools

    region = 'SR'+args.SR
    
    limit_result[region+'_'+uncmode+'__'+str(args.mass)] = combLim(yield_table[region],args.mass)
    

    #--print in a nice table
    if(args.printTable):
        printSignalTable(signalTable)
        printPaperTable(ATLAS_8TeV_paperTable)#6

    
##
print()
print("SIGNAL= ",args.filename.split('hst_')[1])
print("MASS  = ",args.mass)

print(limit_result)



