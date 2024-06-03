import ROOT
import os,json,datetime
import ATLAS_table as ATLAS

##Make a handy dictionary for SR info
ATLAS_SR={}
for item in ATLAS.SignalRegion:    
    ATLAS_SR[item[0]]=item[1]

#print(ATLAS_SR)
def get_timestamp():
    #-------------------------------
    datetime_object = datetime.datetime.now()
    datetime_object = datetime.datetime.strptime(str(datetime_object), "%Y-%m-%d %H:%M:%S.%f")
    timestamp=str(datetime_object.year)+str(("%02d" % (datetime_object.month,)))+str(("%02d" % (datetime_object.day)))
    timestamp+='_'+str(("%02d" % (datetime_object.hour,)))+str(("%02d" % (datetime_object.minute,)))+str(("%02d" % (datetime_object.second))) 
    
    return timestamp
    

timestamp = get_timestamp()

#--------------------------------------------

def get_limits(filename,mass,SR):
    with open(filename,'r') as infile:
        limit_results = json.load(infile)
            
        limits=[]
        median   = limit_results[f'{SR}_nom__Mass{mass}'][f'{mass}'][1]
        theory   = limit_results[f'{SR}_nom__Mass{mass}'][f'{mass}'][2]
        uncUp    = limit_results[f'{SR}_up__Mass{mass}'][f'{mass}'][1]
        uncDown  = limit_results[f'{SR}_down__Mass{mass}'][f'{mass}'][1]        

        limits.extend([median,uncUp,uncDown,theory])
        
        return limits
        
        
def get_UpperLimit(values,SR,model):
    
    labels=[f'M{i}' for i in values]
    
    
    N = len(labels)
    unc     = ROOT.TGraph(2*N)     # unc
    median  = ROOT.TGraph(N)       # median
    theory  = ROOT.TGraph(N)       # theory
        
    
    upperLimit=[]
    
    for i in range(N):
        
        mass = values[i]
        
        limit = get_limits(f'LimitOutput/{model}_M{mass}.json',mass,SR)
        #print(mass, limit)
        
        median.SetPoint(    i,    values[i], limit[0]) # median
        unc.SetPoint(       i,    values[i], limit[1]) # s+ds
        unc.SetPoint( 2*N-1-i,    values[i], limit[2]) # s-ds
        theory.SetPoint(    i,    values[i], limit[3]) #theory
        
        
    upperLimit.extend([median,unc,theory])
    
    
    return upperLimit


def get_color(SR):
    #SR1-SR3  = Htlep == kRed
    #SR3-SR6  = pTmin == kBlue
    #SR7-SR10 = MET_HTb150 == kGreen
    #SR11-SR14= MET_HTa150 == kGeen
    #SR15-SR17= Meff = kOrange
    #SR18-SR20= Meff_METa100 = kOrange
    
    color_palette ={
        'SR1':ROOT.kRed+1,
        'SR2':ROOT.kRed+2,
        'SR3':ROOT.kRed-4,
        'SR4':ROOT.kBlue,
        'SR5':ROOT.kBlue+2,
        'SR6':ROOT.kBlue-4,
        'SR7':ROOT.kGreen+1,
        'SR8':ROOT.kGreen+2,
        'SR9':ROOT.kGreen+3,
        'SR10':ROOT.kGreen-2,
        'SR11':ROOT.kViolet+1,
        'SR12':ROOT.kViolet+2,
        'SR13':ROOT.kViolet+3,
        'SR14':ROOT.kViolet-9,
        'SR15':ROOT.kOrange,
        'SR16':ROOT.kOrange+1,
        'SR17':ROOT.kOrange+2,
        'SR18':ROOT.kOrange-5,
        'SR19':ROOT.kOrange-6,
        'SR20':ROOT.kOrange-7,                
    }
    
    return color_palette[SR]


def plot_limit(limit,model):    
    
    #Canvas
    W = 1000
    H = 800
    T = 0.08*H
    B = 0.12*H
    L = 0.15*W
    R = 0.25*W
    c = ROOT.TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetLogy(1)
    c.cd()
    
    ##frame
    frame = c.DrawFrame(1.4,0.001, 3.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetLabelSize(0.03)
    frame.GetXaxis().SetTitleOffset(1.2)
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.GetXaxis().SetNdivisions(510)
    #frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("#sigma (pb)")
    frame.GetXaxis().SetTitle("m_{VLL} (GeV)")
    frame.SetMinimum(0.005)
    #frame.SetMaximum(max(up2s)*1.05)
    frame.SetMaximum(1000)
    frame.GetXaxis().SetLimits(50,300)
    
    
    x1 = 0.75
    x2 = x1 + 0.24
    y2 = 0.85
    y1 = 0.35
    legend =ROOT.TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0);legend.SetBorderSize(0);legend.SetTextSize(0.02);
    legend.SetTextAlign(12);legend.SetTextFont(42)
    legend.SetHeader(f'{model}: 95% CL')
    #legend.AddEntry(median,"observed",'L')
    #legend.AddEntry(unc   ,"Uncertainty",'f')
    
    for i,SR in enumerate(limit.keys()):
            
        median = limit[SR][0]
        unc    = limit[SR][1]    
    
        ##
        unc.SetFillColor(ROOT.kOrange)
        unc.SetLineColor(ROOT.kOrange)
        unc.SetFillStyle(1001)
    
        median.SetLineColor(get_color(SR))
        median.SetLineWidth(2)
        median.SetLineStyle(1)
        median.SetMarkerStyle(20)
        median.SetMarkerSize(0.6)
        median.SetMarkerColor(i+1)
        
        if(i==0):
            median.Draw('LP')
            #unc.Draw('F')
        else:
            median.Draw('LPSAME')
            #unc.Draw('FSAME')
            
            
        legend.AddEntry(median,f"{SR}:{ATLAS_SR[SR]}",'L')    
    
    theory = limit[SR][2]
    legend.AddEntry(theory,"Theory xsec",'f')
    #theoryxsec
    theory.SetLineColor(ROOT.kRed)
    theory.SetFillColor(ROOT.kRed-9)
    theory.SetLineWidth(2)
    theory.SetLineStyle(1)
    theory.SetMarkerStyle(1)
    theory.Draw('L')
    
    ##
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')
    

    legend.Draw()
    
    
    c.SaveAs(f"{model}_test_limit_{timestamp}.pdf")
    c.SaveAs(f"{model}_test_limit_{timestamp}.png")
    c.Close()


def main(model):
    values=[100,150,200,250]
    SR=[f'SR{i}' for i in range (1,21)]
    SR=['SR11','SR15']
    upperLimit={}
    for sr in SR:
        upperLimit[sr] = get_UpperLimit(values,sr,model)


    print(upperLimit)
    plot_limit(upperLimit,model)


if(__name__=='__main__'):
    main('VLLTau')
    main('VLLEle')
    main('VLLMu')    
