#! /usr/bin/env python
import ROOT as r,sys,math,array

class hist:
    def __init__( self , iVars, iHists):
        self._vals  = iVars
        self._hists = iHists
        self._mass  = r.RooRealVar("x","x",iHists[0].GetXaxis().GetXmin(),iHists[0].GetXaxis().GetXmax())
        self._mass.setBins(iHists[0].GetNbinsX())
        self._dhist = []
        self._hpdf  = []
                
    def histpdf(self,iH,iShift=0):
        lDH    = r.RooDataHist(iH.GetName()+"d0" ,iH.GetName()+"d0" ,r.RooArgList(self._mass),iH)
        if iShift != 0:
            lHPdf  = r.RooHistPdf (iH.GetName()+"dh0",iH.GetName()+"dh0",r.RooArgList(iShift),r.RooArgList(self._mass),lDH)
        else:
            lHPdf  = r.RooHistPdf (iH.GetName()+"dh0",iH.GetName()+"dh0",r.RooArgSet(self._mass),lDH)
        self._dhist.extend([lDH])
        self._hpdf .extend([lHPdf])
        return lHPdf
        
    def morph(self,iValue):
        lMarker=0
        for pVal in self._vals:
            if pVal > iValue:
                break
            lMarker=lMarker+1
        lVal   = iValue-self._vals[lMarker-1]
        lVal  /= float(abs(self._vals[lMarker]-self._vals[lMarker-1]))
        lMed   = r.RooRealVar("med","med",0,1)
        lRH0   = self.histpdf(self._hists[lMarker])
        lRH1   = self.histpdf(self._hists[lMarker-1])
        lMorph = r.RooIntegralMorph("Morph","Morph",lRH0,lRH1,self._mass,lMed)
        lMed.setVal(lVal);
        lOut   = lMorph.createHistogram("tmp"+str(iValue),self._mass)
        lInt   = (self._hists[lMarker].Integral()-self._hists[lMarker-1].Integral())*lVal+self._hists[lMarker-1].Integral()
        lOut.Scale(lInt)
        return lOut

    def shift(self,iH,iScale):
        lDM     = r.RooRealVar("dm","dm", 0.,-10,10)
        lShift  = r.RooFormulaVar("shift","x-dm",r.RooArgList(self._mass,lDM))  
        lHPdf = self.histpdf(iH,lShift)
        lDM.setVal(iScale)
        lUp = lHPdf.createHistogram("x")
        lUp.SetTitle(lHPdf.GetName()+"_scaleUp")
        lUp.SetName (lHPdf.GetName()+"_scaleUp")
        lUp.Scale(iH.Integral())
        
        lDM.setVal(-iScale)
        lDown = lHPdf.createHistogram("x")
        lDown.SetTitle(lHPdf.GetName()+"_scaleDown")
        lDown.SetName (lHPdf.GetName()+"_scaleDown")
        lDown.Scale(iH.Integral())
        return [lUp,lDown]
    
    def smear(self,iH,iScale):
        lDM     = r.RooRealVar("dm","dm", 1.,0.,2.)
        lShift  = r.RooFormulaVar("shift","("+self._mass.GetName()+"-"+str(iH.GetMean())+")/dm+"+str(iH.GetMean()),r.RooArgList(self._mass,lDM))  
        lHPdf = self.histpdf(iH,lShift)
        lDM.setVal(1.+iScale)
        lUp = lHPdf.createHistogram("x")
        lUp.SetTitle(lHPdf.GetName()+"_scaleUp")
        lUp.SetName (lHPdf.GetName()+"_scaleUp")
        lUp.Scale(iH.Integral())
                
        lDM.setVal(1.-iScale)
        lDown = lHPdf.createHistogram("x")
        lDown.SetTitle(lHPdf.GetName()+"_scaleDown")
        lDown.SetName (lHPdf.GetName()+"_scaleDown")
        lDown.Scale(iH.Integral())
        return [lUp,lDown]
