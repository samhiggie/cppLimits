#ifndef plot_h
#define plot_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <boost/program_options.hpp>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TROOT.h>

#include "time.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooFitResult.h"
#include "RooBernstein.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TString.h"
#include "TList.h"
#include "TAxis.h"
#include "TAxis.h"
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <map>
#include <utility>

TPaveText add_lumi(int intyear,bool doRatio){
    string year;
    year = to_string(intyear);
    float lowX=0.65;
    float lowY=0.835;
    if(!doRatio){
        lowX=0.60;
    }
    TPaveText lumi  = TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC");
    lumi.SetBorderSize(   0 );
    lumi.SetFillStyle(    0 );
    lumi.SetTextAlign(   12 );
    lumi.SetTextColor(    1 );
    lumi.SetTextSize(0.04);
    lumi.SetTextFont (   42 );
    if (year == "2016")
        lumi.AddText(((year)+" 35.9 fb^{-1} (13 TeV)").c_str());
    if (year == "2017")
        lumi.AddText(((year)+" 41.8 fb^{-1} (13 TeV)").c_str());
    if (year == "2018")
        lumi.AddText(((year)+" 59.7 fb^{-1} (13 TeV)").c_str());
    return lumi;
}

TPaveText add_CMS(bool doRatio){
    float lowX=0.17;
    float lowY=0.835;
    if(!doRatio){
        lowX=0.17;
        lowY=0.835;
    }
    TPaveText lumi  = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC");
    lumi.SetTextFont(61);
    lumi.SetTextSize(0.06);
    lumi.SetBorderSize(   0 );
    lumi.SetFillStyle(    0 );
    lumi.SetTextAlign(   12 );
    lumi.SetTextColor(    1 );
    lumi.AddText("CMS");
    return lumi;
}

TPaveText add_Preliminary(string channel, bool doRatio){
    float lowX=0.45;
    float lowY=0.835;
    if(!doRatio){
        lowX=0.30;
    }
    TPaveText lumi  = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC");
    lumi.SetTextFont(52);
    lumi.SetTextSize(0.04);
    lumi.SetBorderSize(   0 );
    lumi.SetFillStyle(    0 );
    lumi.SetTextAlign(   12 );
    lumi.SetTextColor(    1 );
    lumi.AddText(("Preliminary "+channel).c_str());
    return lumi;
}

#endif
