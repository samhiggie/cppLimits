#ifndef Office_h
#define Office_h

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
#include <boost/program_options.hpp>

#include "LimitShape.h"
#include "plot.h"

namespace po = boost::program_options;
using namespace std;
using namespace ROOT;
using namespace RooFit;

/////////////////
//unpacking by vector 
//https://stackoverflow.com/questions/11044504/any-solution-to-unpack-a-vector-to-function-arguments-in-c
////////////////
#include <iostream>
#include <utility>
#include <vector>
#include <cassert>

//https://stackoverflow.com/questions/25925290/c-round-a-double-up-to-2-decimal-places
double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}

namespace util {
template <typename ReturnType, typename... Args>
struct function_traits_defs {
  static constexpr size_t arity = sizeof...(Args);

  using result_type = ReturnType;

  template <size_t i>
  struct arg {
    using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
  };
};

template <typename T>
struct function_traits_impl;

template <typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(Args...)>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(*)(Args...)>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...)>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const&&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) volatile>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) volatile&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) volatile&&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const volatile>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const volatile&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const volatile&&>
    : function_traits_defs<ReturnType, Args...> {};

template <typename T, typename V = void>
struct function_traits
    : function_traits_impl<T> {};

template <typename T>
struct function_traits<T, decltype((void)&T::operator())>
    : function_traits_impl<decltype(&T::operator())> {};

template <size_t... Indices>
struct indices {
  using next = indices<Indices..., sizeof...(Indices)>;
};
template <size_t N>
struct build_indices {
  using type = typename build_indices<N - 1>::type::next;
};
template <>
struct build_indices<0> {
  using type = indices<>;
};
template <size_t N>
using BuildIndices = typename build_indices<N>::type;

namespace details {
template <typename FuncType,
          typename VecType,
          size_t... I,
          typename Traits = function_traits<FuncType>,
          typename ReturnT = typename Traits::result_type>
ReturnT do_call(FuncType& func,
                VecType& args,
           indices<I...> ) {
  assert(args.size() >= Traits::arity);
  return func(args[I]...);
}
}  // namespace details

template <typename FuncType,
          typename VecType,
          typename Traits = function_traits<FuncType>,
          typename ReturnT = typename Traits::result_type>
ReturnT unpack_caller(FuncType& func,
                VecType& args) {
  return details::do_call(func, args, BuildIndices<Traits::arity>());
}
}  // namespace util

float findPull(RooDataSet* nominal,RooDataSet* up,RooDataSet* down){
    cout<<"roodata set for pull "<<nominal->GetName()<<endl;
    float val_nom = float(nominal->sumEntries());
    float val_up = float(up->sumEntries());
    float val_down = float(down->sumEntries());
    float s1 = float(val_nom - val_up)/float(val_nom);
    float s2 = float(val_nom - val_down)/float(val_nom);
    return float(1.0000 + sqrt(s1*s1+s2*s2));
}

class Office{
    public:
    //Data members
    po::variables_map shape_vm;
    string           output="default";//string for naming
    string           outputdir="default_folder";
    string           channel="mmmt";
    int           year=2017;
    RooWorkspace * wsp; 
    ofstream txtfile;
    map<string,ofstream> txtfiles;
    vector<string> systematics = {
        "scale_e","scale_m_etalt1p2",
                       "scale_m_eta1p2to2p1","scale_m_etagt2p1",
                       "scale_t_1prong","scale_t_1prong1pizero",
                       "scale_t_3prong","scale_t_3prong1pizero"
    };
    //vector<shape*> signalLimitShapes;
    map<string,LimitShape*> shapes;


    //functions
    Office(po::variables_map shape_vm);
    ~Office();
//void loadSignalShapes(vector<LimitShape*> inputshapes);
    void loadShape(string name, LimitShape * shape);
    void interpolateParameters(string type, map<string,LimitShape*>& inputshapes,string systematic_name);
    void createTxtfile();
    void printDatacard(bool isrealdata);
    void createTxtfilePerMass();
    void printDatacardPerMass(bool isrealdata);
    void setFunctionName(RooAbsPdf * function,string name);


};
Office::Office(po::variables_map shape_vm){
    output = shape_vm["output-file"].as<std::string>();
    outputdir = shape_vm["output-dir"].as<std::string>();
    channel = shape_vm["channel"].as<std::string>();
    year = stoi(output.substr(0,4));
    wsp = new RooWorkspace("w");
    
}

Office::~Office(){}

void Office::loadShape(string name, LimitShape * shape){
    shapes[name]=shape;
    return;
}
void Office::interpolateParameters(string type, map<string,LimitShape*>& inputshapes,string systematic_name){
    if (type.compare("gaussian")==0){
    cout<<"interpolating!"<<endl;
    TF1 * meanfit = new TF1("meanfit","pol1",16,66);
    TF1 * normfit = new TF1("normfit","pol3",16,66);
    TF1 * sigmafit = new TF1("sigmafit","pol3",16,66);
    TGraphErrors * meangraph = new TGraphErrors();
    TGraphErrors * normgraph = new TGraphErrors();
    TGraphErrors * sigmagraph = new TGraphErrors();
    int c = 0;
    //LimitShape * sp = (LimitShape*) inputshapes["a20"];
    cout<<"shape exist? name: "<<inputshapes["a20"]<<endl;
    cout<<"shape exist? name: "<<inputshapes["a20"]->shape_name<<endl;
    for(auto const& x: inputshapes){
        string mass = x.first.substr(1,2); 
        LimitShape * sp = x.second;
        //cout<<"setting points "<<x.first<<endl; 
        //cout<<"point "<<to_string(c)<<endl;
        //cout<<"Mass "<<(mass)<<endl;
        //cout<<"shape name "<<sp->shape_name<<endl;
        //cout<<"var name "<<sp->coeffs["gaussian"+mass]->At(0)->GetName()<<endl;
        //cout<<"value "<<((RooRealVar*)sp->coeffs["gaussian"+mass]->At(0))->getVal()<<endl;

        meangraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["gaussian"+mass]->At(0))->getVal()
        );
        sigmagraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["gaussian"+mass]->At(1))->getVal()
        );
        normgraph->SetPoint(c,
            stof((mass)),
            sp->data->sumEntries()
        );

        c++;
    }
    TCanvas * can = new TCanvas("splinefits","splinefits",600,600);
    can->cd();
    gStyle->SetOptStat(1); 
    gStyle->SetOptFit(1); 
    gStyle->SetStatX(0.6);
    gStyle->SetStatY(0.9);

    bool doRatio = 0;
    TPaveText lumi=add_lumi(year,doRatio);
    TPaveText cms=add_CMS(doRatio);
    TPaveText pre=add_Preliminary(channel, doRatio);

    meangraph->SetName("mean");
    meangraph->SetTitle("");
    meangraph->GetXaxis()->SetTitle("Mass");
    meangraph->GetYaxis()->SetTitle("Mean Fit Parameter");
    meangraph->SetMarkerStyle(8);
    meangraph->Draw("AP");
    meangraph->Fit(meanfit);
    meanfit->Draw("same");
    lumi.Draw();
    cms.Draw();
    pre.Draw();


    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    normgraph->SetName("norm");
    normgraph->SetTitle("");
    normgraph->GetXaxis()->SetTitle("Mass");
    normgraph->GetYaxis()->SetTitle("Norm Fit Parameter");
    normgraph->SetMarkerStyle(8);
    normgraph->Draw("AP");
    normgraph->Fit(normfit);
    normfit->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    sigmagraph->SetName("sigma");
    sigmagraph->SetTitle("");
    sigmagraph->GetXaxis()->SetTitle("Mass");
    sigmagraph->GetYaxis()->SetTitle("Sigma Fit Parameter");
    sigmagraph->SetMarkerStyle(8);
    sigmagraph->Draw("AP");
    sigmagraph->Fit(sigmafit);
    sigmafit->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();
    cout<<"the year is "<<year<<endl;
    
    //doing the interpolation to a spline-like function to save to RooWorkspace 
    //map<string, Roo*> signaltemplates;
    //map<string, Roo*> signalnorms;
    //map<string, Roo*> signaldatasets;
    //map<string, Roo*> x;
    //map<string, Roo*> m;
    //map<string, Roo*> s;

    RooRealVar * MH    = new RooRealVar("MH","MH", 18, 63);
    RooRealVar * Mll   = new RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0);
    //RooRealVar * MHerr = new RooRealVar("MHerr","MHerr", 0, -4 , 4);
    cout<<" mean formula ("+to_string(meanfit->GetParameter(0))+" + "+to_string(meanfit->GetParameter(1))+"*@0)"<<endl;

    RooFormulaVar * intMean = new RooFormulaVar("intMean","intMean",("("+to_string(meanfit->GetParameter(0))+" + "+to_string(meanfit->GetParameter(1))+"*@0)").c_str(),RooArgList(*MH));

    RooFormulaVar * intSigma = new RooFormulaVar("intSigma","intSigma",("("+to_string(sigmafit->GetParameter(0))+" + "+to_string(sigmafit->GetParameter(1))+"*@0+"+to_string(sigmafit->GetParameter(2))+"*@0*@0+"+to_string(sigmafit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));

    //RooFormulaVar * intNorm = new RooFormulaVar("signal_norm","signal_norm",("("+to_string(normfit->GetParameter(0))+" + "+to_string(normfit->GetParameter(1))+"*@0+"+to_string(normfit->GetParameter(2))+"*@0*@0+"+to_string(normfit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));

    RooGaussian * intSignalTemplate = new  RooGaussian(("signal_"+systematic_name+"_"+output).c_str(),   ("signal_"+systematic_name+"_"+output).c_str(),*Mll, *intMean, *intSigma );
    cout<<"made signal template"<<endl;

    //for mass in range(16,66):;
    //    massEval = meanfit.Eval(mass);
    //    normEval = normfit.Eval(mass);
    //    sigmaEval = sigmafit.Eval(mass);
    //    signalnorms[str(mass)] = normEval;
    //    print("evaluation of mean at ",mass," is ",massEval, " generating signal template ");
    //    print("evaluation of norm at ",mass," is ",normEval);
    //    print("evaluation of sigma at ",mass," is ",sigmaEval);
    //    x[str(mass)] = RooRealVar("mll",    "mll",massEval-2.0,massEval+2.0);
    //    m[str(mass)] = RooRealVar("mean",    "mean",massEval,"GeV");
    //    s[str(mass)] = RooRealVar("sigma",    "sigma", sigmaEval,"GeV");

    wsp->import(*intSignalTemplate);
    }
    if(type.compare("voigtian")==0){
    cout<<"interpolating!"<<endl;
    TF1 * meanfitup = new TF1("meanfitup","pol1",16,66);
    TF1 * meanfit = new TF1("meanfit","pol1",16,66);
    TF1 * meanfitdown = new TF1("meanfitdown","pol1",16,66);
    TF1 * normfitup = new TF1("normfitup","pol3",16,66);
    TF1 * normfit = new TF1("normfit","pol3",16,66);
    TF1 * normfitdown = new TF1("normfitdown","pol3",16,66);
    TF1 * alphafitup = new TF1("alphafitup","pol3",16,66);
    TF1 * alphafit = new TF1("alphafit","pol3",16,66);
    TF1 * alphafitdown = new TF1("alphafitdown","pol3",16,66);
    TF1 * sigmafitup = new TF1("sigmafitup","pol3",16,66);
    TF1 * sigmafit = new TF1("sigmafit","pol3",16,66);
    TF1 * sigmafitdown = new TF1("sigmafitdown","pol3",16,66);

    TGraphErrors * meangraphup = new TGraphErrors();
    TGraphErrors * meangraph = new TGraphErrors();
    TGraphErrors * meangraphdown = new TGraphErrors();
    TGraphErrors * normgraphup = new TGraphErrors();
    TGraphErrors * normgraph = new TGraphErrors();
    TGraphErrors * normgraphdown = new TGraphErrors();
    TGraphErrors * sigmagraphup = new TGraphErrors();
    TGraphErrors * sigmagraph = new TGraphErrors();
    TGraphErrors * sigmagraphdown = new TGraphErrors();
    TGraphErrors * alphagraphup = new TGraphErrors();
    TGraphErrors * alphagraph = new TGraphErrors();
    TGraphErrors * alphagraphdown = new TGraphErrors();

    int c = 0;
    //LimitShape * sp = (LimitShape*) inputshapes["a20"];
    cout<<"shape exist? name: "<<inputshapes["a20"]<<endl;
    cout<<"shape exist? name: "<<inputshapes["a20"]->shape_name<<endl;
    for(auto const& x: inputshapes){
        string mass = x.first.substr(1,2); 
        LimitShape * sp = x.second;

        meangraphup->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(0))->getVal()*1.05
        );
        meangraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(0))->getVal()
        );
        meangraphdown->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(0))->getVal()/1.05
        );
        alphagraphup->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(1))->getVal()*1.1
        );
        alphagraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(1))->getVal()
        );
        alphagraphdown->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(1))->getVal()/1.1
        );
        sigmagraphup->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(2))->getVal()*1.2
        );
        sigmagraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(2))->getVal()
        );
        sigmagraphdown->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["voigtian"+mass]->At(2))->getVal()/1.2
        );
        normgraphup->SetPoint(c,
            stof((mass)),
            sp->data->sumEntries()*1.10
        );
        normgraph->SetPoint(c,
            stof((mass)),
            sp->data->sumEntries()
        );
        normgraphdown->SetPoint(c,
            stof((mass)),
            sp->data->sumEntries()/1.10
        );

        c++;
    }
    TCanvas * can = new TCanvas("splinefits","splinefits",600,600);
    can->cd();
    gStyle->SetOptStat(1); 
    gStyle->SetOptFit(1); 
    gStyle->SetStatX(0.6);
    gStyle->SetStatY(0.9);

    bool doRatio = 0;
    TPaveText lumi=add_lumi(year,doRatio);
    TPaveText cms=add_CMS(doRatio);
    TPaveText pre=add_Preliminary(channel, doRatio);

    meanfitup->SetLineColor(kBlue);
    meanfitdown->SetLineColor(kBlue);

    meangraph->SetName("mean");
    meangraph->SetTitle("");
    meangraph->GetXaxis()->SetTitle("Mass");
    meangraph->GetYaxis()->SetTitle("Mean Fit Parameter");
    meangraph->SetMarkerStyle(8);
    meangraph->Draw("AP");
    meangraphup->Fit(meanfitup);
    meangraph->Fit(meanfit);
    meangraphdown->Fit(meanfitdown);
    meanfitup->Draw("same");
    meanfit->Draw("same");
    meanfitdown->Draw("same");
    lumi.Draw();
    cms.Draw();
    pre.Draw();


    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    normfitup->SetLineColor(kBlue);
    normfitdown->SetLineColor(kBlue);

    normgraph->SetName("norm");
    normgraph->SetTitle("");
    normgraph->GetXaxis()->SetTitle("Mass");
    normgraph->GetYaxis()->SetTitle("Norm Fit Parameter");
    normgraph->SetMarkerStyle(8);
    normgraph->Draw("AP");
    normgraphup->Fit(normfitup);
    normgraph->Fit(normfit);
    normgraphdown->Fit(normfitdown);
    normfitup->Draw("same");
    normfit->Draw("same");
    normfitdown->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    alphafitup->SetLineColor(kBlue);
    alphafitdown->SetLineColor(kBlue);

    alphagraph->SetName("alpha");
    alphagraph->SetTitle("");
    alphagraph->GetXaxis()->SetTitle("Mass");
    alphagraph->GetYaxis()->SetTitle("Alpha Fit Parameter");
    alphagraph->SetMarkerStyle(8);
    alphagraph->Draw("AP");
    alphagraphup->Fit(alphafitup);
    alphagraph->Fit(alphafit);
    alphagraphdown->Fit(alphafitdown);
    alphafitup->Draw("same");
    alphafit->Draw("same");
    alphafitdown->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_AlphaConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_AlphaConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    sigmafitup->SetLineColor(kBlue);
    sigmafitdown->SetLineColor(kBlue);

    sigmagraph->SetName("sigma");
    sigmagraph->SetTitle("");
    sigmagraph->GetXaxis()->SetTitle("Mass");
    sigmagraph->GetYaxis()->SetTitle("Sigma Fit Parameter");
    sigmagraph->SetMarkerStyle(8);
    sigmagraph->Draw("AP");
    sigmagraphup->Fit(sigmafitup);
    sigmagraph->Fit(sigmafit);
    sigmagraphdown->Fit(sigmafitdown);
    sigmafitup->Draw("same");
    sigmafit->Draw("same");
    sigmafitdown->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();
    cout<<"the year is "<<year<<endl;
    
    //doing the interpolation to a spline-like function to save to RooWorkspace 
    //map<string, Roo*> signaltemplates;
    //map<string, Roo*> signalnorms;
    //map<string, Roo*> signaldatasets;
    //map<string, Roo*> x;
    //map<string, Roo*> m;
    //map<string, Roo*> s;

    RooRealVar * MH    = new RooRealVar("MH","MH", 18, 63);
    RooRealVar * Mll   = new RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0);
    //RooRealVar * MHerr = new RooRealVar("MHerr","MHerr", 0, -4 , 4);
    cout<<" mean formula ("+to_string(meanfit->GetParameter(0))+" + "+to_string(meanfit->GetParameter(1))+"*@0)"<<endl;

    RooFormulaVar * intMean = new RooFormulaVar(("intMean_"+systematic_name+"_"+output).c_str(),("intMean_"+systematic_name+"_"+output).c_str(),("("+to_string(meanfit->GetParameter(0))+" + "+to_string(meanfit->GetParameter(1))+"*@0)").c_str(),RooArgList(*MH));

    RooFormulaVar * intSigma = new RooFormulaVar(("intSigma_"+systematic_name+"_"+output).c_str(),("intSigma_"+systematic_name+"_"+output).c_str(),("("+to_string(sigmafit->GetParameter(0))+" + "+to_string(sigmafit->GetParameter(1))+"*@0+"+to_string(sigmafit->GetParameter(2))+"*@0*@0+"+to_string(sigmafit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));

    RooFormulaVar * intAlpha = new RooFormulaVar(("intAlpha_"+systematic_name+"_"+output).c_str(),("intAlpha_"+systematic_name+"_"+output).c_str(),("("+to_string(alphafit->GetParameter(0))+" + "+to_string(alphafit->GetParameter(1))+"*@0+"+to_string(alphafit->GetParameter(2))+"*@0*@0+"+to_string(alphafit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));

    RooFormulaVar * intNorm = new RooFormulaVar(("signal_"+systematic_name+"_"+output+"_norm").c_str(),("signal_"+systematic_name+"_"+output+"_norm").c_str(),("("+to_string(normfit->GetParameter(0))+" + "+to_string(normfit->GetParameter(1))+"*@0+"+to_string(normfit->GetParameter(2))+"*@0*@0+"+to_string(normfit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));

    RooVoigtian * intSignalTemplate = new  RooVoigtian(("signal_"+systematic_name+"_"+output).c_str(),   ("signal_"+systematic_name+"_"+output).c_str(),*Mll, *intMean, *intAlpha, *intSigma );
    cout<<"made signal template"<<endl;

    //for mass in range(16,66):;
    //    massEval = meanfit.Eval(mass);
    //    normEval = normfit.Eval(mass);
    //    sigmaEval = sigmafit.Eval(mass);
    //    signalnorms[str(mass)] = normEval;
    //    print("evaluation of mean at ",mass," is ",massEval, " generating signal template ");
    //    print("evaluation of norm at ",mass," is ",normEval);
    //    print("evaluation of sigma at ",mass," is ",sigmaEval);
    //    x[str(mass)] = RooRealVar("mll",    "mll",massEval-2.0,massEval+2.0);
    //    m[str(mass)] = RooRealVar("mean",    "mean",massEval,"GeV");
    //    s[str(mass)] = RooRealVar("sigma",    "sigma", sigmaEval,"GeV");

    wsp->import(*intSignalTemplate);//importing signal template
    wsp->import(*intNorm); //importing normalization


    }
    if(type.compare("doublegaussian")==0){
    cout<<"interpolating!"<<endl;
    TF1 * meanfit = new TF1("meanfit","pol1",16,66);
    TF1 * normfit = new TF1("normfit","pol3",16,66);
    TF1 * sigma1fit = new TF1("sigma1fit","pol3",16,66);
    TF1 * sigma2fit = new TF1("sigma2fit","pol3",16,66);
    TGraphErrors * meangraph = new TGraphErrors();
    TGraphErrors * normgraph = new TGraphErrors();
    TGraphErrors * sigma2graph = new TGraphErrors();
    TGraphErrors * sigma1graph = new TGraphErrors();
    int c = 0;
    //LimitShape * sp = (LimitShape*) inputshapes["a20"];
    cout<<"shape exist? name: "<<inputshapes["a20"]<<endl;
    cout<<"shape exist? name: "<<inputshapes["a20"]->shape_name<<endl;
    for(auto const& x: inputshapes){
        string mass = x.first.substr(1,2); 
        LimitShape * sp = x.second;

        meangraph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["doublegaussian"+mass]->At(0))->getVal()
        );
        sigma1graph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["doublegaussian"+mass]->At(1))->getVal()
        );
        sigma2graph->SetPoint(c,
            stof((mass)),
            ((RooRealVar*)sp->coeffs["doublegaussian"+mass]->At(2))->getVal()
        );
        normgraph->SetPoint(c,
            stof((mass)),
            sp->data->sumEntries()
        );

        c++;
    }
    TCanvas * can = new TCanvas("splinefits","splinefits",600,600);
    can->cd();
    gStyle->SetOptStat(1); 
    gStyle->SetOptFit(1); 
    //gStyle->SetStatX(0.6);
    //gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.5);
    TGaxis::SetMaxDigits(2);

    bool doRatio = 0;
    TPaveText lumi=add_lumi(year,doRatio);
    TPaveText cms=add_CMS(doRatio);
    TPaveText pre=add_Preliminary(channel, doRatio);

    meangraph->SetName("mean");
    meangraph->SetTitle("");
    meangraph->GetXaxis()->SetTitle("Mass");
    meangraph->GetYaxis()->SetTitle("Mean Fit Parameter");
    meangraph->SetMarkerStyle(8);
    meangraph->Draw("AP");
    meangraph->Fit(meanfit);
    meanfit->Draw("same");
    lumi.Draw();
    cms.Draw();
    pre.Draw();


    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    normgraph->SetName("norm");
    normgraph->SetTitle("");
    normgraph->GetXaxis()->SetTitle("Mass");
    normgraph->GetYaxis()->SetTitle("Norm Fit Parameter");
    normgraph->SetMarkerStyle(8);
    normgraph->Draw("AP");
    normgraph->Fit(normfit);
    normfit->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    sigma1graph->SetName("sigma1");
    sigma1graph->SetTitle("");
    sigma1graph->GetXaxis()->SetTitle("Mass");
    sigma1graph->GetYaxis()->SetTitle("Inner Sigma Fit Parameter");
    sigma1graph->SetMarkerStyle(8);
    sigma1graph->Draw("AP");
    sigma1graph->Fit(sigma1fit);
    sigma1fit->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_Sigma1Constraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_Sigma1Constraint_"+systematic_name+"_"+output+".png"));
    can->Clear();

    sigma2graph->SetName("sigma2");
    sigma2graph->SetTitle("");
    sigma2graph->GetXaxis()->SetTitle("Mass");
    sigma2graph->GetYaxis()->SetTitle("Outer Sigma Fit Parameter");
    sigma2graph->SetMarkerStyle(8);
    sigma2graph->Draw("AP");
    sigma2graph->Fit(sigma2fit);
    sigma2fit->Draw("same");

    lumi.Draw();
    cms.Draw();
    pre.Draw();

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_Sigma2Constraint_"+systematic_name+"_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_Sigma2Constraint_"+systematic_name+"_"+output+".png"));
    can->Clear();
    cout<<"the year is "<<year<<endl;
    
    //doing the interpolation to a spline-like function to save to RooWorkspace 
    //map<string, Roo*> signaltemplates;
    //map<string, Roo*> signalnorms;
    //map<string, Roo*> signaldatasets;
    //map<string, Roo*> x;
    //map<string, Roo*> m;
    //map<string, Roo*> s;

    RooRealVar * MH    = new RooRealVar("MH","MH", 18, 63);
    RooRealVar * Mll   = new RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0);
    //RooRealVar * MHerr = new RooRealVar("MHerr","MHerr", 0, -4 , 4);
    cout<<" mean formula ("+to_string(meanfit->GetParameter(0))+" + "+to_string(meanfit->GetParameter(1))+"*@0)"<<endl;

    RooFormulaVar * intMean = new RooFormulaVar("intMean","intMean",("("+to_string(meanfit->GetParameter(0))+" + "+to_string(meanfit->GetParameter(1))+"*@0)").c_str(),RooArgList(*MH));

    RooFormulaVar * intSigma2 = new RooFormulaVar("intSigma2","intSigma2",("("+to_string(sigma2fit->GetParameter(0))+" + "+to_string(sigma2fit->GetParameter(1))+"*@0+"+to_string(sigma2fit->GetParameter(2))+"*@0*@0+"+to_string(sigma2fit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));

    RooFormulaVar * intSigma1 = new RooFormulaVar("intSigma1","intSigma1",("("+to_string(sigma1fit->GetParameter(0))+" + "+to_string(sigma1fit->GetParameter(1))+"*@0+"+to_string(sigma1fit->GetParameter(2))+"*@0*@0+"+to_string(sigma1fit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));

    //RooFormulaVar * intNorm = new RooFormulaVar("signal_norm","signal_norm",("("+to_string(normfit->GetParameter(0))+" + "+to_string(normfit->GetParameter(1))+"*@0+"+to_string(normfit->GetParameter(2))+"*@0*@0+"+to_string(normfit->GetParameter(3))+"*@0*@0*@0)").c_str(),RooArgList(*MH));
    RooGaussModel * temppdf_1 = new RooGaussModel(
                                        (TString)("sigfit_1"),
                                        (TString)("sigfit_1"),
                                        *Mll,
                                        *intMean,
                                        *intSigma1
                                        );
    RooGaussModel * temppdf_2 = new RooGaussModel(
                                        (TString)("sigfit_2"),
                                        (TString)("sigfit_2"),
                                        *Mll,
                                        *intMean,
                                        *intSigma2
                                        );
    RooRealVar * frac_pdfs = new RooRealVar("frac_pdfs","fraction between gauss 1 and 2",0.5);

    RooAddModel * intSignalTemplate = new  RooAddModel(
                ("signal_"+channel).c_str(),   
                ("signal_"+channel).c_str(),
                RooArgList(*temppdf_1,*temppdf_2),
                *frac_pdfs
                );
    cout<<"made signal template"<<endl;

    //for mass in range(16,66):;
    //    massEval = meanfit.Eval(mass);
    //    normEval = normfit.Eval(mass);
    //    sigmaEval = sigmafit.Eval(mass);
    //    signalnorms[str(mass)] = normEval;
    //    print("evaluation of mean at ",mass," is ",massEval, " generating signal template ");
    //    print("evaluation of norm at ",mass," is ",normEval);
    //    print("evaluation of sigma at ",mass," is ",sigmaEval);
    //    x[str(mass)] = RooRealVar("mll",    "mll",massEval-2.0,massEval+2.0);
    //    m[str(mass)] = RooRealVar("mean",    "mean",massEval,"GeV");
    //    s[str(mass)] = RooRealVar("sigma",    "sigma", sigmaEval,"GeV");

    wsp->import(*intSignalTemplate);


    }

    return;
}

void setFunctionName(RooAbsPdf * function,string name){
    function->SetName(name.c_str());
    
}

void Office::createTxtfile(){
    txtfile.open(outputdir+"/"+"datacard_full_"+output+".txt");
    return;
}
void Office::printDatacard(bool isrealdata){
    txtfile << ("# Sams Datacard \n");
    txtfile << ("#contains real data "+to_string(isrealdata)+"\n"); //number of bins - only one category ... no control region
    txtfile << ("imax 1\n"); //number of bins - only one category ... no control region
    txtfile << ("jmax 2\n"); //number of processes minus 1
    txtfile << ("kmax *\n"); //number of nuisance parameters
    txtfile << ("---------------\n");
    txtfile << ("shapes * bin1 HToAAWorkspace_full_"+output+".root w:$PROCESS\n");
    txtfile << ("---------------\n");

    txtfile << ("bin         bin1   \n");
    txtfile << ("observation   -1 \n"); // for parametric fit this needs to be -1

    txtfile << ("------------------------------\n");
    txtfile << ("bin                    ");
    txtfile << (" bin1 ");
    txtfile << ("bin1 bin1 \n");
    txtfile << ("process                ");
    txtfile << (" signal_"+output+"               irBkg_"+output+"      Bkg_"+output+"\n");
    txtfile << ("process                ");
    txtfile << ("0 1 2");
    txtfile << ("\n");
    txtfile << ("rate                   ");
    txtfile << ("1 1 1 \n");
    txtfile << ("------------------------------\n");
    txtfile << ("lumi     lnN              1.016    1.016    1.016\n");



    vector<string> signal_names = {"a15","a20","a25","a30","a35","a40","a45","a50","a55","a60"};


    for (auto const& sys: systematics){
        cout << "finding pull for systematic " <<sys<<endl;
        cout << "working on shape "<<(shapes["Nominal_a40"]->shape_name).c_str()<<endl;
        cout << "data entries "<<to_string((shapes["Nominal_a40"]->data)->sumEntries()).c_str()<<endl;
        cout << "data entries Up"<<to_string((shapes[sys+"Up_a40"]->data)->sumEntries()).c_str()<<endl;

        float signalpullup = findPull(
            (shapes["Nominal_a40"]->data),
            (shapes[sys+"Up_a40"]->data),
            (shapes[sys+"Down_a40"]->data));
        float ZZpullup = findPull(
            (shapes["Nominal_irBkg"]->data),
            (shapes[sys+"Up_irBkg"]->data),
            (shapes[sys+"Down_irBkg"]->data));
        
        txtfile << (sys+output+"   lnN             "+to_string(round_up(signalpullup,5)) +"       "+to_string(round_up(ZZpullup,5)) +"         -    \n");


    }

    for(auto const& x: shapes["Nominal_irBkg"]->coeffs){
        for(const auto&& obj: *(shapes["Nominal_irBkg"]->coeffs)[x.first]){
            RooRealVar * var = (RooRealVar*) obj;
            txtfile << var->GetName()<<"  param "+to_string(round_up(var->getVal(),5))+"   "+to_string(round_up(var->getError(),5)) << endl;
        }
    }
    for(auto const& x: shapes["Nominal_Bkg"]->coeffs){
        for(const auto&& obj: *(shapes["Nominal_Bkg"]->coeffs)[x.first]){
            RooRealVar * var = (RooRealVar*) obj;
            txtfile << var->GetName()<<"  param "+to_string(round_up(var->getVal(),5))+"   "+to_string(round_up(var->getError(),5)) << endl;
        }
    }
    //for(auto const& x: shapes["Nominal_a40"]->coeffs){
    //    for(const auto&& obj: *(shapes["Nominal_a40"]->coeffs)[x.first]){
    //        RooRealVar * var = (RooRealVar*) obj;
    //        txtfile << var->GetName()<<"  param "+to_string(var->getVal())+"   "+to_string(var->getError()) << endl;
    //    }
    //}
    string othername;
    string firstname;
    float error;

    error = 0.10;
    txtfile <<"signal_"+output+"_norm"<<"  param "+to_string(round_up(shapes["Nominal_a40"]->data->sumEntries(),5))+"   "+ to_string(round_up(shapes["Nominal_a40"]->data->sumEntries()*error,5)) << endl;
    for(auto const& x: shapes["Nominal_a40"]->coeffs){
        for(const auto&& obj: *(shapes["Nominal_a40"]->coeffs)[x.first]){
            RooRealVar * var = (RooRealVar*) obj;
            firstname = (string) (var->GetName());
            //if((firstname).find("mean")!=string::npos){othername=("intMean_"+output).c_str();error=0.001;} 
            if((firstname).find("mean")!=string::npos){othername=("intMean_"+output).c_str();error=0.05;} 
            if((firstname).find("sigma")!=string::npos){othername=("intSigma_"+output).c_str();error=0.20;} 
            if((firstname).find("alpha")!=string::npos){othername=("intAlpha_"+output).c_str();error=0.10;} 
            //txtfile << othername<<"  param "+to_string(var->getVal())+"   "+to_string(var->getError()) << endl;
            txtfile << othername<<"  param "+to_string(round_up(var->getVal(),5))+"   "+ to_string(round_up(var->getVal()*error,5))<< endl;
        }
    }

    //txtfile <<" intMean param  1.0 1.05"<< endl;
    //txtfile <<" intSigma param  1.0 1.2"<< endl;
    //txtfile <<" intAlpha param  1.0 1.2"<< endl;

    /*
    for systematic in self.systematics:
        signalpullup   = findPull(self.shapes["Nominal"]["a40"].data,self.shapes[systematic+"Up"]["a40"].data,self.shapes[systematic+"Down"]["a40"].data)
        ZZpullup   = findPull(self.shapes["Nominal"]["irBkg"].data,self.shapes[systematic+"Up"]["irBkg"].data,self.shapes[systematic+"Down"]["irBkg"].data)
        self.txtfile.write(systematic+args.output+"   lnN              {0:.9f}        {1:.9f}          -    \n".format(signalpullup,ZZpullup))


    for coeffkey, coefflist in self.shapes["Nominal"]["irBkg"].coeffs.items():
        for coeff in coefflist:
            #make 0th value the norm?? and 1/sqrt(N) events for the error?? idk....
            self.txtfile.write(coeff.GetName()+" param "+str(coeff.getVal())+" "+str(coeff.getError())+"\n")

    for coeffkey, coefflist in self.shapes["Nominal"]["Bkg"].coeffs.items():
        for coeff in coefflist:
            self.txtfile.write(coeff.GetName()+" param "+str(coeff.getVal())+" "+str(coeff.getError())+"\n")

    */
    return;
}
void Office::createTxtfilePerMass(){
    string mass="a";
    for(int i=18;i<63;i++){ 
    mass="a"+to_string(i);
    txtfiles[mass].open(outputdir+"/"+"datacard_"+mass+"_"+output+".txt");
    }
    return;
}
void Office::printDatacardPerMass(bool isrealdata){
    for(auto const& x: txtfiles){
        txtfiles[x.first] << ("# Sams Datacard \n");
        txtfiles[x.first] << ("#contains real data "+to_string(isrealdata)+"\n"); //number of bins - only one category ... no control region
        txtfiles[x.first] << ("imax 1\n"); //number of bins - only one category ... no control region
        txtfiles[x.first] << ("jmax 2\n"); //number of processes minus 1
        txtfiles[x.first] << ("kmax *\n"); //number of nuisance parameters
        txtfiles[x.first] << ("---------------\n");
        txtfiles[x.first] << ("shapes * bin1 HToAAWorkspace_full_"+output+".root w:$PROCESS\n");
        txtfiles[x.first] << ("---------------\n");

        txtfiles[x.first] << ("bin         bin1   \n");
        txtfiles[x.first] << ("observation   -1 \n"); // for parametric fit this needs to be -1

        txtfiles[x.first] << ("------------------------------\n");
        txtfiles[x.first] << ("bin                    ");
        txtfiles[x.first] << (" bin1 ");
        txtfiles[x.first] << ("bin1 bin1 \n");
        txtfiles[x.first] << ("process                ");
        txtfiles[x.first] << (" signal_"+output+"               irBkg_"+output+"      Bkg_"+output+"\n");
        txtfiles[x.first] << ("process                ");
        txtfiles[x.first] << ("0 1 2");
        txtfiles[x.first] << ("\n");
        txtfiles[x.first] << ("rate                   ");
        txtfiles[x.first] << ("1 1 1 \n");
        txtfiles[x.first] << ("------------------------------\n");
        txtfiles[x.first] << ("lumi     lnN              1.016    1.016    1.016\n");



        vector<string> signal_names = {"a15","a20","a25","a30","a35","a40","a45","a50","a55","a60"};
        vector<int> signal_masses = {20,25,30,35,40,45,50,55,60};

        int closestMass = 0;
        int mass = stoi(((x.first).substr(1,2)));
        if (mass%5 == 0 ){closestMass=mass;}
        for(int i=0; i<8;i++){
            if(mass<20){closestMass=20;break;}
            if(mass>signal_masses[i] && mass<signal_masses[i+1]){closestMass=signal_masses[i];break;}
            if(mass>60){closestMass=60; break;}
             

        }
        
        for (auto const& sys: systematics){

            float signalpullup = findPull(
                (shapes["Nominal_a"+to_string(closestMass)]->data),
                (shapes[sys+"Up_a"+to_string(closestMass)]->data),
                (shapes[sys+"Down_a"+to_string(closestMass)]->data));
            float ZZpullup = findPull(
                (shapes["Nominal_irBkg"]->data),
                (shapes[sys+"Up_irBkg"]->data),
                (shapes[sys+"Down_irBkg"]->data));
            
            txtfiles[x.first] << (sys+output+"   lnN             "+to_string(round_up(signalpullup,5)) +"       "+to_string(round_up(ZZpullup,5)) +"         -    \n");


        }

        for(auto const& cx: shapes["Nominal_irBkg"]->coeffs){
            for(const auto&& obj: *(shapes["Nominal_irBkg"]->coeffs)[cx.first]){
                RooRealVar * var = (RooRealVar*) obj;
                txtfiles[x.first] << var->GetName()<<"  param "+to_string(round_up(var->getVal(),5))+"   "+to_string(round_up(var->getError(),5)) << endl;
            }
        }
        for(auto const& cx: shapes["Nominal_Bkg"]->coeffs){
            for(const auto&& obj: *(shapes["Nominal_Bkg"]->coeffs)[cx.first]){
                RooRealVar * var = (RooRealVar*) obj;
                txtfiles[x.first] << var->GetName()<<"  param "+to_string(round_up(var->getVal(),5))+"   "+to_string(round_up(var->getError(),5)) << endl;
            }
        }

        wsp->var("MH")->setVal(mass);
        txtfiles[x.first] << ("signal_"+output+"_norm")+"   param " << round_up(wsp->function(("signal_"+output+"_norm").c_str())->getVal(),5)<<"   "<<round_up(wsp->function(("signal_"+output+"_norm").c_str())->getVal()*0.10,5) << endl;
        //txtfiles[x.first] << ("intMean_"+output)+"   param " << wsp->function(("intMean_"+output).c_str())->getVal()<<"   "<<(wsp->function(("intMean_"+output).c_str())->getVal()*0.001) << endl;
        txtfiles[x.first] << ("intMean_"+output)+"   param " << round_up(wsp->function(("intMean_"+output).c_str())->getVal(),5)<<"   "<<round_up(wsp->function(("intMean_"+output).c_str())->getVal()*0.05,5) << endl;
        txtfiles[x.first] << ("intSigma_"+output)+"  param " << round_up(wsp->function(("intSigma_"+output).c_str())->getVal(),5)<<"   "<<round_up(wsp->function(("intSigma_"+output).c_str())->getVal()*0.20,5) << endl;
        txtfiles[x.first] << ("intAlpha_"+output)+"  param " << round_up(wsp->function(("intAlpha_"+output).c_str())->getVal(),5)<<"   "<<round_up(wsp->function(("intAlpha_"+output).c_str())->getVal()*0.20,5) << endl;

        } //mass loop
    return;
}

#endif
