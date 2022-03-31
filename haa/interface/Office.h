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

#include "cppLimits/haa/interface/LimitShape.h"

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
    void interpolateParameters(map<string,LimitShape*>& inputshapes);
    void createTxtfile();
    void printDatacard(bool isrealdata);
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
void Office::interpolateParameters(map<string,LimitShape*>& inputshapes){
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


    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_MeanConstraint_"+output+".png"));
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

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_NormConstraint_"+output+".png"));
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

    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+output+".pdf"));
    can->SaveAs((TString)(outputdir+"/"+"DiMuonMass_SigmaConstraint_"+output+".png"));
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

    RooGaussian * intSignalTemplate = new  RooGaussian(("signal_"+channel).c_str(),   ("signal_"+channel).c_str(),*Mll, *intMean, *intSigma );
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


    return;
}

void Office::createTxtfile(){
    txtfile.open(outputdir+"/"+"datacard_full_"+output+".txt");
    return;
}
void setFunctionName(RooAbsPdf * function,string name){
    function->SetName(name.c_str());
    
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
    txtfile << (" signal_"+channel+"               irBkg_"+channel+"      Bkg_"+channel+"\n");
    txtfile << ("process                ");
    txtfile << ("0 1 2");
    txtfile << ("\n");
    txtfile << ("rate                   ");
    txtfile << ("1 1 1 \n");
    txtfile << ("------------------------------\n");
    txtfile << ("lumi     lnN              1.01    1.01    1.01\n");



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
        
        txtfile << (sys+output+"   lnN             "+to_string(signalpullup) +"       "+to_string(ZZpullup) +"         -    \n");


    }

    for(auto const& x: shapes["Nominal_irBkg"]->coeffs){
        for(const auto&& obj: *(shapes["Nominal_irBkg"]->coeffs)[x.first]){
            RooRealVar * var = (RooRealVar*) obj;
            txtfile << var->GetName()<<"  param "+to_string(var->getVal())+"   "+to_string(var->getError()) << endl;
        }
    }
    for(auto const& x: shapes["Nominal_Bkg"]->coeffs){
        for(const auto&& obj: *(shapes["Nominal_Bkg"]->coeffs)[x.first]){
            RooRealVar * var = (RooRealVar*) obj;
            txtfile << var->GetName()<<"  param "+to_string(var->getVal())+"   "+to_string(var->getError()) << endl;
        }
    }
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

#endif
