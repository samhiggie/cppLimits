#ifndef LimitShape_h
#define LimitShape_h

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
#include "TObject.h"
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
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "RooVoigtian.h"
#include "RooFitResult.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"
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


namespace po = boost::program_options;
using namespace std;
using namespace ROOT;
using namespace RooFit;

#include <iostream>
#include <utility>
#include <vector>
#include <cassert>

#include "plot.h"

class LimitShape{
    public:
        po::variables_map shape_vm;
        string           output="default";//string for naming
        string           outputdir="default_folder";
        string           channel = "mmmt";
        //po::positional_options_description *shape_p;
        //po::options_description *shape_desc;
        Int_t            treeNum;
        TTree            *ttree;
        TTree            *newtree;//for normalization rescaling
        string           shape_name;
        string           shape_dist;
        string           type="shapetype";//RooFit Pdf class
        string           systematic="Nominal";//systematic for shapes

        TFile           *tfile;
        RooRealVar *poi;
        RooRealVar *norm;
        RooRealVar *eventweight;
        RooRealVar *AMass;
        RooRealVar *m_vis;
        map<string,RooAbsPdf*> pdf;
        map<string,float> nll;
        map<string,float> deltanll; // dictionary for the change in nlls from the pdf fits
        map<string,TList*> coeffs; // dictionary for RooRealVars
        map<string,TList*> finalcoeffs; // dictionary for RooRealVars
        RooDataSet *data;
        map<string,RooFitResult*> fitresults; // dictionary to hold fit results from the PDFs;
        map<string,RooFitResult*> finalfitresults;// dictionary to hold fit results from the PDFs;
        map<string,RooPlot*> plotframe;
        map<string,TCanvas*> canvas;
        float rangeMin  = 0.0;
        float rangeMax  = 0.0;
        //float coeffMin  = -10000.0;
        //float coeffMax  =  10000.0;
        vector<float> coeffMin;
        vector<float> coeffMax;
        vector<float> coeffvect;
        int maxscans  =  20;
        float shape_binning  = 16;
        //float shape_binning  = 4;
        int recursion = 0;
        int year = 2017;

        //current variables - make this into a class later
        //datatype
        unsigned int     value;

        //branch
        TBranch     b_value;

        //new variables
        float  newvar;

        //new branch
        TBranch b_newvar;

        //////////////////////////////////////
        // function methods and inheritance
        //////////////////////////////////////

        LimitShape(  po::variables_map vm,
                //po::options_description *desc,
                //po::positional_options_description *p,
                string dist);
        ~LimitShape();

        void Init(TTree *ttree);
        void connectFile(string path);
        void setsystematic(string sys);
        int  fillTree(string pathtotree);
        void rescaleFinalweight(float num);
        void fillPDFs(string type,int order,string dist);
        void fillDataset(string cutstring);
        void fitToData();
        void finalFitToData(string type,float sfr);
        void recursiveFitToData(string type);
        void scanFitToData(string type);
        void printParamsNLLs();
        void createPlots(string type,int order,string dist);
        //other containers
};

//constructors
LimitShape::LimitShape(po::variables_map vm,
    //po::options_description *desc,
    //po::positional_options_description *p,
    string dist){
    //setting the program options 

    shape_vm = vm;
    //shape_p = p;
    //shape_desc = desc;

    //arglist = vm;
    shape_dist = dist;
    cout<<"creating shape with distribution "<<dist<<endl;

    output = shape_vm["output-file"].as<std::string>();
    outputdir = shape_vm["output-dir"].as<std::string>();
    year = stoi(output.substr(0,4));
    channel = shape_vm["channel"].as<std::string>();
    cout<<"creating directory "<<outputdir<<endl;
    //output = arglist["output-file"].as<std::string>();
    //creating a directory
    const int dir_err = mkdir(outputdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        printf("Error creating directory!\n");
        //exit(1);
    }
  
  
    //system("dir");

    //setting the vars
    poi = new RooRealVar("mll","m_{#mu#mu}", 18, 62);
    norm = new RooRealVar("norm",   "Normalization",1.0,0.0,1000000.0);
    eventweight = new RooRealVar("weight",   "weight",0.0,5.0);

}
LimitShape::~LimitShape(){}

// Functions
void LimitShape::connectFile(string path)
{

    std::cout<<"trying to open "<<path.c_str()<<std::endl;
    //unique::ptr<TFile> tfile(TFile::Open(path.c_str()));
    //tfile->TFile::Open(path.c_str(),"read");
    tfile = new TFile(path.c_str(),"READ");
    tfile->cd();
    //tfile->ls();

    return;
}
int LimitShape::fillTree(string pathtotree)
{

    TDirectory * dir = (TDirectory*) tfile->Get(pathtotree.c_str());
    bool contains=(dir->GetListOfKeys()->Contains((systematic+"_"+shape_dist).c_str()));
    if (!contains){
    cout<<"Failure in obtaining tree"<<endl;
    return 0;
    }
    else{
    cout<<"trying to get "<<(systematic+"_"+shape_dist)<<endl;
    ttree = (TTree*)tfile->Get((pathtotree+"/"+systematic+"_"+shape_dist).c_str());
    cout<<"ttree entries "<<ttree->GetEntries()<<endl;
    return 1;
    }
}
void LimitShape::rescaleFinalweight(float num){ // this doesn't work memory leak problem 
    cout<<"must fill tree first"<<endl;
    Double_t finalweight;
    ttree->SetBranchAddress("finalweight",&finalweight);

    TTree * newtree = ttree->CloneTree(0);
    ttree->GetEntry(0);
    cout<<"finalweight before 1st entry "<<finalweight;
    for (int ent = 0; ent<ttree->GetEntries();ent++){

        ttree->GetEntry(ent);
        finalweight=finalweight*num;
        newtree->Fill();
    }   
    
    newtree->GetEntry(0);
    cout<<"finalweight after 1st entry "<<finalweight;
    ttree = newtree->CloneTree(0);
    
    
    return;
}
void LimitShape::setsystematic(string sys)
{
    systematic = sys;
    return;
}
void LimitShape::fillPDFs(string type,int order,string dist)
{
     
    shape_name = systematic+"_"+shape_dist+"_"+output;
    cout<< "shapename "<<shape_name<<endl;
    //RooRealVar * v = new RooRealVar("var","var",1.0, -10000.0, 10000.0);
    TList * tempvector_coeffs = new TList(); 
    TList * tempvector_finalcoeffs = new TList(); 

    if(type.compare("poly")==0){
        rangeMin = 18.0;
        rangeMax = 62.0;
        //vector<RooRealVar> tempvector_coeffs; 
        //vector<RooFormulaVar> tempvector_finalcoeffs; 

        for (int cnum=0;cnum < order+1; cnum++){
            tempvector_coeffs->Add(
                    new RooRealVar(
                    (TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+to_string(order)+"_"+systematic+"_"+output),
                    (TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+to_string(order)+"_"+systematic+"_"+output),
                    //(TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+systematic+"_"+output),
                    //(TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+systematic+"_"+output),
                                //1.0, -10000.0, 10000.0)
                                1.0, coeffMin[cnum], coeffMax[cnum])
                    );
        }
        for (int cnum=0;cnum < order+1; cnum++){
            tempvector_finalcoeffs->Add(
                new RooFormulaVar(
                (TString)("c"+to_string(cnum)+"_sq_"+shape_dist+"_"+to_string(order)+"_"+systematic+"_"+output),
                //(TString)("c"+to_string(cnum)+"_sq_"+shape_dist+"_"+systematic+"_"+output),
                "@0*@1",
                RooArgList(*((RooRealVar*)tempvector_coeffs->At(cnum)),*((RooRealVar*)tempvector_coeffs->At(cnum))))
                );
            } 

        //std::pair<std::map<char,int>::iterator,bool> ret;
        coeffs[type+to_string(order)] = tempvector_coeffs; //empty vector of type RooRealVar 

        finalcoeffs[type+to_string(order)] = tempvector_finalcoeffs;//empty vector of type RooRealVar 
        
        pdf[type+to_string(order)] =  new RooPolynomial(
                                        (systematic+"_"+shape_dist+"_"+output+"_fit_"+type+"_"+shape_dist).c_str(),
                                        (systematic+"_"+shape_dist+"_"+output+"_fit_"+type+"_"+shape_dist).c_str(),
                                        //(systematic+"_"+output+"_fit_"+type+"_"+shape_dist).c_str(),
                                        *poi,
                                        //util::unpack_caller(
                                        //RooArgList,tempvector_finalcoeffs.at(type+to_string(order)))
                                        RooArgList(*tempvector_finalcoeffs)
                    );

        //pdf.insert(
        //    pair<string,RooBernstein> (
        //        type+to_string(order),
        //        new RooBernstein(
        //            (TString)(systematic+"_"+output+"_fit_"+type+"_"+shape_dist),
        //            (TString)(systematic+"_"+output+"_fit_"+type+"_"+shape_dist),
        //            &poi,
        //            int an = util::unpack_caller(
        //                RooArgList,tempvector_finalcoeffs.at(type+to_string(order))
        //                )
        //            )
        //        )
        //    );
        }
    if(type.compare("bernstein")==0){
        rangeMin = 18.0;
        rangeMax = 62.0;
        //vector<RooRealVar> tempvector_coeffs; 
        //vector<RooFormulaVar> tempvector_finalcoeffs; 

        for (int cnum=0;cnum < order+1; cnum++){
            tempvector_coeffs->Add(
                    new RooRealVar(
                    (TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+systematic+"_"+output),
                    (TString)("c"+to_string(cnum)+"_"+shape_dist+"_"+systematic+"_"+output),
                                //1.0, -10000.0, 10000.0)
                                1.0, coeffMin[cnum], coeffMax[cnum])
                    );
        }
        for (int cnum=0;cnum < order+1; cnum++){
            tempvector_finalcoeffs->Add(
                new RooFormulaVar(
                (TString)("c"+to_string(cnum)+"_sq_"+shape_dist+"_"+systematic+"_"+output),
                "@0*@1",
                RooArgList(*((RooRealVar*)tempvector_coeffs->At(cnum)),*((RooRealVar*)tempvector_coeffs->At(cnum))))
                );
            } 

        //std::pair<std::map<char,int>::iterator,bool> ret;
        coeffs[type+to_string(order)] = tempvector_coeffs; //empty vector of type RooRealVar 

        finalcoeffs[type+to_string(order)] = tempvector_finalcoeffs;//empty vector of type RooRealVar 
        
        pdf[type+to_string(order)] =  new RooBernstein(
                                        (systematic+"_"+output+"_fit_"+type+"_"+shape_dist).c_str(),
                                        (systematic+"_"+output+"_fit_"+type+"_"+shape_dist).c_str(),
                                        *poi,
                                        //util::unpack_caller(
                                        //RooArgList,tempvector_finalcoeffs.at(type+to_string(order)))
                                        RooArgList(*tempvector_finalcoeffs)
                    );

        //pdf.insert(
        //    pair<string,RooBernstein> (
        //        type+to_string(order),
        //        new RooBernstein(
        //            (TString)(systematic+"_"+output+"_fit_"+type+"_"+shape_dist),
        //            (TString)(systematic+"_"+output+"_fit_"+type+"_"+shape_dist),
        //            &poi,
        //            int an = util::unpack_caller(
        //                RooArgList,tempvector_finalcoeffs.at(type+to_string(order))
        //                )
        //            )
        //        )
        //    );
        }
    if(type.compare("voigtian")==0){
        cout<<"mass num? "<<shape_dist.substr(1,2)<<endl;
        float mass = stoi((shape_dist.substr(1,2)));
        rangeMin = mass - 2.0;
        rangeMax = mass + 2.0;
        if(shape_vm["SigRanges"].as<int>()==1){
        if(mass==15){ rangeMin=14; rangeMax=25;}
        if(mass==20){ rangeMin=18; rangeMax=40;}
        if(mass==25){ rangeMin=18; rangeMax=40;}
        if(mass==30){ rangeMin=25; rangeMax=50;}
        if(mass==35){ rangeMin=25; rangeMax=45;}
        if(mass==40){ rangeMin=38; rangeMax=42;}
        if(mass==45){ rangeMin=35; rangeMax=55;}
        if(mass==50){ rangeMin=40; rangeMax=60;}
        if(mass==55){ rangeMin=40; rangeMax=66;}
        if(mass==60){ rangeMin=40; rangeMax=66;}
        poi->setRange(rangeMin,rangeMax);
        }
        if(shape_vm["SigCoeRanges"].as<int>()==1){
                        //vect0 = alpha vect1=sigma
        if(mass==15){ coeffvect[0]=0.15; coeffvect[1]=0.15;}
        if(mass==20){ coeffvect[0]=0.15; coeffvect[1]=0.15;}
        if(mass==25){ coeffvect[0]=0.18; coeffvect[1]=0.18;}
        if(mass==30){ coeffvect[0]=0.21; coeffvect[1]=0.21;}
        if(mass==35){ coeffvect[0]=0.24; coeffvect[1]=0.24;}
        if(mass==40){ coeffvect[0]=0.27; coeffvect[1]=0.27;}
        if(mass==45){ coeffvect[0]=0.30; coeffvect[1]=0.30;}
        if(mass==50){ coeffvect[0]=0.33; coeffvect[1]=0.33;}
        if(mass==55){ coeffvect[0]=0.36; coeffvect[1]=0.36;}
        if(mass==60){ coeffvect[0]=0.39; coeffvect[1]=0.39;}
        }
        shape_binning = 25;
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                            //stof(shape_dist.substr(1,2)),
                            //stof(shape_dist.substr(1,2))-2.0,
                            //stof(shape_dist.substr(1,2))+2.0)
                            0.5*(rangeMin+rangeMax),
                            rangeMin,
                            rangeMax)
                );
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("alpha_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("alpha_"+shape_dist+"_"+systematic+"_"+output),
                            0.1,
                            0.01,
                            //30.0)
                            coeffvect[0])
                            //8.0)
                            //0.3,
                            //0.1,
                            //0.5)
                );
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("sigma_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("sigma_"+shape_dist+"_"+systematic+"_"+output),
                            1.0,
                            0.01,
                            //30.0)
                            coeffvect[1])
                            //0.3,
                            //0.1,
                            //0.5)
                );
        coeffs[type+to_string(order)] = tempvector_coeffs; //empty vector of type RooRealVar 

        pdf[type+to_string(order)] = new RooVoigtian(
                                            (TString)("sigfit"),
                                            (TString)("sigfit"),
                                            *poi,
                                            *((RooRealVar*) tempvector_coeffs->At(0)),
                                            *((RooRealVar*) tempvector_coeffs->At(1)),
                                            *((RooRealVar*) tempvector_coeffs->At(2))
                                            );

        //pdf->insert(
        //    pair<string,RooGaussian> (
        //        (type+to_string(order)),
        //        new RooGaussian(
        //            (TString)("sigfit"),
        //            (TString)("sigfit"),
        //            &poi,
        //            tempvector_coeffs[0],
        //            tempvector_coeffs[1]
        //            )
        //        )
        //    );


    }
    if(type.compare("doublegaussian")==0){
        cout<<"mass num? "<<shape_dist.substr(1,2)<<endl;
        float mass = stoi((shape_dist.substr(1,2)));
        rangeMin = mass - 2.0;
        rangeMax = mass + 2.0;
        shape_binning = 25;
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                            stof(shape_dist.substr(1,2)),
                            stof(shape_dist.substr(1,2))-2.0,
                            stof(shape_dist.substr(1,2))+2.0)
                );
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("sigma_1_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("sigma_1_"+shape_dist+"_"+systematic+"_"+output),
                            1.0,
                            0.01,
                            //30.0)
                            coeffvect[0])
                            //0.3,
                            //0.1,
                            //0.5)
                );
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("sigma_2_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("sigma_2_"+shape_dist+"_"+systematic+"_"+output),
                            1.0,
                            0.01,
                            //30.0)
                            coeffvect[1])
                            //0.3,
                            //0.1,
                            //0.5)
                );
        coeffs[type+to_string(order)] = tempvector_coeffs; //empty vector of type RooRealVar 
        RooGaussModel * temppdf_1 = new RooGaussModel(
                                            (TString)("sigfit_1"),
                                            (TString)("sigfit_1"),
                                            *poi,
                                            *((RooRealVar*) tempvector_coeffs->At(0)),
                                            *((RooRealVar*) tempvector_coeffs->At(1))
                                            );
        RooGaussModel * temppdf_2 = new RooGaussModel(
                                            (TString)("sigfit_2"),
                                            (TString)("sigfit_2"),
                                            *poi,
                                            *((RooRealVar*) tempvector_coeffs->At(0)),
                                            *((RooRealVar*) tempvector_coeffs->At(2))
                                            );
        RooRealVar * frac_pdfs = new RooRealVar("frac_pdfs","fraction between gauss 1 and 2",0.5);
        RooRealVar * norm_1 = new RooRealVar("norm_1","events in sig_1",1,10000);
        RooRealVar * norm_2 = new RooRealVar("norm_2","events in sig_2",1,10000);

        pdf[type+to_string(order)] = new RooAddModel(
                                            (TString)("sigfit"),
                                            (TString)("sigfit"),
                                            RooArgList(*temppdf_1,*temppdf_2),
                                            *frac_pdfs
                                            );

        //pdf[type+to_string(order)] = new RooAddPdf(
        //                                    (TString)("sigfit"),
        //                                    (TString)("sigfit"),
        //                                    RooArgList(*temppdf_1,*temppdf_2),
        //                                    RooArgList(*norm_1, *norm_2)
        //                                    );


    }

    if(type.compare("gaussian")==0){
        cout<<"mass num? "<<shape_dist.substr(1,2)<<endl;
        float mass = stoi((shape_dist.substr(1,2)));
        rangeMin = mass - 2.0;
        rangeMax = mass + 2.0;
        shape_binning = 25;
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("mean_"+shape_dist+"_"+systematic+"_"+output),
                            stof(shape_dist.substr(1,2)),
                            stof(shape_dist.substr(1,2))-2.0,
                            stof(shape_dist.substr(1,2))+2.0)
                );
        tempvector_coeffs->Add(
                new RooRealVar(
                (TString)("sigma_"+shape_dist+"_"+systematic+"_"+output),
                (TString)("sigma_"+shape_dist+"_"+systematic+"_"+output),
                            0.5,
                            0.0,
                            2.0)
                            //0.3,
                            //0.1,
                            //0.5)
                );
        coeffs[type+to_string(order)] = tempvector_coeffs; //empty vector of type RooRealVar 

        pdf[type+to_string(order)] = new RooGaussian(
                                            (TString)("sigfit"),
                                            (TString)("sigfit"),
                                            *poi,
                                            *((RooRealVar*) tempvector_coeffs->At(0)),
                                            *((RooRealVar*) tempvector_coeffs->At(1))
                                            );

        //pdf->insert(
        //    pair<string,RooGaussian> (
        //        (type+to_string(order)),
        //        new RooGaussian(
        //            (TString)("sigfit"),
        //            (TString)("sigfit"),
        //            &poi,
        //            tempvector_coeffs[0],
        //            tempvector_coeffs[1]
        //            )
        //        )
        //    );


    }
    return;
}
void LimitShape::fillDataset(string cutstring)
{
    //if (!ttree) return;
    eventweight = new RooRealVar("finalweight","finalweight",-10000.0,10000.0);
    AMass = new RooRealVar("AMass","AMass",-10000.0,10000.0);
    m_vis = new RooRealVar("m_vis","m_vis",-10000.0,10000.0);
    //eventweight = new RooRealVar("finalweight","finalweight",0.0,10000.0);
    //data = new RooDataSet(shape_dist.c_str(),shape_dist.c_str(),RooArgSet(*poi),Import(*ttree));
    data = new RooDataSet(shape_dist.c_str(),shape_dist.c_str(),
                        RooArgSet(*poi,*eventweight,*AMass,*m_vis),
                        Import(*ttree),Cut(cutstring.c_str()),WeightVar("finalweight")
                    );
    cout<<"filled dataset, norm is "<<data->sumEntries()<<endl;
    norm = new RooRealVar((shape_dist+"_norm").c_str(),(shape_dist+"_norm").c_str(),data->sumEntries(),0.0,10*data->sumEntries());
    cout<<"normvar is : "<<norm->getVal()<<endl;
    
    return;
}
void LimitShape::fitToData()
{
    for(auto const& x: pdf){
        //used to be at line 482 in the fit ToData function
        cout<<"creating plotting spaces"<<endl;
        canvas[x.first] = new TCanvas(("canvas_"+shape_name).c_str(),("canvas_"+shape_name).c_str(),600,600);
        poi->setRange(rangeMin,rangeMax);
        plotframe[x.first] = poi->frame();
        plotframe[x.first]->SetAxisRange(rangeMin,rangeMax);
        cout<<"working on pdf "<<(x.first).c_str()<<endl;

        //data->plotOn(plotframe[x.first],Range(rangeMin,rangeMax));
        fitresults[x.first] = x.second->fitTo(
                                *data,Range(rangeMin,rangeMax), 
                                Minimizer("Minuit2","migrad"),
                                Save()
                            );
        //fitresults[x.first]->Print();
    }
    

    return;
}


void LimitShape::finalFitToData(string type,float sfr)
{
    //if(recursion>1000){return;}
    for(auto const& x: pdf){
        cout<<"working on pdf final fit"<<(x.first).c_str()<<endl;
        //plotframe[x.first] = poi->frame();
        //data->plotOn(plotframe[x.first],Binning(shape_binning));
        //TIter next(coeffs[x.first]); 
        //TObject * obj = 0; //while(obj = next()){
        //cout<<"recursion level "<<to_string(recursion)<<endl;
        //finalfitresults[x.first] = x.second->fitTo(
        //                        *data,Range(rangeMin,rangeMax), 
        //                        Minimizer("Minuit2","migrad"),
        //                        Save()
        //                    );

        //if(type.compare("gaussian")==0 || type.compare("voigtian")==0){
        if(type.compare("gaussian")==0 || type.compare("voigtian")==0 || type.compare("doublegaussian")==0){
        cout<<"running signal final fit"<<endl;
        for(const auto&& obj: *coeffs[x.first]){
            cout<<"object name "<<obj->GetName()<<endl; 
            cout<<"object class name "<<obj->ClassName()<<endl; 
            obj->Print();
            //lesson ... learned ...type  cast in paranthesis!!
            string varname = obj->GetName();
            cout<<"variable name "<<varname<<endl;
                if(varname.find("sigma") != string::npos)  {
                    cout<<"setting variable ranges for sigma"<<endl;
                    ((RooRealVar*) obj)->setRange( 
                        //((RooRealVar*)obj)->getValV()/2.0,
                        //((RooRealVar*)obj)->getValV()*2.0
                        ((RooRealVar*)obj)->getValV()/sfr,
                        ((RooRealVar*)obj)->getValV()*sfr
                        ); 
                    }
                }
            }
        //if(type.compare("Bkg")==0 || type.compare("irBkg")==0)//
        else{
        //if(type.compare("gaussian")!=0){
        for(const auto&& obj: *coeffs[x.first]){
                //cout<<"object name "<<obj->GetName()<<endl; 
                //cout<<"object class name "<<obj->ClassName()<<endl; 
                //obj->Print();
                //lesson ... learned ...type  cast in paranthesis!!
                //string varname = obj->GetName();
                //if(( abs(((RooRealVar*)obj)->getValV()) < abs(((RooRealVar*)obj)->getError())) &&(((RooRealVar*)obj)->getError()!=0.0) ){
                //cout<<"variable name "<<varname<<endl;
                    ((RooRealVar*) obj)->setRange( 
                        //((RooRealVar*)obj)->getValV()/1.50,
                        //((RooRealVar*)obj)->getValV()*1.50
                        ((RooRealVar*)obj)->getValV()/sfr,
                        ((RooRealVar*)obj)->getValV()*sfr
                        //((RooRealVar*)obj)->getValV()-1.0*recursion,
                        //((RooRealVar*)obj)->getValV()+1.0*recursion
                        ); 
                
                    //recursion++;
                    //LimitShape::finalFitToData(type);
                    
                 //   }
                //else{
                 //       cout<<"found best fit!!!!"<<endl;
                 //       cout<<"value "<<((RooRealVar*)obj)->getValV()<<endl;
                 //       cout<<"error "<<((RooRealVar*)obj)->getError()<<endl;
                 //   }
                }
                //return;
            }


        finalfitresults[x.first] = x.second->fitTo(
                                *data,Range(rangeMin,rangeMax), 
                                Minimizer("Minuit2","migrad"),
                                Save()
                            );
        finalfitresults[x.first]->Print();

    }
    

    return;
}
void LimitShape::printParamsNLLs()
{
    ofstream nllresults;
    nllresults.open(outputdir+"/nllresults_"+shape_dist+"_"+shape_name+"_"+systematic+".txt");
    nllresults << "fit results for "<<shape_name<<endl; 
    nllresults << "normalization "<<data->sumEntries()<<endl;
    //for(auto& [key,value]: coeffs){ this is the future .... c++17
    for(auto const& x: coeffs){
        nllresults<<"finalfit values for the coeff of "<<x.first<<endl;
        //for(RooRealVar* var: *value){
        for(const auto&& obj: *coeffs[x.first]){
            RooRealVar * var = (RooRealVar*) obj;
            cout<<var->GetName()<<"  "<<to_string(var->getValV()); 
            cout<<" err "<<to_string(var->getError())<<endl; 
            nllresults<<var->GetName()<<"  "<<to_string(var->getValV())<<" err "<<to_string(var->getError())<<endl; 
        }
    } 
    vector<float> nllvals;
    for(auto const& x: finalfitresults){
        nllresults<<"nll for type "<<x.first<<"    "<<x.second->minNll()<<endl;
        cout<<"nll for type "<<x.first<<"    "<<x.second->minNll()<<endl;
        nllvals.push_back(x.second->minNll());
        
    }
    int nllsize = (int) nllvals.size();
    for (int i=0; i<(nllsize-1);i++){
        cout<<"Delta nll between "<<to_string(i)<<" and "<<to_string(i+1)<<"     "<<to_string(nllvals[i] - nllvals[i+1])<<endl;
        nllresults<<"Delta nll between "<<to_string(i)<<" and "<<to_string(i+1)<<"     "<<to_string(nllvals[i] - nllvals[i+1])<<endl;
    }

    //printing chi2 and ndof
    for(auto const& x:plotframe){
        nllresults<<"chi square "<<x.second->chiSquare()<<endl;
        nllresults<<"ndof "<<coeffs[x.first]->GetSize()<<endl;
        cout<<"chi square "<<x.second->chiSquare()<<endl;
        cout<<"ndof "<<coeffs[x.first]->GetSize()<<endl;
    }
    return;
}
void LimitShape::createPlots(string type,int order,string dist)
{
    string key = type+to_string(order);//key for the map
    TCanvas * canvas = new TCanvas(("canvas_"+shape_name).c_str(),("canvas_"+shape_name).c_str(),600,600);
    canvas->cd();
    gPad->SetLeftMargin(0.15);
    bool doRatio = 0;
    TPaveText lumi=add_lumi(year,doRatio);
    TPaveText cms=add_CMS(doRatio);
    TPaveText pre=add_Preliminary(channel, doRatio);
    //TPaveText id = TPaveText(0.7,0.5,0.7+0.30,0.5+0.16,"NDC"); // mid right
    TPaveText id = TPaveText(0.3,0.7,0.3+0.30,0.7+0.16,"NDC");
    
    //calculating the corrected binned normalization - required for CMSSW 10 version of root
    float normcor = ((float)rangeMax-(float)rangeMin)/(float)shape_binning;

    //data->plotOn(plotframe[key],Range(rangeMin,rangeMax),Binning(shape_binning));
    //data->plotOn(plotframe[key],Range(rangeMin,rangeMax));

    cout<<"type order ? for plots "<<key<<endl;
    plotframe[key]->SetAxisRange(rangeMin,rangeMax);
    //plotframe[key]->SetTitle(("#mu#mu for "+dist+" "+output).c_str());
    plotframe[key]->SetTitle("");
    plotframe[key]->GetYaxis()->SetTitleOffset(1.6);
    plotframe[key]->GetXaxis()->SetTitle("M_{#mu#mu}");
    TGaxis::SetMaxDigits(2);
    //data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
    
    //pdf[key]->plotOn(plotframe[key],
    //        VisualizeError(*finalfitresults[key],1.0,kFALSE),
    //        Range(rangeMin,rangeMax),
    //        //DrawOption("F"),
    //        FillColor(kOrange)
    //        );
    //if(type.compare("gaussian")==0||type.compare("voigtian")==0){
    if(type.compare("gaussian")==0 || type.compare("voigtian")==0 || type.compare("doublegaussian")==0){
        id.SetBorderSize(   0 );
        id.SetFillStyle(    0 );
        id.SetTextAlign(   12 );
        id.SetTextColor(    1 );
        id.SetTextSize(0.04);
        id.SetTextFont (   42 );
        id.AddText((dist).c_str());

        //doesn't work with normalization?
        cout<<"normalization "<<normcor<<endl;
        //data->plotOn(plotframe[key],Range(rangeMin,rangeMax),Normalization(norm->getVal()*normcor));
        //data->plotOn(plotframe[key],Range(rangeMin,rangeMax),Normalization(norm->getVal()));
        data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
        gStyle->SetOptStat(0); 
        gStyle->SetOptFit(0); 
        cout<<"Plotting signal"<<endl;
        //pdf[key]->plotOn(plotframe[key],
        //        VisualizeError(*fitresults[key],1.0,kFALSE),
        //        Range(rangeMin,rangeMax),
        //        //DrawOption("F"),
        //        FillColor(kOrange)
        //        );
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax),Normalization(norm->getVal()/normcor));
        pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
        //data->plotOn(plotframe[key],Range(rangeMin,rangeMax),Normalization(norm->getVal()/normcor));
        //data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
        plotframe[key]->Draw();
        id.Draw();
        //data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
        plotframe[key]->Draw();
        lumi.Draw();
        cms.Draw();
        pre.Draw();
        cout<<"things drawn"<<endl;
        canvas->SaveAs((outputdir+"/"+to_string(order)+"_"+shape_name+"_"+systematic+".png").c_str()); 
        return;
    }
    if(shape_name.compare("SS_relaxed_data")==0){
        //gStyle->SetStatX(0.9);
        //gStyle->SetStatY(0.5);
        normcor = ((float)rangeMax-(float)rangeMin)/(float)shape_binning;
        cout<<"normalization "<<normcor<<endl;
        //data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax),Normalization(norm->getVal()));
        data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax));
        //cout<<"Plotting poly with normalization "<<norm->getVal()<<endl;
        //pdf[key]->plotOn(plotframe[key],
        //        VisualizeError(*finalfitresults[key],1.0,kFALSE),
        //        Range(rangeMin,rangeMax),
        //        Normalization(norm->getVal()),
        //        //DrawOption("F"),
        //        FillColor(kOrange)
        //        );
        //doesn't work with normalization?
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax),Normalization(norm->getVal()));
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax),Normalization(norm->getVal()/normcor));
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
        pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
        //cout<<"Plotted pdf "<<endl;
        //data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax),Normalization(norm->getVal()));
        //pdf[key]->paramOn(plotframe[key],Layout(0.2,0.7,0.8));
        //pdf[key]->paramOn(plotframe[key]);
        //plotframe[key]->getAttText()->SetTextSize(0.02);
        //data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
        cout<<"params set "<<endl;
        plotframe[key]->Draw("same");
        lumi.Draw("same");
        cms.Draw("same");
        pre.Draw("same");
        cout<<"things drawn"<<endl;
        canvas->SaveAs((outputdir+"/"+to_string(order)+"_"+shape_name+"_"+systematic+".png").c_str()); 
        cout<<"things saved"<<endl;
        return;
    }
    else{
        //gStyle->SetStatX(0.9);
        //gStyle->SetStatY(0.5);
        cout<<"normalization "<<normcor<<endl;
        //data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax),Normalization(norm->getVal()));
        data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax));
        //cout<<"Plotting poly with normalization "<<norm->getVal()<<endl;
        //pdf[key]->plotOn(plotframe[key],
        //        VisualizeError(*finalfitresults[key],1.0,kFALSE),
        //        Range(rangeMin,rangeMax),
        //        Normalization(norm->getVal()),
        //        //DrawOption("F"),
        //        FillColor(kOrange)
        //        );
        //doesn't work with normalization?
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax),Normalization(norm->getVal()/normcor));
        pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
        //pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
        //cout<<"Plotted pdf "<<endl;
        //data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax),Normalization(norm->getVal()));
        //pdf[key]->paramOn(plotframe[key],Layout(0.2,0.7,0.8));
        //pdf[key]->paramOn(plotframe[key]);
        //plotframe[key]->getAttText()->SetTextSize(0.02);
        //data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
        cout<<"params set "<<endl;
        plotframe[key]->Draw("same");
        lumi.Draw("same");
        cms.Draw("same");
        pre.Draw("same");
        cout<<"things drawn"<<endl;
        canvas->SaveAs((outputdir+"/"+to_string(order)+"_"+shape_name+"_"+systematic+".png").c_str()); 
        cout<<"things saved"<<endl;
        return;
    }
}








#endif
