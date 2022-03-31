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

/////////////////
//unpacking by vector 
//https://stackoverflow.com/questions/11044504/any-solution-to-unpack-a-vector-to-function-arguments-in-c
////////////////
#include <iostream>
#include <utility>
#include <vector>
#include <cassert>

//ignore RooFit mass output

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
        string           shape_name;
        string           shape_dist;
        string           type="shapetype";//RooFit Pdf class
        string           systematic="Nominal";//systematic for shapes

        TFile           *tfile;
        RooRealVar *poi;
        RooRealVar *norm;
        RooRealVar *eventweight;
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
        float coeffMin  = -10000.0;
        float coeffMax  =  10000.0;
        int maxscans  =  20;
        float shape_binning  = 16;
        int recursion = 0;

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
        void fillPDFs(string type,int order,string dist);
        void fillDataset();
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
  
  
    system("dir");

    //setting the vars
    poi = new RooRealVar("mll","m_{#mu#mu}", 18, 62);
    norm = new RooRealVar("norm",   "Normalization",1.0,0.0,10);
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
    tfile->ls();

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
                                1.0, coeffMin, coeffMax)
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
                                1.0, coeffMin, coeffMax)
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
void LimitShape::fillDataset()
{
    //if (!ttree) return;
    eventweight = new RooRealVar("finalweight","finalweight",-10000.0,10000.0);
    //eventweight = new RooRealVar("finalweight","finalweight",0.0,10000.0);
    //data = new RooDataSet(shape_dist.c_str(),shape_dist.c_str(),RooArgSet(*poi),Import(*ttree));
    data = new RooDataSet(shape_dist.c_str(),shape_dist.c_str(),
                        RooArgSet(*poi,*eventweight),
                        Import(*ttree),WeightVar("finalweight")
                    );
    cout<<"filled dataset dataset"<<endl;
    norm = new RooRealVar((shape_dist+"_norm").c_str(),(shape_dist+"_norm").c_str(),data->sumEntries(),0.0,10*data->sumEntries());
    
    return;
}
void LimitShape::fitToData()
{
    for(auto const& x: pdf){
        cout<<"working on pdf "<<(x.first).c_str()<<endl;
        canvas[x.first] = new TCanvas(("canvas_"+shape_name).c_str(),("canvas_"+shape_name).c_str(),600,600);
        poi->setRange(rangeMin,rangeMax);
        plotframe[x.first] = poi->frame();
        plotframe[x.first]->SetAxisRange(rangeMin,rangeMax);
        //data->plotOn(plotframe[x.first],Range(rangeMin,rangeMax),Binning(shape_binning));
        data->plotOn(plotframe[x.first],Range(rangeMin,rangeMax));
        fitresults[x.first] = x.second->fitTo(
                                *data,Range(rangeMin,rangeMax), 
                                Minimizer("Minuit2","migrad"),
                                Save()
                            );
        //fitresults[x.first]->Print();
    }
    

    return;
}

void LimitShape::scanFitToData(string type)
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    //ofstream scanfile;
    //scanfile.open(outputdir+"/fitscan_"+shape_dist+"_"+shape_name+"_"+systematic+".csv");
    TFile * scanfile = new TFile((outputdir+"/fitscan_"+shape_dist+"_"+shape_name+"_"+systematic+".root").c_str(),"recreate");
    scanfile->cd();
    
    //vector<float> 
    
    
    for(auto const& x: pdf){
        cout<<"working on pdf final fit"<<(x.first).c_str()<<endl;
        TTree * scantree = new TTree(("scantree_"+x.first).c_str(),("scantree_"+x.first).c_str());
        vector<vector<float>> vect_coeffval;
        vector<vector<float>> vect_coefferr;
        vector<vector<int>>   vect_scanmin1;
        vector<vector<int>>   vect_scanmax1;
        vector<vector<int>>   vect_scanmin2;
        vector<vector<int>>   vect_scanmax2;
        string coeffname;
        float coeffval;
        float coefferr;
        string coeffname2;
        float coeffval2;
        float coefferr2;
        int   scanmin1;
        int   scanmax1;
        int   scanmin2;
        int   scanmax2;
        int cc1 = 0;
        int cc2 = 0;
        scantree->Branch(("var_name")              ,&coeffname);
        scantree->Branch(("var_val")              ,&coeffval);
        scantree->Branch(("var_error")   ,&coefferr);
        scantree->Branch(("var_name2")              ,&coeffname2);
        scantree->Branch(("var_val2")              ,&coeffval2);
        scantree->Branch(("var_error2")   ,&coefferr2);
        scantree->Branch(("scanmin1"),&scanmin1);
        scantree->Branch(("scanmax1"),&scanmax1);
        scantree->Branch(("scanmin2"),&scanmin2);
        scantree->Branch(("scanmax2"),&scanmax2);

        //creating the branches based on the number of coeffs
        //for(const auto&& obj1: *coeffs[x.first]){
        //    vect_coeffval.push_back({}); 
        //    vect_coefferr.push_back({}); 
        //    vect_scanmin1.push_back({}); 
        //    vect_scanmax1.push_back({}); 
        //    vect_scanmin2.push_back({}); 
        //    vect_scanmax2.push_back({}); 
        //for(const auto&& obj2: *coeffs[x.first]){
        //    string varname = obj2->GetName();
        //    vect_coeffval[cc1][cc2]=0.0; 
        //    vect_coefferr[cc1][cc2]=0.0; 
        //    vect_scanmin1[cc1][cc2]=0.0; 
        //    vect_scanmax1[cc1][cc2]=0.0; 
        //    vect_scanmin2[cc1][cc2]=0.0; 
        //    vect_scanmax2[cc1][cc2]=0.0; 
        //    //scantree->Branch(varname.c_str()              ,&(coeffval[cc1][cc2]));
        //    //scantree->Branch((varname+"_error").c_str()   ,&(coefferr[cc1][cc2]));
        //    //scantree->Branch(("scanmin1_"+varname).c_str(),&(scanmin1[cc1][cc2]));
        //    //scantree->Branch(("scanmax1_"+varname).c_str(),&(scanmax1[cc1][cc2]));
        //    //scantree->Branch(("scanmin2_"+varname).c_str(),&(scanmin2[cc1][cc2]));
        //    //scantree->Branch(("scanmax2_"+varname).c_str(),&(scanmax2[cc1][cc2]));
        //    cc2++;
        //}
        //    cc1++;
        //}

        //conducting the scans
        cout<<"number of parameters in fit model "<<coeffs[x.first]->GetSize()<<endl;
        //if(type.compare("gaussian")!=0){
        int coeffcount = (int) coeffs[x.first]->GetSize();
        cc1 = 0;
        cc2 = 0;
        for(const auto&& obj1: *coeffs[x.first]){
            for(int scannum1=0; scannum1 < maxscans; scannum1++){
                    string varname = obj1->GetName();
                    ((RooRealVar*) obj1)->setRange( 
                        coeffMin+scannum1,//shinking the fit space
                        coeffMax-scannum1
                        ); 
                    coeffname = varname;
                    coeffval=((RooRealVar*)obj1)->getValV();
                    coefferr=((RooRealVar*)obj1)->getError();
                    scanmin1=coeffMin+scannum1;
                    scanmax1=coeffMax-scannum1;
            for(const auto&& obj2: *coeffs[x.first]){
                for(int scannum=0; scannum < maxscans; scannum++){
                    //cout<<"working on coeff "<<coeffcount<<endl;
                    string varname2 = obj2->GetName();
                    //cout<<"variable name "<<varname<<endl;
                    ((RooRealVar*) obj2)->setRange( 
                        coeffMin+scannum,//shinking the fit space
                        coeffMax-scannum
                        ); 
                    //plotframe[x.first] = poi->frame();
                    //data->plotOn(plotframe[x.first],Binning(shape_binning));
                    //finalfitresults[x.first] = x.second->fitTo(
                    //                        *data,Range(rangeMin,rangeMax), 
                    //                        Minimizer("Minuit2","migrad"),
                    //                        Save()
                    //                    );
                    x.second->fitTo(
                                            *data,Range(rangeMin,rangeMax), 
                                            Minimizer("Minuit2","migrad"),
                                            Save()
                                        );
                    //filling the tree
                    //vect_coeffval[cc1][cc2]=((RooRealVar*)obj2)->getValV();
                    //vect_coefferr[cc1][cc2]=((RooRealVar*)obj2)->getError();
                    //vect_scanmin1[cc1][cc2]=coeffMin+scannum1;
                    //vect_scanmax1[cc1][cc2]=coeffMax-scannum1;
                    //vect_scanmin2[cc1][cc2]=coeffMin+scannum;
                    //vect_scanmax2[cc1][cc2]=coeffMax-scannum;
                    coeffname2 = varname2;
                    coeffval2=((RooRealVar*)obj2)->getValV();
                    coefferr2=((RooRealVar*)obj2)->getError();
                    scanmin2=coeffMin+scannum;
                    scanmax2=coeffMax-scannum;
                    cout<<"coef "<<varname<<" val "<<((RooRealVar*)obj1)->getValV();
                    cout<<" err "<<((RooRealVar*)obj1)->getError()<<"  ";
                    cout<<"coef2 "<<varname2<<" val2 "<<((RooRealVar*)obj2)->getValV();
                    cout<<" err "<<((RooRealVar*)obj1)->getError()<<"  ";
                    cout<<" range1 "<<to_string(coeffMin+scannum1)<<" "<<to_string(coeffMax-scannum1);
                    cout<<" range2 "<<to_string(coeffMin+scannum)<<" "<<to_string(coeffMax-scannum)<<endl;
                    scantree->Fill();
                    
                    }//inner scan loop
                
                //scanfile<<shape_name+"_"+systematic+","<<x.first<<",";
                //int truefit=1;
                //for(const auto&& obj: *coeffs[x.first]){
                //    string varname = obj->GetName();
                //    //scanfile<<varname<<","<<coeffMin+1.0*scannum<<","<<coeffMax-1.0*scannum<<",";
                //    //scanfile<<((RooRealVar*)obj)->getValV()<<",";
                //    //scanfile<<((RooRealVar*)obj)->getError()<<",";
                //    if(( abs(((RooRealVar*)obj)->getValV()) > abs(((RooRealVar*)obj)->getError())) &&(((RooRealVar*)obj)->getError()!=0.0) ){
                //        //scanfile<<"1";
                //        truefit=truefit*1; 
                //        }
                //    else{
                //        //scanfile<<"0";
                //        truefit=truefit*0; 
                //        }
                //    //scanfile<<",";
                //    }//for loop for coeffs
                //scanfile<<to_string(truefit);
                //scanfile<<endl;

                cc2++;
                }//for loop coeff
                }//outer scan loop
                cc1++;
                }//for loop coeff

            //if(type.compare("gaussian")==0){
            //cout<<"running gaussian final fit"<<endl;
            //for(const auto&& obj: *coeffs[x.first]){
            //    cout<<"object name "<<obj->GetName()<<endl; 
            //    cout<<"object class name "<<obj->ClassName()<<endl; 
            //    obj->Print();
            //    //lesson ... learned ...type  cast in paranthesis!!
            //    string varname = obj->GetName();
            //    cout<<"variable name "<<varname<<endl;
            //        if(varname.find("sigma") != string::npos)  {
            //            cout<<"setting variable ranges for sigma"<<endl;
            //            ((RooRealVar*) obj)->setRange( 
            //                ((RooRealVar*)obj)->getValV()/2.0,
            //                ((RooRealVar*)obj)->getValV()*2.0
            //                ); 
            //            }
            //        }
            //    }
            //finalfitresults[x.first] = x.second->fitTo(
            //                        *data,Range(rangeMin,rangeMax), 
            //                        Minimizer("Minuit2","migrad"),
            //                        Save()
            //                    );
            //finalfitresults[x.first]->Print();
        //}//gaussian? 
    scantree->Write(scantree->GetName(),TObject::kOverwrite);
    }//loop over pdf

    scanfile->Write();
    scanfile->Close();
    return;
}
void LimitShape::recursiveFitToData(string type)
{
    if(recursion>((int)coeffMax)){cout<<"hit max recursion"<<endl; return;}
    for(auto const& x: pdf){
        cout<<"working on pdf final fit"<<(x.first).c_str()<<endl;
        plotframe[x.first] = poi->frame();
        //data->plotOn(plotframe[x.first],Binning(shape_binning));
        data->plotOn(plotframe[x.first]);
        TIter next(coeffs[x.first]); 
        TObject * obj = 0; //while(obj = next()){
        cout<<"recursion level "<<to_string(recursion)<<endl;
        finalfitresults[x.first] = x.second->fitTo(
                                *data,Range(rangeMin,rangeMax), 
                                Minimizer("Minuit2","migrad"),
                                Save()
                            );
        cout<<"number of parameters in fit model "<<coeffs[x.first]->GetSize()<<endl;
        if(type.compare("gaussian")!=0){
            int coeffcorrectcount = 0;
            int coeffcount = (int) coeffs[x.first]->GetSize();
            
            for(const auto&& obj: *coeffs[x.first]){
                cout<<"working on coeff "<<coeffcount<<endl;
                //cout<<"object name "<<obj->GetName()<<endl; 
                //cout<<"object class name "<<obj->ClassName()<<endl; 
                //obj->Print();
                //lesson ... learned ...type  cast in paranthesis!!
                string varname = obj->GetName();
                cout<<"variable name "<<varname<<endl;
                if(( abs(((RooRealVar*)obj)->getValV()) < abs(((RooRealVar*)obj)->getError())) &&(((RooRealVar*)obj)->getError()!=0.0) ){
                    cout<<"Error too large... shrinking fit"<<endl;
                    ((RooRealVar*) obj)->setRange( 
                        //((RooRealVar*)obj)->getValV()/1.50,
                        //((RooRealVar*)obj)->getValV()*1.50
                        //((RooRealVar*)obj)->getValV()-1.0*recursion,
                        //((RooRealVar*)obj)->getValV()+1.0*recursion
                        coeffMin+1.0*recursion,//shinking the fit space
                        coeffMax-1.0*recursion
                        ); 
                
                    
                    }
                else{
                        coeffcorrectcount+=1;
                        cout<<"found good value for parameter!!!!"<<endl;
                        cout<<"value "<<((RooRealVar*)obj)->getValV()<<endl;
                        cout<<"error "<<((RooRealVar*)obj)->getError()<<endl;
                    }
                }//for loop for coeffs

            if(coeffcorrectcount==coeffcount){
                cout<<"found best fit for all params!!!!"<<endl;
                cout<<"exiting recusion ... "<<endl;
                return;
                //break;
                }
            else{ //recurison block
                recursion++;
                LimitShape::recursiveFitToData(type);
                
                }
            }//if for gaussian

        //if(type.compare("gaussian")==0){
        //cout<<"running gaussian final fit"<<endl;
        //for(const auto&& obj: *coeffs[x.first]){
        //    cout<<"object name "<<obj->GetName()<<endl; 
        //    cout<<"object class name "<<obj->ClassName()<<endl; 
        //    obj->Print();
        //    //lesson ... learned ...type  cast in paranthesis!!
        //    string varname = obj->GetName();
        //    cout<<"variable name "<<varname<<endl;
        //        if(varname.find("sigma") != string::npos)  {
        //            cout<<"setting variable ranges for sigma"<<endl;
        //            ((RooRealVar*) obj)->setRange( 
        //                ((RooRealVar*)obj)->getValV()/2.0,
        //                ((RooRealVar*)obj)->getValV()*2.0
        //                ); 
        //            }
        //        }
        //    }
        //finalfitresults[x.first] = x.second->fitTo(
        //                        *data,Range(rangeMin,rangeMax), 
        //                        Minimizer("Minuit2","migrad"),
        //                        Save()
        //                    );
        //finalfitresults[x.first]->Print();
    }//loop over pdf

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
        if(type.compare("gaussian")!=0){
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

        if(type.compare("gaussian")==0){
        cout<<"running gaussian final fit"<<endl;
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

    return;
}
void LimitShape::createPlots(string type,int order,string dist)
{
    gPad->SetLeftMargin(0.15);
    string key = type+to_string(order);//key for the map
    cout<<"type order ? for plots "<<key<<endl;
    plotframe[key]->SetAxisRange(rangeMin,rangeMax);
    plotframe[key]->SetTitle(("#mu#mu for "+dist+" "+output).c_str());
    plotframe[key]->GetYaxis()->SetTitleOffset(1.6);
    plotframe[key]->GetXaxis()->SetTitle("M_{#mu#mu}");
    TGaxis::SetMaxDigits(2);
    //data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax));
    data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
    
    //pdf[key]->plotOn(plotframe[key],
    //        VisualizeError(*finalfitresults[key],1.0,kFALSE),
    //        Range(rangeMin,rangeMax),
    //        //DrawOption("F"),
    //        FillColor(kOrange)
    //        );
    if(type.compare("gaussian")==0){
    cout<<"Plotting gaussian"<<endl;
    pdf[key]->plotOn(plotframe[key],
            VisualizeError(*fitresults[key],0.4,kFALSE),
            Range(rangeMin,rangeMax),
            //DrawOption("F"),
            FillColor(kOrange)
            );
    }
    else{
    cout<<"Plotting poly?"<<endl;
    pdf[key]->plotOn(plotframe[key],
            //VisualizeError(*finalfitresults[key],1.0,kFALSE),
            Range(rangeMin,rangeMax),
            DrawOption("F"),
            FillColor(kOrange)
            );
    }
    pdf[key]->plotOn(plotframe[key],LineColor(kBlue),Range(rangeMin,rangeMax));
    //data->plotOn(plotframe[key],Binning(shape_binning),Range(rangeMin,rangeMax));
    data->plotOn(plotframe[key],Range(rangeMin,rangeMax));
    pdf[key]->paramOn(plotframe[key],data);
    plotframe[key]->getAttText()->SetTextSize(0.02);
    plotframe[key]->Draw("same");
    canvas[key]->SaveAs((outputdir+"/"+to_string(order)+"_"+shape_name+"_"+systematic+".png").c_str()); 
    return;
}








#endif
