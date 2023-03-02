#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include <boost/program_options.hpp>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "time.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"

///////////////////////////////////////////
//
//supporting header files 
///////////////////////////////////////////

#include "cppLimits/haa/interface/LimitShape.h"
#include "cppLimits/haa/interface/Office.h"



int main(int ac, char* av[]) {

    //namespaces
    using namespace std;
    using namespace ROOT; 
    using namespace RooFit;
    namespace po = boost::program_options;
    gErrorIgnoreLevel = kInfo;
    //RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    // command line options? 
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "help message")
        ("quiet,q", "suppress messages to stdout")
        ("output-file,o",  po::value< std::string >(), "output file path: /path/to/file")
        ("input-file", po::value< std::string >(), "Path to data file: /path/to/data/file")
        ("output-dir,d", po::value< std::string >(), "Path to data file: /path/to/data/file")
        ("irbkgo,irbkgo", po::value< int >(), "order poly")
        ("irbkgt,irbkgt", po::value< string >(), "type")
        ("bkgo,bkgo", po::value< int >(), "order bkg poly")
        ("bkgt,bkgt", po::value< string >(), " bkg type")
        ("channel,channel", po::value< string >(), "channel")
        ("runs,runs", po::value< int >()->default_value(0), "run full optmizer with all orders for shapes")
        ("basic,basic", po::value< int >()->default_value(0), "create the datacards and RooWorkspaces in the old 2016 way")
        ("sfr,sfr", po::value< float >()->default_value(1.5), "second fit range based on multiple of nominal value")
        ("sfrIr,sfrIr", po::value< float >()->default_value(1.5), "second fit range based on multiple of nominal value for the irreducible background")
        ("coeRange,coeRange", po::value< float >()->default_value(10000), "initial range on the coeffs for the fit")
        ("coeRangeIr,coeRangeIr", po::value< float >()->default_value(10000), "initial range on the coeffs for the irreducible background fit")
        ("sigt,sigt", po::value< string >(), " signal type gaussian, voigtian, double gauss")
        ("coeRangeSig,coeRangeSig", po::value< vector<float>>()->multitoken(), "initial range on the alpha and sigma")
        ("SigRanges,SigRanges", po::value< int >()->default_value(0), "per signal template ranges for poi")
        ("SigCoeRanges,SigCoeRanges", po::value< int >()->default_value(0), "per signal template ranges for poi")
        ("plots,plots", po::value< int >()->default_value(0), "create pltos?")
        ("ssrelaxed,ssrelaxed", po::value< int >()->default_value(0), "use Same Sign Relaxed Dataset for Bkg-datadriven shape")
        ("masses,masses", po::value< int >()->default_value(0), "create datacards for all mass points")
        ;

    po::positional_options_description p;
    p.add("input-file", 1);
    po::variables_map vm;
    //p.add("output-file", 1);

    po::store(po::command_line_parser(ac, av).
            options(desc).positional(p).run(), vm);
    po::notify(vm);
    std::string inputFile;
    float sfr = vm["sfr"].as<float>();
    float initialRange = vm["coeRange"].as<float>();
    float sfrIr = vm["sfrIr"].as<float>();
    float initialRangeIr = vm["coeRangeIr"].as<float>();
    vector<float> initialRangeSig = vm["coeRangeSig"].as<vector<float>>();

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("input-file"))
    {
        
        inputFile = vm["input-file"].as<std::string>();
        cout << "Input files are: " 
             << inputFile << "\n";
    }
    string channel = vm["channel"].as<string>();
    string output = vm["output-file"].as<string>();
    string AMasscut = "AMass<120.0";
    if((channel).find("mmmt")!=string::npos){AMasscut = "AMass<120.0";}
    if((channel).find("mmtt")!=string::npos){AMasscut = "AMass<130.0";}
    if((channel).find("mmet")!=string::npos){AMasscut = "AMass<120.0";}
    if((channel).find("mmem")!=string::npos){AMasscut = "AMass<110.0";}

    //vector<string> systematics = { "Nominal","scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
    //           "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
    //           "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
    //           "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"};
    vector<string> systematics = {
    "Nominal",
    "scale_tUp",
    "scale_mUp",
    "scale_eUp",
    //mmtt
    "scale_t_1prong_TauLoose_MuoVLoose_EleVLooseUp",
    "scale_t_1prong1pizero_TauLoose_MuoVLoose_EleVLooseUp",
    "scale_t_3prong_TauLoose_MuoVLoose_EleVLooseUp",
    "scale_t_3prong1pizero_TauLoose_MuoVLoose_EleVLooseUp",
    //mmmt
    "scale_t_1prong_TauLoose_MuoTight_EleVLooseUp",
    "scale_t_1prong1pizero_TauLoose_MuoTight_EleVLooseUp",
    "scale_t_3prong_TauLoose_MuoTight_EleVLooseUp",
    "scale_t_3prong1pizero_TauLoose_MuoTight_EleVLooseUp",
    //mmet
    "scale_t_1prong_TauMedium_MuoVLoose_EleVTightUp",
    "scale_t_1prong1pizero_TauMedium_MuoVLoose_EleVTightUp",
    "scale_t_3prong_TauMedium_MuoVLoose_EleVTightUp",
    "scale_t_3prong1pizero_TauMedium_MuoVLoose_EleVTightUp",
    "scale_tDown",
    "scale_mDown",
    "scale_eDown",
    //mmtt
    "scale_t_1prong_TauLoose_MuoVLoose_EleVLooseDown",
    "scale_t_1prong1pizero_TauLoose_MuoVLoose_EleVLooseDown",
    "scale_t_3prong_TauLoose_MuoVLoose_EleVLooseDown",
    "scale_t_3prong1pizero_TauLoose_MuoVLoose_EleVLooseDown",
    //mmmt
    "scale_t_1prong_TauLoose_MuoTight_EleVLooseDown",
    "scale_t_1prong1pizero_TauLoose_MuoTight_EleVLooseDown",
    "scale_t_3prong_TauLoose_MuoTight_EleVLooseDown",
    "scale_t_3prong1pizero_TauLoose_MuoTight_EleVLooseDown",
    //mmet
    "scale_t_1prong_TauMedium_MuoVLoose_EleVTightDown",
    "scale_t_1prong1pizero_TauMedium_MuoVLoose_EleVTightDown",
    "scale_t_3prong_TauMedium_MuoVLoose_EleVTightDown",
    "scale_t_3prong1pizero_TauMedium_MuoVLoose_EleVTightDown"
    };

    int opt = vm["runs"].as<int>();
    int old = vm["basic"].as<int>();
    
    if(opt==0){

        LimitShape * data = new LimitShape(vm,"data_obs");
        data->setsystematic("Nominal");

        int foundtree = 0;
        data->connectFile(inputFile);
        foundtree = data->fillTree(channel+"_inclusive");

        if(foundtree==0){
            cout<<"Likely that data obs is empty ... might want to loosen cuts? idk"<<endl;
            cout<<"saving Bkg as data_obs for now"<<endl;
            data->shape_dist = "Bkg";
            data->connectFile(inputFile);
            data->fillTree(channel+"_inclusive");
            cout<<"renaming the distribution back to data_obs"<<endl;
        }
        data->shape_dist = "data_obs";
        data->fillDataset("mll-m_vis>0&&"+AMasscut);

        map<string,LimitShape*> ZZshapes;
        for(auto const& sys: systematics){
            ZZshapes[sys] = new LimitShape(vm,"irBkg");
            ZZshapes[sys]->setsystematic(sys);
            ZZshapes[sys]->connectFile(inputFile);
            ZZshapes[sys]->fillTree(channel+"_inclusive");
            ZZshapes[sys]->coeffMin = -1.0 * initialRangeIr;
            ZZshapes[sys]->coeffMax = 1.0 * initialRangeIr;
            ZZshapes[sys]->fillPDFs(vm["irbkgt"].as<string>(),vm["irbkgo"].as<int>(),"irBkg");
            ZZshapes[sys]->fillDataset("mll-m_vis>0&&"+AMasscut);
            //if(sys.compare("Nominal")==0){
            ZZshapes[sys]->fitToData();
            ZZshapes[sys]->finalFitToData("irBkg",sfrIr);
            if(vm["plots"].as<int>()==1){
            ZZshapes[sys]->createPlots(vm["irbkgt"].as<string>(),vm["irbkgo"].as<int>(),"irBkg");
            ZZshapes[sys]->printParamsNLLs();
            //}
            }
        }    
        
        cout<<"done with ZZ ..."<<endl;
        systematics = { "Nominal" };

        map<string,LimitShape*> FFshapes;
        LimitShape * tempForNorm;
        for(auto const& sys: systematics){
            if(vm["plots"].as<int>()==0){
                FFshapes[sys] = new LimitShape(vm,"Bkg");
                FFshapes[sys]->setsystematic("Nominal");
                FFshapes[sys]->connectFile(inputFile);
                FFshapes[sys]->fillTree(channel+"_inclusive");
                FFshapes[sys]->coeffMin = -1.0 * initialRange;
                FFshapes[sys]->coeffMax = 1.0 * initialRange;
                FFshapes[sys]->fillPDFs(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
                FFshapes[sys]->fillDataset("mll-m_vis>0&&"+AMasscut);
                FFshapes[sys]->fitToData();
                FFshapes[sys]->finalFitToData("Bkg",sfr);
            if(vm["plots"].as<int>()==1){
                FFshapes[sys]->createPlots(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
                FFshapes[sys]->printParamsNLLs();
            }
            }
            else{
                //get normalization from Signal Region 
                tempForNorm = new LimitShape(vm,"Bkg");
                tempForNorm->setsystematic("Nominal");
                tempForNorm->connectFile(inputFile);
                tempForNorm->fillTree(channel+"_inclusive");
                tempForNorm->fillDataset("((mll-m_vis>0)&&"+AMasscut+")");
                tempForNorm->norm->setVal(tempForNorm->data->sumEntries());
                cout<<" Normalization from tight Bkg "<<tempForNorm->norm->getVal()<<endl;


                FFshapes[sys] = new LimitShape(vm,"SS_relaxed_data");
                FFshapes[sys]->setsystematic("Nominal");
                FFshapes[sys]->connectFile(inputFile);
                FFshapes[sys]->fillTree(channel+"_inclusive");
                FFshapes[sys]->coeffMin = -1.0 * initialRange;
                FFshapes[sys]->coeffMax = 1.0 * initialRange;
                FFshapes[sys]->fillPDFs(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
                FFshapes[sys]->fillDataset("mll-m_vis>0&&"+AMasscut);
                cout<<" Normalization from loose Bkg "<<FFshapes[sys]->data->sumEntries()<<endl;
                FFshapes[sys]->norm = tempForNorm->norm;
                cout<<" Correct Normalization ? "<<FFshapes[sys]->norm->getVal()<<endl;
                FFshapes[sys]->fitToData();
                FFshapes[sys]->finalFitToData("SS_relaxed_data",sfr);
            if(vm["plots"].as<int>()==1){
                FFshapes[sys]->createPlots(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
                FFshapes[sys]->printParamsNLLs();
            }
            }
        }    
         
        cout<<"done with FF ..."<<endl;

        //systematics = { "Nominal","scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
        //           "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
        //           "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
        //           "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"};
        systematics = {
        "Nominal",
        "scale_tUp",
        "scale_mUp",
        "scale_eUp",
        //mmtt
        "scale_t_1prong_TauLoose_MuoVLoose_EleVLooseUp",
        "scale_t_1prong1pizero_TauLoose_MuoVLoose_EleVLooseUp",
        "scale_t_3prong_TauLoose_MuoVLoose_EleVLooseUp",
        "scale_t_3prong1pizero_TauLoose_MuoVLoose_EleVLooseUp",
        //mmmt
        "scale_t_1prong_TauLoose_MuoTight_EleVLooseUp",
        "scale_t_1prong1pizero_TauLoose_MuoTight_EleVLooseUp",
        "scale_t_3prong_TauLoose_MuoTight_EleVLooseUp",
        "scale_t_3prong1pizero_TauLoose_MuoTight_EleVLooseUp",
        //mmet
        "scale_t_1prong_TauMedium_MuoVLoose_EleVTightUp",
        "scale_t_1prong1pizero_TauMedium_MuoVLoose_EleVTightUp",
        "scale_t_3prong_TauMedium_MuoVLoose_EleVTightUp",
        "scale_t_3prong1pizero_TauMedium_MuoVLoose_EleVTightUp",
        "scale_tDown",
        "scale_mDown",
        "scale_eDown",
        //mmtt
        "scale_t_1prong_TauLoose_MuoVLoose_EleVLooseDown",
        "scale_t_1prong1pizero_TauLoose_MuoVLoose_EleVLooseDown",
        "scale_t_3prong_TauLoose_MuoVLoose_EleVLooseDown",
        "scale_t_3prong1pizero_TauLoose_MuoVLoose_EleVLooseDown",
        //mmmt
        "scale_t_1prong_TauLoose_MuoTight_EleVLooseDown",
        "scale_t_1prong1pizero_TauLoose_MuoTight_EleVLooseDown",
        "scale_t_3prong_TauLoose_MuoTight_EleVLooseDown",
        "scale_t_3prong1pizero_TauLoose_MuoTight_EleVLooseDown",
        //mmet
        "scale_t_1prong_TauMedium_MuoVLoose_EleVTightDown",
        "scale_t_1prong1pizero_TauMedium_MuoVLoose_EleVTightDown",
        "scale_t_3prong_TauMedium_MuoVLoose_EleVTightDown",
        "scale_t_3prong1pizero_TauMedium_MuoVLoose_EleVTightDown"
        };



        //vector<LimitShape*> signals;
        map<string,LimitShape*> signal;
        map<string,map<string,LimitShape*>> systematicsignals;
        map<string,LimitShape*> nominal_signals;
        //map<<string,map<string,LimitShapes*>> allsignals;
        vector<string> masses = {
                //"a15",
                "a20",
                "a25",
                "a30",
                "a35",
                "a40",
                "a45",
                "a50",
                "a55",
                "a60"
                };


        Office * off = new Office(vm); 


        int mass;
        for(auto const& x: masses){
            for(auto const& sys: systematics){
                cout<<"working on signal "<<x<<endl;
                //signals.push_back(new LimitShape(vm,x));
                signal[x] = new LimitShape(vm,x);
                //auto s = sh.second;
                mass = stoi((signal[x]->shape_dist).substr(1,2)); 
                signal[x]->coeffvect = initialRangeSig;
                signal[x]->setsystematic(sys);
                signal[x]->connectFile(inputFile);
                signal[x]->fillTree(channel+"_inclusive");
                //signal[x]->rescaleFinalweight(10.0); // this doesn't work... memory leak problems
                //signal[x]->fillPDFs("gaussian",mass,signal[x]->shape_dist);
                signal[x]->fillPDFs(vm["sigt"].as<string>(),mass,signal[x]->shape_dist);
                signal[x]->fillDataset("mll-m_vis>0&&"+AMasscut);
                //don't need to fit the systematics... they are rate parameters 
                signal[x]->fitToData();
                //signal[x]->finalFitToData("gaussian",2.0);
                signal[x]->finalFitToData(vm["sigt"].as<string>(),2.0);
                if(vm["plots"].as<int>()==1){
                    signal[x]->createPlots(vm["sigt"].as<string>(),mass,signal[x]->shape_dist);
                }
                signal[x]->printParamsNLLs();

                off->loadShape(sys+"_"+signal[x]->shape_dist,signal[x]);
                systematicsignals[sys][x] = signal[x];
                if (sys=="Nominal"){
                    nominal_signals[x] = new LimitShape(vm,x);
                    nominal_signals[x]->coeffvect = initialRangeSig;
                    nominal_signals[x]->setsystematic(sys);
                    nominal_signals[x]->connectFile(inputFile);
                    nominal_signals[x]->fillTree(channel+"_inclusive");
                    //nominal_signals[x]->rescaleFinalweight(10.0); // this doesn't work... memory leak problems
                    //nominal_signals[x]->fillPDFs("gaussian",mass,nominal_signals[x]->shape_dist);
                    nominal_signals[x]->fillPDFs(vm["sigt"].as<string>(),mass,nominal_signals[x]->shape_dist);
                    nominal_signals[x]->fillDataset("mll-m_vis>0&&"+AMasscut);
                    nominal_signals[x]->fitToData();
                    //nominal_signals[x]->finalFitToData("gaussian",2.0);
                    nominal_signals[x]->finalFitToData(vm["sigt"].as<string>(),2.0);
                    if(vm["plots"].as<int>()==1){
                    //nominal_signals[x]->createPlots("gaussian",mass,nominal_signals[x]->shape_dist);
                    nominal_signals[x]->createPlots(vm["sigt"].as<string>(),mass,nominal_signals[x]->shape_dist);
                    nominal_signals[x]->printParamsNLLs();
                    }
                    off->loadShape(sys+"_"+nominal_signals[x]->shape_dist,nominal_signals[x]);
                }
            }
        } 

        //loading the shapes 
        
        for(auto const& sys: systematics){
            off->loadShape(sys+"_irBkg",ZZshapes[sys]);
            //off->interpolateParametersSystematics(vm["sigt"].as<string>(),systematicsignals[sys],sys); //this will import signal spline function to workspace
        }
        off->loadShape("Nominal_Bkg",FFshapes["Nominal"]);




        off->wsp->import(*(data->data));

        //for(auto const& sys: systematics){
        //    off->wsp->import(*(ZZshapes[sys]->pdf[vm["irbkgt"].as<string>()+to_string(vm["irbkgo"].as<int>())]));
        //}

        ZZshapes["Nominal"]->pdf[vm["irbkgt"].as<string>()+to_string(vm["irbkgo"].as<int>())]->SetName(("irBkg_"+output).c_str());
        FFshapes["Nominal"]->pdf[vm["bkgt"].as<string>()+to_string(vm["bkgo"].as<int>())]->SetName(("Bkg_"+output).c_str());
        ZZshapes["Nominal"]->norm->SetName(("irBkg_"+output+"_norm").c_str());
        FFshapes["Nominal"]->norm->SetName(("Bkg_"+output+"_norm").c_str());

        off->wsp->import(*(ZZshapes["Nominal"]->pdf[vm["irbkgt"].as<string>()+to_string(vm["irbkgo"].as<int>())]));
        off->wsp->import(*(FFshapes["Nominal"]->pdf[vm["bkgt"].as<string>()+to_string(vm["bkgo"].as<int>())]));

        //off->interpolateParameters("gaussian",nominal_signals); //this will import signal spline function to workspace
        //off->interpolateParameters(vm["sigt"].as<string>(),nominal_signals); //this will import signal spline function to workspace
        //off->interpolateParameters(vm["sigt"].as<string>(),nominal_signals); //this will import signal spline function to workspace
        off->interpolateParameters2016Style(vm["sigt"].as<string>(),nominal_signals); //this will import signal spline function to workspace
        if(old==0){
        off->wsp->import(*(ZZshapes["Nominal"]->norm));
        off->wsp->import(*(FFshapes["Nominal"]->norm));
        }

        off->wsp->writeToFile(
            ((off->outputdir)+"/"+"HToAAWorkspace_full_"+(off->output)+".root").c_str()
            );
        off->decorrelateParameters("Nominal_irBkg",vm["irbkgt"].as<string>()+to_string(vm["irbkgo"].as<int>()));
        off->decorrelateParameters("Nominal_Bkg",vm["bkgt"].as<string>()+to_string(vm["bkgo"].as<int>()));
        off->wsp->writeToFile(
            ((off->outputdir)+"/"+"HToAAWorkspace_decorrelated_"+(off->output)+".root").c_str()
            );

        off->createTxtfile();
        //off->printDatacard(false);
        off->print2016StyleDatacard(false);
        if(vm["masses"].as<int>()==1){
            off->createTxtfilePerMass();
            //off->printDatacardPerMass(false);
            off->print2016StyleDatacardPerMass(false);
        }
        //for(auto const& pdf: FF){
        //    off->wsp->import(pdf.second);
        //}
        //for(auto const& pdf: ZZ){
        //    off->wsp->import(pdf.second);
        //}
    }
    else{
    
        //commenting out other templates 

        LimitShape * ZZ = new LimitShape(vm,"irBkg");
        ZZ->setsystematic("Nominal");
        ZZ->connectFile(inputFile);
        ZZ->fillTree(channel+"_inclusive");
        ZZ->coeffMin = -1.0 * initialRangeIr;
        ZZ->coeffMax = 1.0 * initialRangeIr;
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),1,"irBkg");
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),2,"irBkg");
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),3,"irBkg");
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),4,"irBkg");
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),5,"irBkg");
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),6,"irBkg");
        ZZ->fillDataset("mll-m_vis>0&&"+AMasscut);
        ZZ->fitToData();
        ZZ->finalFitToData("irBkg",sfrIr);
        //ZZ->scanFitToData("irBkg");
        ZZ->printParamsNLLs();
        ZZ->createPlots(vm["irbkgt"].as<string>(),1,"irBkg");
        ZZ->createPlots(vm["irbkgt"].as<string>(),2,"irBkg");
        ZZ->createPlots(vm["irbkgt"].as<string>(),3,"irBkg");
        ZZ->createPlots(vm["irbkgt"].as<string>(),4,"irBkg");
        ZZ->createPlots(vm["irbkgt"].as<string>(),5,"irBkg");
        ZZ->createPlots(vm["irbkgt"].as<string>(),6,"irBkg");
        
        cout<<"done with ZZ ..."<<endl;

        LimitShape * tempForNorm;
        tempForNorm = new LimitShape(vm,"Bkg");
        tempForNorm->setsystematic("Nominal");
        tempForNorm->connectFile(inputFile);
        tempForNorm->fillTree(channel+"_inclusive");
        tempForNorm->fillDataset("mll-m_vis>0&&"+AMasscut);


        LimitShape * FF = new LimitShape(vm,"SS_relaxed_data");
        FF->setsystematic("Nominal");
        FF->connectFile(inputFile);
        FF->fillTree(channel+"_inclusive");
        FF->coeffMin = -1.0 * initialRange;
        FF->coeffMax = 1.0 * initialRange;
        FF->fillPDFs(vm["bkgt"].as<string>(),1,"SS_relaxed_data");
        FF->fillPDFs(vm["bkgt"].as<string>(),2,"SS_relaxed_data");
        FF->fillPDFs(vm["bkgt"].as<string>(),3,"SS_relaxed_data");
        FF->fillPDFs(vm["bkgt"].as<string>(),4,"SS_relaxed_data");
        FF->fillPDFs(vm["bkgt"].as<string>(),5,"SS_relaxed_data");
        FF->fillPDFs(vm["bkgt"].as<string>(),6,"SS_relaxed_data");
        FF->fillDataset("mll-m_vis>0&&"+AMasscut);
        cout<<"Normalization of SS relaxed "<<FF->norm->getVal()<<endl;
        FF->norm = tempForNorm->norm;
        cout<<"Normalization of FF in signal region "<<FF->norm->getVal()<<endl;
        FF->fitToData();
        FF->finalFitToData("SS_relaxed_data",sfr);
        //FF->scanFitToData("SS_relaxed_data");
        //FF->recursiveFitToData("SS_relaxed_data");
        FF->printParamsNLLs();
        FF->createPlots(vm["bkgt"].as<string>(),1,"SS_relaxed_data");
        FF->createPlots(vm["bkgt"].as<string>(),2,"SS_relaxed_data");
        FF->createPlots(vm["bkgt"].as<string>(),3,"SS_relaxed_data");
        FF->createPlots(vm["bkgt"].as<string>(),4,"SS_relaxed_data");
        FF->createPlots(vm["bkgt"].as<string>(),5,"SS_relaxed_data");
        FF->createPlots(vm["bkgt"].as<string>(),6,"SS_relaxed_data");

        cout<<"done with FF ..."<<endl;
        //end comment

        LimitShape * signal = new LimitShape(vm,"a40");
        signal->setsystematic("Nominal");
        signal->connectFile(inputFile);

        //not ready yet
        signal->coeffvect = initialRangeSig;

        signal->fillTree(channel+"_inclusive");
        //signal->fillPDFs("gaussian",40,"a40");
        signal->fillPDFs(vm["sigt"].as<string>(),40,"a40");
        signal->fillDataset("mll-m_vis>0&&"+AMasscut);
        signal->fitToData();
        //signal->finalFitToData("gaussian",2.0);
        signal->finalFitToData(vm["sigt"].as<string>(),2.0);
        signal->printParamsNLLs();
        //signal->createPlots("gaussian",40,"a40");
        signal->createPlots(vm["sigt"].as<string>(),40,"a40");

        LimitShape * signal_a20 = new LimitShape(vm,"a20");
        signal_a20->setsystematic("Nominal");
        signal_a20->connectFile(inputFile);

        //not ready yet
        signal_a20->coeffvect = initialRangeSig;

        signal_a20->fillTree(channel+"_inclusive");
        //signal_a20->fillPDFs("gaussian",20,"a20");
        signal_a20->fillPDFs(vm["sigt"].as<string>(),20,"a20");
        signal_a20->fillDataset("mll-m_vis>0&&"+AMasscut);
        signal_a20->fitToData();
        //signal_a20->finalFitToData("gaussian",2.0);
        signal_a20->finalFitToData(vm["sigt"].as<string>(),2.0);
        signal_a20->printParamsNLLs();
        //signal_a20->createPlots("gaussian",20,"a20");
        signal_a20->createPlots(vm["sigt"].as<string>(),20,"a20");

        LimitShape * signal_a60 = new LimitShape(vm,"a60");
        signal_a60->setsystematic("Nominal");
        signal_a60->connectFile(inputFile);
        //not ready yet
        signal_a60->coeffvect = initialRangeSig;

        signal_a60->fillTree(channel+"_inclusive");
        //signal_a60->fillPDFs("gaussian",60,"a60");
        signal_a60->fillPDFs(vm["sigt"].as<string>(),60,"a60");
        signal_a60->fillDataset("mll-m_vis>0&&"+AMasscut);
        signal_a60->fitToData();
        //signal_a60->finalFitToData("gaussian",2.0);
        signal_a60->finalFitToData(vm["sigt"].as<string>(),2.0);
        signal_a60->printParamsNLLs();
        //signal_a60->createPlots("gaussian",60,"a60");
        signal_a60->createPlots(vm["sigt"].as<string>(),60,"a60");

    }
   

    return 0;
}
