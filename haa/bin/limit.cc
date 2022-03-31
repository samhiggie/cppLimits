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
    gErrorIgnoreLevel = kWarning;
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

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
        ("sfr,sfr", po::value< float >()->default_value(1.5), "second fit range based on multiple of nominal value")
        ("sfrIr,sfrIr", po::value< float >()->default_value(1.5), "second fit range based on multiple of nominal value for the irreducible background")
        ("coeRange,coeRange", po::value< float >()->default_value(10000), "initial range on the coeffs for the fit")
        ("coeRangeIr,coeRangeIr", po::value< float >()->default_value(10000), "initial range on the coeffs for the irreducible background fit")
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

    vector<string> systematics = { "Nominal","scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
               "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
               "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
               "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"};

    int opt = vm["runs"].as<int>();
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
        data->fillDataset();

        map<string,LimitShape*> ZZshapes;
        for(auto const& sys: systematics){
            ZZshapes[sys] = new LimitShape(vm,"irBkg");
            ZZshapes[sys]->setsystematic("Nominal");
            ZZshapes[sys]->connectFile(inputFile);
            ZZshapes[sys]->fillTree(channel+"_inclusive");
            ZZshapes[sys]->coeffMin = -1.0 * initialRangeIr;
            ZZshapes[sys]->coeffMax = 1.0 * initialRangeIr;
            ZZshapes[sys]->fillPDFs(vm["irbkgt"].as<string>(),vm["irbkgo"].as<int>(),"irBkg");
            ZZshapes[sys]->fillDataset();
            ZZshapes[sys]->fitToData();
            ZZshapes[sys]->finalFitToData("irBkg",sfrIr);
            ZZshapes[sys]->printParamsNLLs();
            ZZshapes[sys]->createPlots(vm["irbkgt"].as<string>(),vm["irbkgo"].as<int>(),"irBkg");
        }    
        
        cout<<"done with ZZ ..."<<endl;
        systematics = { "Nominal" };

        map<string,LimitShape*> FFshapes;
        for(auto const& sys: systematics){
            FFshapes[sys] = new LimitShape(vm,"Bkg");
            FFshapes[sys]->setsystematic("Nominal");
            FFshapes[sys]->connectFile(inputFile);
            FFshapes[sys]->fillTree(channel+"_inclusive");
            FFshapes[sys]->coeffMin = -1.0 * initialRange;
            FFshapes[sys]->coeffMax = 1.0 * initialRange;
            FFshapes[sys]->fillPDFs(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
            FFshapes[sys]->fillDataset();
            FFshapes[sys]->fitToData();
            FFshapes[sys]->finalFitToData("Bkg",sfr);
            FFshapes[sys]->printParamsNLLs();
            FFshapes[sys]->createPlots(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
        }    
         
        cout<<"done with FF ..."<<endl;

        systematics = { "Nominal","scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
                   "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
                   "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
                   "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"};

        //vector<LimitShape*> signals;
        map<string,LimitShape*> signals;
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
                signals[sys+"_"+x] = new LimitShape(vm,x);
                //auto s = sh.second;
                mass = stoi((signals[sys+"_"+x]->shape_dist).substr(1,2)); 
                signals[sys+"_"+x]->setsystematic(sys);
                signals[sys+"_"+x]->connectFile(inputFile);
                signals[sys+"_"+x]->fillTree(channel+"_inclusive");
                signals[sys+"_"+x]->fillPDFs("gaussian",mass,signals[sys+"_"+x]->shape_dist);
                signals[sys+"_"+x]->fillDataset();
                signals[sys+"_"+x]->fitToData();
                signals[sys+"_"+x]->finalFitToData("gaussian",2.0);
                signals[sys+"_"+x]->printParamsNLLs();
                signals[sys+"_"+x]->createPlots("gaussian",mass,signals[sys+"_"+x]->shape_dist);
                off->loadShape(sys+"_"+signals[sys+"_"+x]->shape_dist,signals[sys+"_"+x]);
                if (sys=="Nominal"){
                    nominal_signals[x] = new LimitShape(vm,x);
                    nominal_signals[x]->setsystematic(sys);
                    nominal_signals[x]->connectFile(inputFile);
                    nominal_signals[x]->fillTree(channel+"_inclusive");
                    nominal_signals[x]->fillPDFs("gaussian",mass,nominal_signals[x]->shape_dist);
                    nominal_signals[x]->fillDataset();
                    nominal_signals[x]->fitToData();
                    nominal_signals[x]->finalFitToData("gaussian",2.0);
                    nominal_signals[x]->printParamsNLLs();
                    nominal_signals[x]->createPlots("gaussian",mass,nominal_signals[x]->shape_dist);
                }
            }
        } 

        //loading the shapes 
        
        for(auto const& sys: systematics){
            off->loadShape(sys+"_irBkg",ZZshapes[sys]);
        }
        off->loadShape("Nominal_Bkg",FFshapes["Nominal"]);


        off->createTxtfile();
        off->printDatacard(false);

        off->wsp->import(*(data->data));

        //for(auto const& sys: systematics){
        //    off->wsp->import(*(ZZshapes[sys]->pdf[vm["irbkgt"].as<string>()+to_string(vm["irbkgo"].as<int>())]));
        //}

        ZZshapes["Nominal"]->pdf[vm["irbkgt"].as<string>()+to_string(vm["irbkgo"].as<int>())]->SetName(("irBkg_"+channel).c_str());
        FFshapes["Nominal"]->pdf[vm["bkgt"].as<string>()+to_string(vm["bkgo"].as<int>())]->SetName(("Bkg_"+channel).c_str());

        off->wsp->import(*(ZZshapes["Nominal"]->pdf[vm["irbkgt"].as<string>()+to_string(vm["irbkgo"].as<int>())]));
        off->wsp->import(*(FFshapes["Nominal"]->pdf[vm["bkgt"].as<string>()+to_string(vm["bkgo"].as<int>())]));
        off->interpolateParameters(nominal_signals); //this will import signal spline function to workspace

        off->wsp->writeToFile(
            ((off->outputdir)+"/"+"HToAAWorkspace_full_"+(off->output)+".root").c_str()
            );

        //for(auto const& pdf: FF){
        //    off->wsp->import(pdf.second);
        //}
        //for(auto const& pdf: ZZ){
        //    off->wsp->import(pdf.second);
        //}
    }
    else{
    
        LimitShape * ZZ = new LimitShape(vm,"irBkg");
        ZZ->setsystematic("Nominal");
        ZZ->connectFile(inputFile);
        ZZ->fillTree(channel+"_inclusive");
        ZZ->coeffMin = -1.0 * initialRangeIr;
        ZZ->coeffMax = 1.0 * initialRangeIr;
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),1,"irBkg");
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),2,"irBkg");
        ZZ->fillPDFs(vm["irbkgt"].as<string>(),3,"irBkg");
        //ZZ->fillPDFs(vm["irbkgt"].as<string>(),4,"irBkg");
        //ZZ->fillPDFs(vm["irbkgt"].as<string>(),5,"irBkg");
        //ZZ->fillPDFs(vm["irbkgt"].as<string>(),6,"irBkg");
        ZZ->fillDataset();
        ZZ->fitToData();
        ZZ->finalFitToData("irBkg",sfrIr);
        //ZZ->scanFitToData("irBkg");
        ZZ->printParamsNLLs();
        //ZZ->createPlots(vm["irbkgt"].as<string>(),1,"irBkg");
        //ZZ->createPlots(vm["irbkgt"].as<string>(),2,"irBkg");
        //ZZ->createPlots(vm["irbkgt"].as<string>(),3,"irBkg");
        //ZZ->createPlots(vm["irbkgt"].as<string>(),4,"irBkg");
        //ZZ->createPlots(vm["irbkgt"].as<string>(),5,"irBkg");
        //ZZ->createPlots(vm["irbkgt"].as<string>(),6,"irBkg");

        LimitShape * FF = new LimitShape(vm,"Bkg");
        FF->setsystematic("Nominal");
        FF->connectFile(inputFile);
        FF->fillTree(channel+"_inclusive");
        FF->coeffMin = -1.0 * initialRange;
        FF->coeffMax = 1.0 * initialRange;
        FF->fillPDFs(vm["bkgt"].as<string>(),1,"Bkg");
        FF->fillPDFs(vm["bkgt"].as<string>(),2,"Bkg");
        FF->fillPDFs(vm["bkgt"].as<string>(),3,"Bkg");
        //FF->fillPDFs(vm["bkgt"].as<string>(),4,"Bkg");
        //FF->fillPDFs(vm["bkgt"].as<string>(),5,"Bkg");
        //FF->fillPDFs(vm["bkgt"].as<string>(),6,"Bkg");
        FF->fillDataset();
        FF->fitToData();
        FF->finalFitToData("Bkg",sfr);
        //FF->scanFitToData("Bkg");
        //FF->recursiveFitToData("Bkg");
        FF->printParamsNLLs();
        //FF->createPlots(vm["bkgt"].as<string>(),1,"Bkg");
        //FF->createPlots(vm["bkgt"].as<string>(),2,"Bkg");
        //FF->createPlots(vm["bkgt"].as<string>(),3,"Bkg");
        //FF->createPlots(vm["bkgt"].as<string>(),4,"Bkg");
        //FF->createPlots(vm["bkgt"].as<string>(),5,"Bkg");
        //FF->createPlots(vm["bkgt"].as<string>(),6,"Bkg");

        LimitShape * signal = new LimitShape(vm,"a40");
        signal->setsystematic("Nominal");
        signal->connectFile(inputFile);
        signal->fillTree(channel+"_inclusive");
        signal->fillPDFs("gaussian",40,"a40");
        signal->fillDataset();
        signal->fitToData();
        signal->finalFitToData("gaussian",2.0);
        signal->printParamsNLLs();
        signal->createPlots("gaussian",40,"a40");

    }
   

    return 0;
}
