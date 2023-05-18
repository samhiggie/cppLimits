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
#include "TList.h"
#include "TPad.h"

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
        ("fitmodeldistribution,fitmodeldistribution", po::value< string >(), "name of distribution in input file that you want to model")
        ("channel,channel", po::value< string >(), "channel")
        ("normalization,normalization", po::value< int >()->default_value(0), "run full optmizer with all orders for shapes")
        ("fitmodeldistributionNorm,fitmodeldistributionNorm", po::value< string >(), "if norm is run then name of distribution in input file that you want to use for a differenent norm")
        ("sfr,sfr", po::value< float >()->default_value(1.5), "second fit range based on multiple of nominal value")
        ("div,div", po::value< int >()->default_value(1000), "initial range on the coeffs for the fit")
        ("sfrdiv,sfrdiv", po::value< int >()->default_value(100), "how many scan points in the second fit range")
        ("coeRange,coeRange", po::value< float >()->default_value(10000), "initial range on the coeffs for the fit")
        ("coeRanges,coeRanges", po::value< vector<float>>()->multitoken(), "vector of coefficient ranges for more flexible fit")
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
    //float initialRange = vm["coeRange"].as<float>();
    vector<float> initialRange = vm["coeRanges"].as<vector<float>>();
    cout<<"ranges ";
    for(int coe=0; coe<initialRange.size(); coe++){
        cout<<to_string(initialRange[coe])<<" ";
    }
    cout<<endl;

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
    string outputdir = vm["output-dir"].as<string>();


    ofstream txtfile;
    txtfile.open("model_results_"+outputdir+"_sfr.csv");
    string fitmodeldistribution = vm["fitmodeldistribution"].as<string>();
    txtfile<<"#csv file containing data on fit results for "<< fitmodeldistribution <<endl;


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
        "Nominal"
    };

    cout<<"norms?"<<endl;
    int opt = vm["normalization"].as<int>();
    int order = vm["bkgo"].as<int>();
    vector<string> line;
    
    cout<<"starting to scan"<<endl;
    LimitShape * fitmodel = new LimitShape(vm,vm["fitmodeldistribution"].as<string>().c_str());
    fitmodel->setsystematic("Nominal");
    fitmodel->connectFile(inputFile);
    fitmodel->fillTree(channel+"_inclusive");
    if(opt==1){
        LimitShape * tempForNorm;
        tempForNorm = new LimitShape(vm,vm["fitmodeldistributionNorm"].as<string>().c_str());
        tempForNorm->setsystematic("Nominal");
        tempForNorm->connectFile(inputFile);
        tempForNorm->fillTree(channel+"_inclusive");
        tempForNorm->fillDataset("mll-m_vis>0&&"+AMasscut);
        cout<<"Normalization of SS relaxed "<<fitmodel->norm->getVal()<<endl;
        fitmodel->norm = tempForNorm->norm;
        cout<<"Normalization of fitmodel in signal region "<<fitmodel->norm->getVal()<<endl;
    }

    //for all the coefficients
    //vector<float> coeRangeVector = vm["coeRanges"].as<vector<float>>();
    //int iterations=0;
    //vector<float> coeRangeVectorBegin;

    //for(auto const& x:coeRangeVector){
    //    coeRangeVectorBegin.push_back(-1.0*x);
    //} 
    //fitmodel->coeRanges = coeRangeVectorBegin;
    //cout<<"Beginning Parameters "<<endl;
    //for(auto const& x:fitmodel->coeRanges){
    //    cout<<x<<",";
    //}
    //cout<<endl;

    for(int coe=0; coe<initialRange.size(); coe++){
        fitmodel->coeffMin.push_back( -1.0 * initialRange[coe]);
        fitmodel->coeffMax.push_back( 1.0 * initialRange[coe]);
    }
    //fitmodel->coeffMin = -1.0 * initialRange;
    //fitmodel->coeffMax = 1.0 * initialRange;
    fitmodel->fillPDFs(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
    fitmodel->fillDataset("mll-m_vis>0&&"+AMasscut);
    fitmodel->fitToData();
    fitmodel->finalFitToData(vm["fitmodeldistribution"].as<string>().c_str(),sfr);
    fitmodel->createPlots(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),"Bkg");
    fitmodel->printParamsNLLs();
    return 0;

    //setting up output datastructure
    //TFile * outScans = new TFile(("outScans_"+outputdir+"_sfr.root").c_str(),"RECREATE");
    ////outScans->cd();
    //TTree * outTree = new TTree("scanpts","scanpts");
    //vector<float> treeCoes;
    //vector<string> treeCoesNames;
    //vector<float> treeCoesErrors;
    //vector<string> treeCoesErrorsNames;
    //vector<int> treeMagCheck;
    //vector<string> treeMagCheckNames;
    //cout<<"Initialized tree"<<endl;
    //
    //for(int ov=0; ov<coeRangeVector.size(); ov++){
    //    treeCoes.push_back(0.0);
    //    treeCoesNames.push_back("coe_"+to_string(ov));
    //    treeCoesErrors.push_back(0.0);
    //    treeCoesErrorsNames.push_back("coeError_"+to_string(ov));
    //    treeMagCheck.push_back(0);
    //    treeMagCheckNames.push_back("coeError_"+to_string(ov));
    //}
    //cout<<"set initial values for vector vars "<<endl;
    //for(int ov=0; ov<coeRangeVector.size(); ov++){
    //    outTree->Branch((treeCoesNames[ov]).c_str(),&(treeCoes[ov]),(treeCoesNames[ov]+"/F").c_str());
    //    outTree->Branch((treeCoesErrorsNames[ov]).c_str(),&(treeCoesErrors[ov]),(treeCoesErrorsNames[ov]+"/F").c_str());
    //    outTree->Branch((treeMagCheckNames[ov]).c_str(),&(treeMagCheck[ov]),(treeMagCheckNames[ov]+"/I").c_str());
    //}
    //cout<<"Initialized branches "<<endl;
    
    
    //for(auto const& x:coeRangeVector) // add begin bracket











    /*




    for(int ov=0; ov<coeRangeVector.size(); ov++){
        for(int i=-1.0*coeRangeVector[ov]; i<coeRangeVector[ov]+1; i++){
            //for(auto const& y:coeRangeVector) // add begin bracket
            for(int iv=0; iv<coeRangeVector.size(); iv++){
                if(iv==ov){break;}
                for(int j=-1.0*coeRangeVector[iv];j<coeRangeVector[iv]+1; j++){
                    cout<<"c0 "<<ov<<" c1 "<<iv<<"        "<<i<<" "<<j<<endl;

                    
                    //fitmodel->coeRanges[ov] = i;
                    //fitmodel->coeRanges[iv] = j;
                    fitmodel->coeffMin[ov] = i;
                    fitmodel->coeffMin[iv] = j;
                    fitmodel->coeffMax[ov] = i;
                    fitmodel->coeffMax[iv] = j;

                    cout<<"coeRanges for iteration "<<iterations<<"   ";
                    for(auto const& x:fitmodel->coeRanges){
                        cout<<x<<",";
                    }
                    cout<<endl;
            

                    fitmodel->fillPDFs(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),vm["fitmodeldistribution"].as<string>().c_str());
                    fitmodel->fillDataset("mll-m_vis>0&&"+AMasscut);

                    fitmodel->fitToData();
                    fitmodel->finalFitToData(vm["fitmodeldistribution"].as<string>().c_str(),sfr);

                    for(auto const& x:fitmodel->coeffs){
                    //map<string, TList*>::iterator mapit;
                    //for(mapit = fitmodel->coeffs.begin(); mapit != fitmodel->coeffs.end(); mapit++){
                        for(const auto&& obj: *fitmodel->coeffs[x.first]){
                        //TList * t = mapit->second;
                        //TIter tlistiter(t);
                        //while (TObject *obj = tlistiter()) // add begin bracket
                            RooRealVar * var = (RooRealVar*) obj;

                            //txtfile<<var->GetName()<<"  "<<to_string(var->getValV())<<" err "<<to_string(var->getError())<<endl; 
                            cout<<var->GetName()<<"  "<<to_string(var->getValV())<<" err "<<to_string(var->getError())<<endl; 

                            //for output tree
                            //treeCoesNames[ov] = var->GetName();
                            //treeCoes[ov] = var->getValV();
                            //treeCoesErrors[ov] = var->getError();
                    
                            
                            line.push_back(var->GetName());
                            line.push_back(to_string(var->getValV()));
                            line.push_back(to_string(var->getError()));
                            line.push_back("error less? ");
                        
                            if(var->getError()<var->getValV()){
                                line.push_back("1");
                                //treeMagCheck[ov] = 1;
                            }
                            else{
                                line.push_back("0");
                                //treeMagCheck[ov] = 1;
                            }
                            //outTree->Fill();
                        }
                        
                    }

                    vector<float> nllvals;
                    for(auto const& x: fitmodel->finalfitresults){
                        //txtfile<<"nll for type "<<x.first<<"    "<<x.second->minNll()<<endl;
                        line.push_back("nll for type "+x.first);
                        line.push_back(to_string(x.second->minNll()));
                        nllvals.push_back(x.second->minNll());
                        
                    }
                    int nllsize = (int) nllvals.size();
                    //for (int i=0; i<(nllsize-1);i++){
                    //    txtfile<<"Delta nll between "<<to_string(i)<<" and "<<to_string(i+1)<<"     "<<to_string(nllvals[i] - nllvals[i+1])<<endl;
                    //}

                    for(auto  const& subline:line){
                        txtfile<<subline<<",";
                        cout<<subline<<",";
                    }
                    //cout<<line<<endl;
                    txtfile<<endl;
                    cout<<endl;
                    line.clear();
            
        
                    iterations++;
                    cout<<"iteration "<<to_string(iterations)<<endl;
                }
            }
        }
    }







    //initialize
    //for(int cnum=0; cnum < order+1; cnum++){
    //    fitmodel->coeRanges.push_back(vm["div"].as<int>());
    //}

    ////editing internal coeRange to be most flexible in fit 

    //for(int cnumout=0; cnumout < order+1; cnumout++){
    //for(int cnumin=0; cnumin < order+1; cnumin++){


    //for(int i=vm["div"].as<int>(); i<(initialRange + vm["div"].as<int>()) ; i = i + vm["div"].as<int>() ){
    //for(int k=vm["div"].as<int>(); k<(initialRange + vm["div"].as<int>()) ; k =  + vm["div"].as<int>() ){


    //    fitmodel->coeRanges[cnumout]=i;
    //    fitmodel->coeRanges[cnumin]=k;

    //    sfr = vm["sfr"].as<float>();

    //    for(int j=0; j<vm["sfrdiv"].as<int>()  ;j++ ){
    //        sfr = sfr*(1.0 - 1.0/((float) vm["sfrdiv"].as<int>()));

    //        cout<<"iteration "<<to_string(i)<<" "<<to_string(j)<<endl;


    //        //fitmodel->coeffMin = -1.0 * i;
    //        //fitmodel->coeffMax =  1.0 * i;
    //        //fitmodel->coeffMin = -1.0 * initialRange;
    //        //fitmodel->coeffMax =  1.0 * initialRange;
    //        for(int cn=0; cn < order+1; cn++){
    //            
    //            //fitmodel->coeRanges[cn]=i;

    //            line.push_back("coe");
    //            line.push_back(to_string(cn));
    //            line.push_back("range");
    //            line.push_back(to_string(fitmodel->coeRanges[cn]));
    //        }

    //        line.push_back(to_string(sfr));
    //        //line.push_back(to_string(fitmodel->coeffMin));
    //        //line.push_back(to_string(fitmodel->coeffMax));


    //        fitmodel->fillPDFs(vm["bkgt"].as<string>(),vm["bkgo"].as<int>(),vm["fitmodeldistribution"].as<string>().c_str());
    //        fitmodel->fillDataset("mll-m_vis>0&&"+AMasscut);

    //        fitmodel->fitToData();
    //        fitmodel->finalFitToData(vm["fitmodeldistribution"].as<string>().c_str(),sfr);
    //        for(auto const& x:fitmodel->coeffs){
    //            for(const auto&& obj: *fitmodel->coeffs[x.first]){
    //                RooRealVar * var = (RooRealVar*) obj;
    //                //txtfile<<var->GetName()<<"  "<<to_string(var->getValV())<<" err "<<to_string(var->getError())<<endl; 
    //                line.push_back(var->GetName());
    //                line.push_back(to_string(var->getValV()));
    //                line.push_back(to_string(var->getError()));
    //                line.push_back("error less? ");
    //                if(var->getError()<var->getValV()){
    //                    line.push_back("1");
    //                }
    //                else{
    //                    line.push_back("0");
    //                }
    //            }
    //        }

    //        vector<float> nllvals;
    //        for(auto const& x: fitmodel->finalfitresults){
    //            //txtfile<<"nll for type "<<x.first<<"    "<<x.second->minNll()<<endl;
    //            line.push_back("nll for type "+x.first);
    //            line.push_back(to_string(x.second->minNll()));
    //            nllvals.push_back(x.second->minNll());
    //            
    //        }
    //        int nllsize = (int) nllvals.size();
    //        //for (int i=0; i<(nllsize-1);i++){
    //        //    txtfile<<"Delta nll between "<<to_string(i)<<" and "<<to_string(i+1)<<"     "<<to_string(nllvals[i] - nllvals[i+1])<<endl;
    //        //}

    //        for(auto  const& subline:line){
    //            txtfile<<subline<<",";
    //        }
    //        txtfile<<endl;
    //        line.clear();
    //    }
    //}
    //}
    //fitmodel->coeRanges.clear();
    //for(int cnum=0; cnum < order+1; cnum++){
    //    fitmodel->coeRanges[cnum]=vm["div"].as<int>();
    //}
    //}
    //}


    //fitmodel->finalFitToData(vm["fitmodeldistribution"].as<string>().c_str(),sfr);
    //fitmodel->scanFitToData(vm["fitmodeldistribution"].as<string>().c_str());
    //fitmodel->recursiveFitToData(vm["fitmodeldistribution"].as<string>().c_str());

    //fitmodel->printParamsNLLs();
    //fitmodel->createPlots(vm["bkgt"].as<string>(),1,vm["fitmodeldistribution"].as<string>().c_str());
    //fitmodel->createPlots(vm["bkgt"].as<string>(),2,vm["fitmodeldistribution"].as<string>().c_str());
    //fitmodel->createPlots(vm["bkgt"].as<string>(),3,vm["fitmodeldistribution"].as<string>().c_str());
    //fitmodel->createPlots(vm["bkgt"].as<string>(),4,vm["fitmodeldistribution"].as<string>().c_str());
    //fitmodel->createPlots(vm["bkgt"].as<string>(),5,vm["fitmodeldistribution"].as<string>().c_str());
    //fitmodel->createPlots(vm["bkgt"].as<string>(),6,vm["fitmodeldistribution"].as<string>().c_str());
    //end comment

    cout<<"done with FF ..."<<endl;
    cout<<"number of iterations "<<to_string(iterations)<<endl;
    txtfile.close();

    //outScans->cd();
    //outTree->Write();
    //outScans->Write();
    //outScans->Close();


    return 0;



    */
}
