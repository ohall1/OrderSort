//http://www-numi.fnal.gov/offline_software/srt_public_context/WebDocs/Companion/cxx_crib/objects.html

// My code include.
#include "DataSource.h"
#include "Unpacker.h"
#include "Calibrator.h"
#include "Analysis.h"

// ROOT include.
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TROOT.h>
//#include <TStopwatch.h> //Create objects to use as benchmark tests of performance

// C++ include.
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <ctime>
#include <sys/time.h>

// ** Includes for system interaction ** //
//#include <signal.h>
//#include "TSystem.h"
//#include "kbhit.h"
//#include <TSysEvtHandler.h>



void Usage(char *progname) {
    fprintf(stderr,"This will at some point become a useage instruction...");
    fprintf(stderr,"Usage\n%s \n\t-O (for online data analysis, source_option=2)\n",progname);
    fprintf(stderr,"\t-i [id=0] (DataSpy id for online data)\n");
    fprintf(stderr,"\t-F [RunName] (for data analysis from file, source_option=1)\n");
    fprintf(stderr,"\t-d [DataDir] (directory path for data file)\n");
    fprintf(stderr,"\t-x (enable data Xfer to remote DataSink, options in ./config directory)\n");
    fprintf(stderr,"\t-U [unpacker_option=0] (Unpacker output level; 0: no ouput, 1: histograms, 2: TTree, 3: histograms+TTree)\n");
    fprintf(stderr,"\t-C [calibrator_option=1] (Calibrator output level; 0: no ouput, 1: histograms, 2: TTree, 3: histograms+TTree)\n");
    fprintf(stderr,"\t-A [analysis_option=1] (Analysis output level; 0: no ouput, 1: histograms, 2: TTree, 3: histograms+TTree)\n");
    fprintf(stderr,"\t-w [time_window=3202] (time_window for event clustering)\n");
    fprintf(stderr,"\t-R (enable writing output to Root file)\n");
    fprintf(stderr,"\t-L [RLrun=-1, RLfirst, RLnum] (loop over list of MIDAS data files, overrides option -F if RLrun>0)\n");
    fprintf(stderr,"\t-v (verbose mode)\n");
    fprintf(stderr,"\t-o [Output file name]\n");

    exit(1);
}

double getRealTime() {
  struct timeval time;
  gettimeofday(&time, NULL);
  return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}

double getCPUtime() {
  return (double)clock() / CLOCKS_PER_SEC;
}

int main  (int argc, char **argv){
  
  double start_t = getRealTime();
  double start_clk = getCPUtime();
  
  int i,j;

  std::string RunName;
  std::string userOutFile;
  bool b_userOutFile = false;

  bool b_verbose=false;

  int time_window= 3202;

  int source_option=2; //1: file, 2: online

  int unpacker_option=0;
  int calibrator_option=1; //3: histograms+tree, 1: only histograms, 2: only tree
  int analysis_option=1; //3; //check implamantation of TTree bit!!!

  int id= 0;

  bool b_Xfer= false;

  std::string DataDir = "/scratch2/ohall/AIDA/data/RIBF128/AIDA/Raw/";  //"/home/data/rootFiles/";// Input data directory.
  std::string OutDir = "/scratch2/ohall/AIDA/data/RIBF128/AIDA/root/";  //"/home/data/rootFiles/"; 

  bool b_root_file= false;

  int RLrun= -1; int RLfirst= 0; int RLnum= 0;

  int loops;


  // Loop to read command line arguments
  if (argc >1) {

    std::cout << " command line arguments: " <<argc<< std::endl;

    for(i=1;i <argc;i++) {

      std::cout << argv[i] << std::endl;

      if ( (argv[i][0] == '-') || (argv[i][0] == '/') ) {
	
	switch( argv[i][1] ) {

	  //analyze online  data
        case 'O':  
	  source_option = 2;
	  break;

	  //analyze data from file
	case 'F': 
	  RunName=  argv[++i];
	  source_option = 1;
	  std::cout << "  F= " << RunName << std::endl;
	  break;

	case 'i':
	  id = atoi(argv[++i]);
	  break;

	case 'x':
	  b_Xfer = true;
	  break;

	case 'd':
	  DataDir=  argv[++i];
	  break;

	case 'U':
	  j = atoi(argv[++i]); std::cout << "  U= " << j << std::endl;
	  if(j>=0 && j<=3) unpacker_option= j;
	  else{
	    std::cout << " Invalid unpacker_option: -U " << j << std::endl;
	    Usage(argv[0]);
	  }
	  break;


	case 'C':
	  j = atoi(argv[++i]); std::cout << "  C= " << j << std::endl;
	  if(j>=0 && j<=3) calibrator_option= j;
	  else{
	    std::cout << " Invalid calibrator_option: -C " << j << std::endl;
	    Usage(argv[0]);
	  }
	  break;

	case 'A':
	  j = atoi(argv[++i]); std::cout << "  A= " << j << std::endl;
	  if(j>=0 && j<=3) analysis_option= j;
	  else{
	    std::cout << " Invalid analysis_option: -A " << j << std::endl;
	    Usage(argv[0]);
	  }
	  break;

	case 'w':
	  time_window = atoi(argv[++i]);
	  break;

	case 'R':
	  b_root_file = true;
	  break;

	case 'L':
	  RLrun = atoi(argv[++i]);
	  RLfirst = atoi(argv[++i]);
	  RLnum = atoi(argv[++i]);
	  source_option = 1;
	  break;

	case 'v':
	  b_verbose = true;
	  break;

	case 'o':
	  userOutFile = argv[++i];
	  b_userOutFile = true;
	  std::cout << userOutFile <<std::endl;
	  break;
 	
	default:
	  Usage(argv[0]);
	  break;
	}
      }
    }
  }


  /* System interaction prompt
  TApplication *rootapp= new TApplication("AIDAonline",&argc,argv);
  gSystem->AddSignalHandler( new TSignalHandler( kSigInterrupt ) );
  gSystem->AddSignalHandler( new TSignalHandler( kSigTermination ) );
  */

  std::cout << "AIDA be like: in da F11!!" << std::endl;

  std::string FileNameData="";

  char *FileCalibParameters= (char*)"config/parameters_N82_test.txt";
  std::string FileNameRoot=""; //= "aida_sort_" + RunName +".root";
  TFile * fMain; //= new TFile(FileNameRoot.data(),"RECREATE");

  char RLname[64];

  //open Root output file (TTrees, histograms, etc...)
  if(b_root_file){
    if(source_option==2){
      FileNameRoot = OutDir+"aida_sort_online.root";
    }
    else if(RLrun>0){
      if(!b_userOutFile){	
      	sprintf(RLname,"R%d_list",RLrun);   //*********Change here for AIDA/R run prefixes******************//
      	FileNameRoot = OutDir+"aida_sort_"+ std::string(RLname) +".root";
      }
      else if(b_userOutFile){
      	FileNameRoot = OutDir +userOutFile;
      	std::cout <<"here"<<std::endl;
      }
    }
    else{
    	if(!b_userOutFile){
    		FileNameRoot = OutDir+"aida_sort_" + RunName +".root";
    	}
    	else if(b_userOutFile){
    		FileNameRoot = OutDir+userOutFile;
    	}
    }
    fMain = new TFile(FileNameRoot.data(),"RECREATE");
  }


  unsigned long n_entry=0;
  //  unsigned long n_write= 10000000;
  unsigned long n_update= 100000000;

  std::cout << " declaring my class objects"<<std::endl;
  DataSource midas_data;
  Unpacker unpacker_data;

  Calibrator calibrator_data;
  Analysis analysis_data;

  std::cout << " initializing things"<<std::endl;

  unpacker_data.InitUnpacker(unpacker_option, FileCalibParameters);
  /*if(calibrator_option || analysis_option)*/ calibrator_data.InitCalibrator(calibrator_option, FileCalibParameters);
  /*if(analysis_option)*/ analysis_data.InitAnalysis(analysis_option, FileCalibParameters);

  analysis_data.SetEventTimeWindow(1.*time_window); 
  std::cout << " *** ANALYSIS: TIME WINDOW=  " << analysis_data.GetEventTimeWindow() << std::endl;

  if(b_verbose){
    std::cout << "------------------------------ DEBUG ENABLED------------------------"<<std::endl;
    unpacker_data.SetBDebug(true);
    calibrator_data.SetBDebug(true);
    analysis_data.SetBDebug(true);
  }
  else{
    unpacker_data.SetBDebug(false);
    calibrator_data.SetBDebug(false);
    analysis_data.SetBDebug(false);
  }

  std::cout<< " +++ Histograming:\n"
	   << "                   U: "<<unpacker_data.GetBHistograms() << "\n"    
	   << "                   C: "<<calibrator_data.GetBHistograms() << "\n"    
	   << "                   A: "<<analysis_data.GetBHistograms()  << "\n"<<std::endl;


  int Nloop=0;

  //loop over several input files;
  do{
    
    if(source_option==1){
      
      if(RLrun>0){ //if looping over a list of runs
	sprintf(RLname,"R%d_%d",RLrun,RLfirst+Nloop);    //*********Change here for AIDA/R run prefixes******************//
	FileNameData= DataDir + std::string(RLname);
      }
      else FileNameData= DataDir+RunName;

      std::cout << "\nSorting file "<<Nloop+1<<": "<< FileNameData << std::endl;
    }

    midas_data.InitDataSource(source_option, id, FileNameData, b_Xfer);
    if(b_verbose) midas_data.SetBDebug(true);
    else midas_data.SetBDebug(false);
    
    if(midas_data.GetBSourceOpen()){
      
      for(;;){
	
	/* Dealing with keyboard input
	//to pause loop and take a moment to reflect on the beauty of our planet and histograms 
	if(kbhit()==1){
	  int if_continue;
	  std::cout<<"online monitor is paused;"<<std::endl;
	  std::cout<<"press 0 <Enter> to stop, 1 <Enter> to restart, 2 <Enter> to clear histograms;"<<std::endl;
	  std::cin >> if_continue;
	  
	  char *choose;
	  if(if_continue==0){
	    choose=(char*)"stop";
	  }
	  else if(if_continue==1){
	    choose=(char*)"restart";
	  }
	  if(if_continue==2){
	    choose=(char*)"clear";
	  }
	  
	  std::cout<<"You chose "<<choose <<std::endl;
	  
	  if(if_continue == 0){
	    char my_ch;
	    std::cout<<"  0: exit AIDA sort (y/n)?"<<std::endl;
	    std::cin >> my_ch;
	    if(my_ch == 'y' || my_ch == 'Y'){
	      break;
	    }
	    else if_continue=1;
	  }
	  if(if_continue==1){
	    unpacker_data.UpdateHistograms();
	    calibrator_data.UpdateHistograms();
	    analysis_data.UpdateHistograms();
	    std::cout<<"monitor restarted, press <Space> to pause;"<<std::endl;
	  }
	  if(if_continue==2){
	    unpacker_data.ResetHistograms();
	    calibrator_data.ResetHistograms();
	    analysis_data.ResetHistograms();
	    std::cout<<"monitor reset, press <Space> to pause;"<<std::endl;
	  }
	}//end of statement controling keyboard halt of excecution: if(kbhit()==1)
	*/
	
	midas_data.Process();
	loops++;	

	if(midas_data.GetBPushData()){
	  unpacker_data.Process(midas_data);
	  
	  if(unpacker_data.GetBPushData()){
	    calibrator_data.Process(unpacker_data);
	    
	    if(calibrator_data.GetBPushData()){
	      analysis_data.Process(midas_data, calibrator_data);
	    }
	    
	    n_entry++;
	    
	    if( (n_entry%n_update)==0 ){
	      std::cout << n_entry << " entries unpacked..."<<std::endl;
	      unpacker_data.UpdateHistograms();
	      calibrator_data.UpdateHistograms();
	      //analysis_data.UpdateHistograms();
	      
	    }
	  }
	}
	
	// if(loops > 10000) {std::cout << "10000 entries processed. " << std::endl; break;}

	//how to break from online?
	if(midas_data.GetBEndOfData()){
	  std::cout << "\n\n ---- Reached end of input file --- " << std::endl;
	  break;
	}
      }
    } // If SourceData.IsOpen()
    
    midas_data.Close();  
    
    Nloop++;
  } while(Nloop<RLnum); //if list of runs, loop over RLnum data files.
  
  if(b_root_file){  //If write option selected - close analyses and write output files.
    
    std::cout << "\n   + SAVING HISTOGRAMS AND ROOT TREES TO FILE: " <<FileNameRoot<<std::endl;
    
    unpacker_data.Close();
    calibrator_data.Close();
    analysis_data.Close();
    
    std::cout << "Closing analysis sections. " << std::endl;

    fMain->Close();

    std::cout << "Closed output file. " << std::endl;
  }

  std::cout << "All done :)" << std::endl;
  std::cout << std::endl << "AIDA anlysis finished." << std::endl 
	    << "\t Elapsed real time: " << getRealTime()-start_t << "s" << std::endl
	    << "\t Elapsed CPU time: " << getCPUtime()-start_clk << "s" << std::endl;
  
  /* Close system interaction
  rootapp->Terminate();
  rootapp->Run();
  std::cout << "\n ++++ And this happens after exiting from rootapp ++++"<<std::endl;
  */

  return 0;
}
