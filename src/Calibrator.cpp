#include "Calibrator.h"
//#include "Common.h"
//#include "DutClass.h"


#include <bitset>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <string>
#include <ctime>

const unsigned char order[64]={62, 63, 59, 60, 61, 56, 57, 58, 52, 53, 54, 55, 49, 50, 51, 45,
			       46, 47, 48, 42, 43, 44, 38, 39, 40, 41, 35, 36, 37, 31, 32, 33,
			       34, 28, 29, 30, 24, 25, 26, 27, 21, 22, 23, 17, 18, 19, 20, 14,
			       15, 16, 10, 11, 12,  7,  3,  0,  8,  4,  1,  9,  5,  2, 13,  6};

/*! \fn Close()
 * Close the root file opened.
 *
 *
 */
void Calibrator::Close(){
  //close root file if opened
  Write();
  
  /*std::cout << std::endl << "Calibrator processed " << n_calib << " data items." 
	    <<"\t Elapsed time: " << common::getRealTime()-start_t << "s \t CPU time:" << common::getCPUtime()-start_clk << "s" 
	    <<"\t Data items per (real time) second: " << (double)n_calib/(common::getRealTime()-start_t) << std::endl;
  */
}

/*! \fn Process(Unpacker& my_unp_data)
 *
 *
 *
 */
void Calibrator::Process(Unpacker & my_unp_data){

  ResetData(); //set values of unp_data to defaults (zero)

  //channel will be set to default values for some bits of data where they don't apply (SYNC pulse)
  SetModule( my_unp_data.GetFee64Id() );
  SetChannel( my_unp_data.GetChId() );

  if( IsValidChannel(GetModule(), GetChannel() )){ //should already be checked if valid in Unpacker?
    
    SetBValidCh( true ); //not used! can be commented out ...
    SetBEnabled( b_ch_enabled[GetModule()-1][GetChannel()] );
    
  }
  else return; //now ensure it will be valid, no need to check again in each step....?

  SetDataType( my_unp_data.GetDataType() );
  //SetCorrFlag(my_unp_data.GetCorrFlag());
  //SetSyncFlag(my_unp_data.GetSyncFlag());

  //------------------------------------------------------------------
  //   if INFO DATA: should only be DISC. and CORR SCALARS
  //-----------------------------------------------------------------
  if(GetDataType()==2){ //Information word

    if(my_unp_data.GetInfoCode()==6){ //FAST DISCRIMINATOR hit

      SetTmStpDisc(my_unp_data.GetTmStp()); //set function will check validity of channel

    }

    else if(my_unp_data.GetInfoCode()==8){   //If correlation scalar

  //    if(  GetModule()==13 ||GetModule()==14 ||GetModule()==16  ){

	int64_t offset=  4*my_unp_data.GetCorrScaler() - my_unp_data.GetTmStp(); //AIDA= 100MHz Corr scaler= 25MHz

  if(tm_stp_corr_offset != offset){
    std::cout << "Old offset = " << tm_stp_corr_offset << " New offset = " << offset <<std::endl;
    SetTmStpOffset(offset);

  }

	SetBCorrStatus(true);

	if(b_debug){
	  std::cout << " OFFFFFFFFFFFFFFFFFFFFF: " << offset<<std::endl;

	  std::cout << "MODULE:"<<int(my_unp_data.GetFee64Id())<<"::  TS(aida), CORR: " << my_unp_data.GetTmStp() << ",  "<< my_unp_data.GetCorrScaler(); //<<",  "<< my_unp_data.GetCorrScaler()<<std::endl;
	  printf("  INFO 0x%lX \n", my_unp_data.GetInfoField());
	}

      //}

      SetInfoCode(my_unp_data.GetInfoCode());
      //SetBPushData(true);
      SetBFillTree(true);
    }

  } //End DataType == 2

  //------------------------------------------------------------------
  //   if ADC DATA: 
  //-----------------------------------------------------------------
  if(GetDataType()==3){

    SetADCdata(my_unp_data.GetADCdata());
    SetADCrange(my_unp_data.GetADCrange()) ;

    if(SetGeometry()){ //calculates DSSD, side, strip, ...

      CalibrateADC(); //energy... range...

      //time discriminator set after calculation of time-stamp
      if( GetBEnabled() && (double)GetADCenergy() > 150) SetBPushData(true); //std::cout << "Calibrator ADC at t:   " << my_unp_data.GetTmStp() << std::endl;}//if not a dissabled chanel (e.g. noisy ch)

      SetBFillTree(true);
    }
  }

  //if we're interested in this event
  if(GetBPushData() || GetBFillTree() ){

    //timing only if it makes sense...
    SetTimeAIDA(my_unp_data.GetTmStp());

    if(GetDataType()==3) SetDiscFlag(SetTimeDisc()); // check if thre is a fast disc. value for this ADC hit

    SetCorrFlag();
    SetTimeExternal();

    if(b_debug){
      if(GetDataType()==2 && (my_unp_data.GetInfoCode()==8 || my_unp_data.GetInfoCode()==4 )){
	if(GetModule()==13 ||GetModule()==14 ||GetModule()==16 ){
	  std::cout << " TIME AIDA: "<< GetTimeAIDA() << "  TIME EXT: "<< GetTimeExternal() 
		    <<  "   offset: "<<tm_stp_corr_offset <<std::endl;
	}
      }
    }
  }

  n_calib++;

  if( GetBRootTree() && GetBFillTree()){
    //save this entry to TTree
    out_root_tree->Fill();
  }

  if(GetBHistograms()) {   //If histogram output set to true fill histograms
    //Do we need any of the time-stamp variables calculated only for GetBPushData() or GetBFillTree()?
    FillHistograms(); 
  } 
}

/*********************************
 * Calibrator::FillHistograms()
 * Fill histograms
 *********************************/
void Calibrator::FillHistograms(){
  //fee 1 to N_FEE64:
  if(GetModule()>0 && GetModule()<=common::N_FEE64){

    if(GetADCrange()==0){ //Decay energy range
      hCalElCh[GetModule()-1]->Fill(GetChannel(),GetADCdata());
      hCalStripCount[GetModule()-1]->Fill(GetChannel());
    }
    else if(GetADCrange()==1){ //Implant energy range
      hCalEhCh[GetModule()-1]->Fill(GetChannel(),GetADCdata());
    }

    hist_fill_count++;
    if( ( hist_fill_count%n_update)==0 ) UpdateHistograms();
  }
}

/***********************************
 * Calibrator:: UpdateHistograms() *
 *   Update histograms if active   *
 ***********************************/
void Calibrator::UpdateHistograms(){
  if(GetBHistograms()){
    
    for(int i=0;i<common::N_FEE64;i++){
      cCal1->cd(i+1); hCalEhCh[i]->Draw("colz");
      cCal2->cd(i+1); hCalElCh[i]->Draw("colz");
    }
    cCal1->Update();
    cCal2->Update();
    if(GetBDebug()){
      std::cout << "   Calibrator::UpdateHistograms() at cnt " << hist_fill_count << std::endl;
    }
  }
}

/********************************************************
 *            Calibrator::ResetHistograms()             * 
 * Reset Histograms to empty (not currently implemented *
 ********************************************************/
void Calibrator::ResetHistograms(){
  if(GetBHistograms()){
    std::cout << ".... Calibrator:: and here we reset histograms; not implemented yet"<<std::endl;
  }
}

/**************************************************************************
 *                      Calibrator::IsValidChannel()                      *
 * Check if a channel is valid, i.e. between known FEE64 and strip limits *
 **************************************************************************/
bool Calibrator::IsValidChannel(int module, int channel){
  if(module > 0 && module <= common::N_FEE64){
    if(channel >=0 && channel < common::N_CHANNEL){
      return true; //common:n_fee64?
    }
  }
  return false;
}

/********************************************************************************
 *                         Calibrator::InitCalibrator()                         *
 * Initialise calibration step - load params, enable outputs, create histograms *
 ********************************************************************************/
void Calibrator::InitCalibrator(int opt, char *file_name){

  //start_t = common::getRealTime();
  //start_clk = common::getCPUtime();

  // if first of 'opt' is 1, enable histogramming
  if( (opt & 0x01) == 1){

    hist_fill_count = 0;

    //Initialise canvases and histograms
    SetBHistograms(true);
    
    if(common::N_FEE64%2==0) {
    cCal1= new TCanvas("cCal1","cCal1",50,50,900,900); cCal1->Divide(common::N_FEE64/2,2);
    cCal2= new TCanvas("cCal2","cCal2",90,90,900,900); cCal2->Divide(common::N_FEE64/2,2);
    }
    else if (common::N_FEE64%2==1) {
      cCal1= new TCanvas("cCal1","cCal1",50,50,900,900); cCal1->Divide((common::N_FEE64/2)+1,2);
      cCal2= new TCanvas("cCal2","cCal2",90,90,900,900); cCal2->Divide((common::N_FEE64/2)+1,2);
    }

    char hname[256];
    char htitle[256];
    for(int i=0;i<common::N_FEE64;i++){
      sprintf(hname,"hCalEhCh_FEE%i",i+1);
      sprintf(htitle,"NNAIDA%i [20GeV range];ch;ADC data",i+1);
      hCalEhCh[i]= new TH2I(hname, htitle, 64, 0, 64, 512, 0, 65536);

      cCal1->cd(i+1); hCalEhCh[i]->Draw("colz");  //here, or specific Draw function?

      sprintf(hname,"hCalElCh%i",i+1);
      sprintf(htitle,"NNAIDA%i [20MeV range];ch;ADC data",i+1);
      hCalElCh[i]= new TH2I(hname, htitle, 64, 0, 64, 512, 0, 65536);

      cCal2->cd(i+1); hCalElCh[i]->Draw("colz");  //here, or specific Draw function?
      
      sprintf(hname,"hCalStrips_FEE%i",i+1);
      sprintf(htitle, "Total number of hits in each strip - FEE#%i; strip; counts",i+1);
      hCalStripCount[i] = new TH1I(hname, htitle, 64, 0, 64);
      
    }
  
    if(b_debug) std::cout << "db    Calibrator.cpp: Initialized histograms and canvas for this step"<<std::endl;
  }

  //if second bit of 'opt' is 1, enable TTree for output
  if( ( (opt & 0x02) >> 1)  == 1){

    std::cout << " ***     Calibrator::InitCalibrator(): initializing TTree" << std::endl;
    //Initialize TTree
    out_root_tree = new TTree("AIDA_calibrate","AIDA_calibrate");
    out_root_tree->Branch("entry_calibrate",&cal_data,"time_aida/D:time_disc/D:time_external/D:adc_energy/D:adc_data/I:dssd/I:strip/I:adc_range/b:side/b:module/b:channel/b:data_type/b:info_code/b:corr_flag/O:disc_flag/O");
                                
    out_root_tree->SetAutoSave(-500000000);
    out_root_tree->SetAutoFlush(-50000000);

    std::cout << "  !!!!! Calibrator:: TTrree->GetAutoSave(): "<< out_root_tree->GetAutoSave() <<" !!! "<<std::endl;

    SetBRootTree(true);
  }

  LoadParameters(file_name);

  ResetData(); //sets values of unp_data structure to default values

}

/******************************************
 *     Calibrator::LoadParameters()       *
 * Load system parameters from input file *
 ******************************************/
void Calibrator::LoadParameters(char * file_name){

  std::ifstream param_file;

  param_file.open(file_name,std::ios::in); //need to specify as text?


  const int MAX_LINES= 40000; //4176 number of parameters used so far...
  int line_count=0;
  int par_count=0;

  std::string my_line;    // Line read in from file
  std::string par_name;   // Name of parameter
  int module, channel;    // Module and channel numbers
  double data;            // Data value to be set as parameter
  double x, y, z;         // Variables read in from file

  if(param_file.is_open()){  //If file open

    for(;;){

      if(std::getline(param_file, my_line)){

	if(GetBDebug()) std::cout << " * LoadParam() - line: " << my_line << std::endl;

	//skip empty lines and comments
	if(!my_line.empty() && '\n' != my_line[0] && '#' != my_line[0]){
	
	  std::istringstream iss(my_line);  
	  y = -1; //?? o a valid number
	  x = -1;
	  z= -1;

	  iss >> par_name >> x >> y >> z;

	  if(GetBDebug()) std::cout << " * LoadParameter() - values read: " << par_name << " "<< x << y << z << std::endl;

	  /******************
	  if(par_name=="b_mod_enabled"){
	    module= x; b_data= int(y);
	    if(IsValidChannel(module, 0)){
	      par_count++;
	      b_mod_enabled[module]= data;
	    }
	  }
	  **********************/
	  if(par_name=="map_dssd"){
	    module= x; data= y;
	    if(IsValidChannel(module, 0)){
	      par_count++;
	      map_dssd[module-1]= int(data);
	    }
	    else if(!IsValidChannel(module, 0)) {
	      std::cout << "Invalid channel for map_dssd: mod " << module << " does not exist." << std::endl;
	      std::cout << "Module enabled by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="map_side"){
	    module= x; data= y;
	    if(IsValidChannel(module, 0)){
	      par_count++;
	      map_side[module-1]= int(data);
	    }
	    else if(!IsValidChannel(module, 0)) {
	      std::cout << "Invalid channel for map_side: mod " << module << " does not exist." << std::endl;
	      std::cout << "Module enabled by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="map_strip"){
	    module= x; data= y;
	    if(IsValidChannel(module, 0)){
	      par_count++;
	      map_strip[module-1]= int(data);
	    }
	    else if(!IsValidChannel(module, 0)) {
	      std::cout << "Invalid channel for map_strip: mod " << module << " does not exist." << std::endl;
	      std::cout << "Module enabled by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="adc_polarity"){
	    module= x; data= y;
	    if(IsValidChannel(module, 0)){
	      par_count++;
	      adc_polarity[module-1]= int(data);
	    }
	    else if(!IsValidChannel(module, 0)) {
	      std::cout << "Invalid channel for adc_polarity: mod " << module << " does not exist." << std::endl;
	      //std::cout << "ADC polarity -1 by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="adc_gain_highE"){
	    module= x; data= y;
	    if(IsValidChannel(module, 0)){
	      par_count++;
	      adc_gain_highE[module-1]= data;
	    }
	    else if(!IsValidChannel(module, 0)) {
	      std::cout << "Invalid channel for adc_high_gainE: mod " << module << " does not exist." << std::endl;
	      std::cout << "adc_gain_highE set to default!!" << std::endl;
	    }
	  }
	  else if(par_name=="b_ch_enabled"){
	    module= x; channel= y;
	    data= z;
	    if(IsValidChannel(module, channel)){
	      par_count++;
	      b_ch_enabled[module-1][channel]= int( data );
	    }
	    else if(!IsValidChannel(module, channel)) {
	      std::cout << "Invalid channel for b_ch_enabled: mod " << module << " ch " << channel << " does not exist." << std::endl;
	      std::cout << "b_ch_enabled set to true by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="adc_offset"){
	    module= x; channel= y;
	    data= z;
	    if(IsValidChannel(module, channel)){
	      par_count++;
	      adc_offset[module-1][channel]= data;
	    }
	    else if(!IsValidChannel(module, channel)) {
	      std::cout << "Invalid channel for adc_offset: mod " << module << " ch " << channel << " does not exist." << std::endl;
	      std::cout << "adc_offset set to 0 by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="adc_thresh"){
	    module= x; channel= y;
	    data= z;
	    if(IsValidChannel(module, channel)){
	      par_count++;
	      adc_thresh[module-1][channel]= data;
	    }
	    else if(!IsValidChannel(module, channel)) {
	      std::cout << "Invalid channel for adc_thresh: mod " << module << " ch " << channel << " does not exist." << std::endl;
	      std::cout << "adc_thresh set to 400 by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="adc_gain"){
	    module= x; channel= y;
	    data= z;
	    if(IsValidChannel(module, channel)){
	      par_count++;
	      adc_gain[module-1][channel]= data;
	    }
	    else if(!IsValidChannel(module, channel)) {
	      std::cout << "Invalid channel for adc_gain: mod " << module << " ch " << channel << " does not exist." << std::endl;
	      std::cout << "adc_gain set to 0.7keV/ch by default!!" << std::endl;
	    }
	  }
	  else if(par_name=="disc_time_window"){
	    data= x;
	    disc_time_window= x;
	    par_count++;
	  }
	  else if(par_name=="aida_time_calib"){
	    data= x;
	    aida_time_calib= x;
	    par_count++;
	  }
	  
	  else{
	    if(b_debug) std::cout << " Calibrator.cpp:: parameter name not recognized: " << par_name << "." <<std::endl;
	  }
	  
	} //End if my_line.empty()
	
	
      } //End if getline

      else{
	std::cout << " Calibrator.cpp:: LoadParameters() reached end of parameter file. Loaded new value for "<< par_count << " parameters." << std::endl;
	break;
      }

      //check forever loop not running amok
      line_count++;
      if(line_count>MAX_LINES){

	std::cout << " Calibrator.cpp:: count of lines read by LoadParameters() exceeds maximum allowed value (" << MAX_LINES << "). Parameters loaded= " << par_count<<std::endl;
	std::cout << "                  Break out of loop and continue... *****************"  << std::endl << std::endl;
	break;
      }
    }

  }
  else{
    std::cout << " Calibrator.cpp:: Error opening parameter file: " << file_name << " *****************"  << std::endl << std::endl;
  }
  return;

}

/*******************************************************************
 *                    Calibrator::ResetData()                      *
 * Reset data struct to default values before moving to next event *
 *******************************************************************/
void Calibrator::ResetData(){
  
  b_push_data = false;
  b_fill_tree = false;

  //reset this?
  b_enabled  = false;
  b_valid_ch = false;

  cal_data.time_aida     = 0;
  cal_data.time_disc     = 0;
  cal_data.time_external = 0;

  cal_data.adc_energy    = -999;
  cal_data.dssd          = -999;
  cal_data.strip         = -999;

  cal_data.adc_range     = 255;
  cal_data.side          = 255;
  cal_data.module        = 255;
  cal_data.channel       = 255;
  cal_data.data_type     = 255;
  cal_data.info_code     = 255;

  //cal_data.sync_flag= false;
  cal_data.corr_flag = false; 
  cal_data.disc_flag = false; 
}

/************************************
 *        Calibrator::Write()       *
 * Write output tree and histograms *
 ************************************/
void Calibrator::Write(){
  if(GetBRootTree()){

    out_root_tree->Write(0,TObject::kOverwrite);
    std::cout << "\n *** writing Calibrator TTree .... "<< std::endl;
    out_root_tree->Print();
  }
  if(GetBHistograms()){

    for(int i=0;i<common::N_FEE64;i++){
      hCalEhCh[i]->Write();
      hCalElCh[i]->Write();
      hCalStripCount[i]->Write();
    }
  }
}

/***********************************************************************************
 *                           Calibrator::SetGeometry()                             *
 * Set DSSD#,n/p sdie and strip number of data event accprding to input parameters *
 ***********************************************************************************/
bool Calibrator::SetGeometry(){

  //if(GetBValidCh()){

  cal_data.dssd = map_dssd[GetModule()-1];
  if(cal_data.dssd < 1) std::cout << "&&&&& CALIB DSSD# < 1 &&&&&" << std::endl;
  
  // higher number channels for strips at center of detector (FEE for each half)
  if(map_strip[GetModule()-1] == 1){
    cal_data.strip = OrderChannel( GetChannel() );
  } else if(map_strip[GetModule()-1]== 2){
    cal_data.strip = 127 -  OrderChannel( GetChannel() ); 
  } else if(map_strip[GetModule()-1] != 1 && map_strip[GetModule()-1] != 2) {
    cal_data.strip = -1;
    return false;
  }
  
  cal_data.side = map_side[GetModule()-1];
  
  return true;
  //}
}

/****************************************************************
 *                  Calibrator::CalibrateADC()                  *
 * Calibrate ADC value to account for ADC offset and zero point *
 ****************************************************************/
void Calibrator::CalibrateADC(){

  if(GetADCrange() == 0){  //LowE decay range
    SetADCenergy( ( GetADCdata() - common::ADC_ZERO - adc_offset[GetModule()-1][GetChannel()] )  //ADC value adjusted for zero-point and offset
		  * adc_polarity[GetModule()-1] * adc_gain[GetModule()-1][GetChannel()] );       //Multiplied by side parity and gain.
  }

  else if(GetADCrange() == 1){  //HighE implant range
    SetADCenergy( (GetADCdata() - common::ADC_ZERO)
		  *  adc_polarity[GetModule()-1] * adc_gain_highE[GetModule()-1] );
  }
  
  else SetADCenergy(-88888);
}

/*******************************************************************
 *                   Calibrator::OrderChannel()                    *
 * Return correct channel after correction for adapter PCB routing *
 *******************************************************************/
unsigned char Calibrator::OrderChannel( unsigned char ch ){
  //assumes we have already checked channel is within valid range!... likely in enabled channels check.
  return order[ch];
}

/*************************************************
 *           Calibrator::SetBDebug()             *
 * Set debug bool flag to receive debug comments *
 *************************************************/
void Calibrator::SetBDebug(bool flag){
  b_debug= flag;
}

/****************************************************
 *          Calibrator::SetBHistograms()            *
 * Set histogram flag to print histograms if wanted *
 ****************************************************/
void Calibrator::SetBHistograms(bool flag){
  b_histograms= flag;
}

/*************************************************
 *           Calibrator::SetBPushData()          *
 * Set flag to push data...?                     *
 *************************************************/
void Calibrator::SetBPushData(bool flag){
  b_push_data= flag;
}

/********************************************************************
 * Calibrator::SetBFillTree()
 * Set BFillTree flag is we are interested in writing event to tree *
 ********************************************************************/
void Calibrator::SetBFillTree(bool flag){
  b_fill_tree= flag;
}

void Calibrator::SetBRootTree(bool flag){
  if(flag){
    if( out_root_tree != NULL){
      b_root_tree= true;
      std::cout << " ***       Calibrator::SetBRootTree(): set to TRUE"<<std::endl;
    }
    else {
      std::cout << " ** WARNING ** Calibrator.cpp: attempted to set b_root_tree as valid when TTree is not initialized! *** " << std::endl;
      b_root_tree= false;
    }
  }
  else b_root_tree= false;
}

/**************************************
 * Calibrator::SetBEnabled()
 *
 **************************************/
void Calibrator::SetBEnabled(bool flag){
  b_enabled= flag;
}

void Calibrator::SetBValidCh(bool flag){
  b_valid_ch= flag;
}

void Calibrator::SetBCorrStatus(bool flag){
  b_corr_status= flag;
}

//void Calibrator::SetBCorrOffset(bool flag){
// b_corr_offset= flag;
//}


void Calibrator::SetTmStpDisc(unsigned long value){
  //  if(GetBValidCh()){
    tm_stp_disc[GetModule()-1][GetChannel()]= value;
    //  }

}

void Calibrator::SetTmStpOffset(int64_t value){ //diff AIDA->EXTERNAL tm-stp
  //double offset=  aida_time_calib*(my_unp_data.GetCorrScaler() - my_unp_data.GetTmStp());

  tm_stp_corr_offset= value; 
  //tm_stp_corr_offset= 1.*t_EXT - aida_time_calib*t_AIDA;
}

void Calibrator::SetTimeAIDA(double value){
  cal_data.time_aida= value;
}

bool Calibrator::SetTimeDisc(){
  /********************
    //check we have detected previous value for DISC for this channel
    if(tm_stp_disc[GetModule()][GetChannel()]>0){

      cal_data.time_disc= tm_stp_disc[GetModule()][GetChannel()];

      tm_stp_disc[GetModule()][GetChannel()]=0; //clear value of stored DISC timestamp for this channel

      if( (GetTimeAIDA() - cal_data.time_disc < disc_time_window) && (GetTimeAIDA()>=cal_data.time_disc) ){
	return true;
      }
      else return false; 
    }
  return false;
  ***********************/
  return false;
}


void Calibrator::SetTimeExternal(){

  //  if(GetCorrFlag() && GetSyncFlag()){  
  if(GetCorrFlag()){  

    cal_data.time_external = (cal_data.time_aida + tm_stp_corr_offset) /  aida_time_calib; 
    //    cal_data.time_external= aida_time_calib*(cal_data.time_aida + tm_stp_corr_offset); 
  }
  else cal_data.time_external= 0;
}

void Calibrator::SetADCenergy(double value){
  cal_data.adc_energy= value;
}

void Calibrator::SetADCdata(int value){
  cal_data.adc_data= value;
}
void Calibrator::SetDSSD(int value){
  cal_data.dssd= value;
}

void Calibrator::SetStrip(int value){
  cal_data.strip= value;
}

void Calibrator::SetADCrange(unsigned char value){
  cal_data.adc_range= value;
}

void Calibrator::SetSide(unsigned char value){
  cal_data.side= value;
}

void Calibrator::SetModule(unsigned char value){
  cal_data.module= value;
}

void Calibrator::SetChannel(unsigned char value){
  cal_data.channel= value;
}

void Calibrator::SetDataType(unsigned char value){
  cal_data.data_type= value;
}

void Calibrator::SetInfoCode(unsigned char value){
  cal_data.info_code= value;
}

//void Calibrator::SetSyncFlag(bool value){
//  cal_data.sync_flag= value;
//}

void Calibrator::SetCorrFlag(){

  cal_data.corr_flag= GetBCorrStatus() ;
}

void Calibrator::SetDiscFlag(bool value){
  cal_data.disc_flag= value;
}

// ********************************************
// ***************  Getters  ******************
// ********************************************
bool Calibrator::GetBDebug(){ return b_debug; }
bool Calibrator::GetBHistograms(){ return b_histograms; }
bool Calibrator::GetBPushData(){ return b_push_data; }
bool Calibrator::GetBFillTree(){ return b_fill_tree; }
bool Calibrator::GetBRootTree(){ return b_root_tree; }

bool Calibrator::GetBEnabled(){ return b_enabled; }
bool Calibrator::GetBValidCh(){ return b_valid_ch; }
bool Calibrator::GetBCorrStatus(){ return b_corr_status; }

unsigned long Calibrator::GetTmStpDisc(){ return tm_stp_disc[GetModule()-1][GetChannel()]; }

int64_t Calibrator::GetTmStpOffset(){ return tm_stp_corr_offset; }

//Getters... for cal_data structure
double Calibrator::GetTimeAIDA(){ return cal_data.time_aida; }
double Calibrator::GetTimeDisc(){ return cal_data.time_disc; }
double Calibrator::GetTimeExternal(){ return cal_data.time_external; }
double Calibrator::GetADCenergy(){ return cal_data.adc_energy; }

int Calibrator::GetADCdata(){ return cal_data.adc_data; }
int Calibrator::GetDSSD(){ return cal_data.dssd; }
int Calibrator::GetStrip(){ return cal_data.strip; }

unsigned char Calibrator::GetADCrange(){ return cal_data.adc_range; }
unsigned char Calibrator::GetSide(){ return cal_data.side; }
unsigned char Calibrator::GetModule(){ return cal_data.module; }
unsigned char Calibrator::GetChannel(){ return cal_data.channel; }
unsigned char Calibrator::GetDataType(){ return cal_data.data_type; }
unsigned char Calibrator::GetInfoCode(){ return cal_data.info_code; }

//bool Calibrator::GetSyncFlag(){ return cal_data.sync_flag; }
bool Calibrator::GetCorrFlag(){ return cal_data.corr_flag; }
bool Calibrator::GetDiscFlag(){ return cal_data.disc_flag; }


Calibrator::Calibrator(){

  out_root_tree = NULL;

  b_debug       = false;

  b_root_tree   = false;
  b_histograms  = false;

  b_push_data   = false;
  b_fill_tree   = false;

  for(int i=0;i<common::N_FEE64;i++){
    for(int j=0;j<common::N_CHANNEL;j++){
      tm_stp_disc[i][j]=0;    
    }
  }
  
  tm_stp_corr_offset=0;

  //...Parameters....
  //by default, all channels enabled
  for(int i=0;i<common::N_FEE64;i++){
    for(int j=0;j<common::N_CHANNEL;j++){
      b_ch_enabled[i][j]= true;
    }
  }

  //Assign default dssd/side/strip numbers
  for(int i=0;i<common::N_FEE64;i++){
    map_dssd[i] = i/4; // FEE# = 1-4 > 0, 5-8 > 1, 9-12 > 2 etc. This is incorrect assignment and should be changed with config parameters.

    if( (i%4)<2 ) map_side[i] = 1; // i = 0,1,4,5,8,9 ... => FEE# = 1, 2, 5, 6, 9, 10, etc... (front/p-side)
    else map_side[i] = 0;          // i = 2,3,6,7,10,11 ... => FEE# = 3, 4, 7, 8, 11, 12 etc... (back/n-side)

    map_strip[i] = i%2; // FEE# even = 0, odd = 1. This is incorrect assignment and should be changed with config parameters.
  }

  for(int i=0;i<common::N_FEE64;i++){
    for(int j=0;j<common::N_CHANNEL;j++){
      adc_offset[i][j] = 0;
      adc_gain[i][j] = 0.7; // keV/ch
    }
    if( (i%4)<2 ) adc_polarity[i] = 1;
    else adc_polarity[i] = -1;
    adc_gain_highE[i] = 0.7; // MeV/ch
  }

  disc_time_window = 50000;
  aida_time_calib = 2; //     ???????AIDA: 25MHz, BRIKEN: 50Hz

  b_enabled  = false;
  b_valid_ch = false;

  b_corr_status = false;

  cal_data.time_aida     = 0;
  cal_data.time_disc     = 0;
  cal_data.time_external = 0;
  cal_data.adc_energy    = 0;

  cal_data.adc_data      = -999;
  cal_data.dssd          = -999;
  cal_data.strip         = -999;

  cal_data.adc_range     = 255;
  cal_data.side          = 255;
  cal_data.module        = 255;
  cal_data.channel       = 255;
  cal_data.data_type     = 255;
  cal_data.info_code     = 255;

  //cal_data.sync_flag= false;
  cal_data.corr_flag = false; 
  cal_data.disc_flag = false; 

}
