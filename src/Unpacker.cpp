#include "Unpacker.h"

#include <bitset>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <string>

void Unpacker::Close(){

  if(GetBRootTree()){
    out_root_tree->Write(0,TObject::kOverwrite);
    std::cout << "\n *** writting Unpacker TTree .... "<< std::endl;
    out_root_tree->Print();
  }

  if(GetBHistograms())  WriteHistograms();

}

//the heart of the matter
void Unpacker::Process(DataSource & my_src_data){

  ResetData(); //set values of unp_data to defaults (zero)

  unsigned char my_char;
  unsigned long my_long;

  unsigned int word_0= my_src_data.GetWord0();
  unsigned int word_1= my_src_data.GetWord1();

  //Top 2 bits of first word - bits 31:30
  SetDataType( (word_0 >> 30) & 0x3 );       //Determines structure of following information.

  //ADC data
  if(GetDataType() == 3){
    unsigned int chIdent = (word_0 >> 16) & 0x0FFF;  //Word 0, bits 16:27
    SetFee64Id( (chIdent >> 6) & 0x003F );           //Top 6 bits of chIdent (22:27) are FEE64 number
    SetChId( chIdent & 0x003F );                     //Bottom 6 bits of chIdent(16:21) are channel number
    SetADCrange( (word_0 >> 28) & 0x0001 );          //Word 0, bit 28 - Veto bit used as ADC range
    SetTmStpLsb( word_1 & 0x0FFFFFFF );              //Word 1, bits 0:27 - timestamp LSB

    SetADCdata( word_0 & 0xFFFF );                   //Bits 0:15 of first word - ADC data

    if(IsValidFee64Id() && b_mod_enabled[GetFee64Id()-1]) SetBPushData(true);  //Push data to Calibrator.cpp is valid mod and ch

    SetBFillTree(true);
  
    ++counter;  
    #ifdef DEBUG_UNP
    if( int((chIdent >> 6) & 0x003F) < 1)  std::cout << "*********************DSSD# < 1 *************************" << std::endl;
    //std::cout << "%%%%%ADC WORD%%%%%  t_lsb: " << (int)(word_1 & 0x0FFFFFFF) << "\t t_msb: " << tm_stp_msb << "\t t: " << unp_data.tm_stp << std::endl;
    //if(counter < 20 ) std::cout << unp_data.tm_stp << std::endl;
    #endif

  } 
  
  //Information data
  else if(GetDataType() == 2){
    SetInfoField( word_0 & 0x000FFFFF );   //Word 0, bits 0:19 - information field
    SetInfoCode( (word_0 >> 20) & 0x000F); //Word 0, bits 20:23 - information code
    SetFee64Id( (word_0 >> 24) & 0x003F);  //Word 0, bits 24:29 - FEE64 module number

    SetTmStpLsb( word_1 & 0x0FFFFFFF );    //Word 1, bits 0:27 - timestamp LSB

    // 2== PAUSE
    if(GetInfoCode()==2){
      tm_stp_msb = GetInfoField();
      //if(!GetBSyncStatus()) SetBSyncStatus(true); //data is synchronized!
    }

    // 3== RESUME
    else if(GetInfoCode()==3){
      tm_stp_msb = GetInfoField();
      if(!GetBSyncStatus()) SetBSyncStatus(true); //data is synchronized!
    }

    // 4==SYNC100 pulse
    else if(GetInfoCode()==4){
      //#ifdef DEBUG_UNP
      //std::cout << "*******SYNC100 PULSE********* t_msb: " << GetInfoField() << std::endl;
      //#endif
      tm_stp_msb = GetInfoField();
      if(!GetBSyncStatus()) SetBSyncStatus(true); //data is synchronized!

      //  SetBPushData(true);
    }

    // 6==DISC data
    else if(GetInfoCode()==6){
      // !!!!!!! BUG !!!!!!!
      // Discriminator format changed - this is no longer correct.
      // InfoField now contains hit pattern per ASIC - lower 16 bits is hit pattern, top 4 bits are ASIC ID (not sure if 0->3 or 1->4)
      // !!!!!!!!!!!!!!!!!!! 05/10/16
      //my_char = GetInfoField() & 0xFFFF;  //channel number is encoded in info-field
      //SetChId( GetInfoField() & 0xFFFF );   //Why is this the case? [CG 02/09/16]
  
      

      //SetBPushData(true);
    }

    // 8==Correlation Scaler
    else if(GetInfoCode()==8){
      my_char = ( GetInfoField() & 0x000F0000 ) >> 16; //index of corr scaler
      my_long = ( GetInfoField() & 0x0000FFFF ) ;      //bits with timestamp

      //assumes numbering FEE64 modules starts from 1 (and now goes 1->common::N_DSSD)
      if(IsValidFee64Id()) {
	if(my_char==0)      corr_scaler_data0[GetFee64Id()-1]= my_long;
	else if(my_char==1) corr_scaler_data1[GetFee64Id()-1]= my_long;
	else if(my_char==2){
	  unsigned long scaler= my_long << 32 | corr_scaler_data1[GetFee64Id()-1] << 16 | corr_scaler_data0[GetFee64Id()-1];

	  if(b_debug){
	    if(GetFee64Id()==13 ||GetFee64Id()==14 ||GetFee64Id()==16 ){  //Why these FEE64 modules? Arbitrary? [CG 06/09/16]
	      
	      std::cout << " CORR SCALER: "<< 1.*scaler <<std::endl;
	    }
	  }
	  
	  //updates state of b_corr_status, if we have already received SYNC100 pulses
	  //also updates value of corr_scaler_offset
	  SetCorrScaler(scaler);
	  //	  SetBCorrStatus(true);
	  SetBPushData(true);
	}
      }
    }
    
    n_unp++;
    SetBFillTree(true); //Fill TTree for all Info Data Codes
    // SetPushData(true); !!only push data for selected info types: now DISC and CORR SCALAR
  } //End of DataType==2

  //Sample trace: Sample Length
  else if(GetDataType() == 1){
    unsigned int chIdent = (word_0 >> 16) & 0x0FFF; //Word 0, bits 16:27 - channel identifier
    
    SetSampleLength( word_0 & 0xFFFF );             //Word 0, bits 0:15 - sample length
    SetChId( chIdent & 0x003F );                    //Word 0, bits 16:21 - channel number
    SetFee64Id( (chIdent >> 6) & 0x003F );          //Word 0, bits 22:27 - FEE64 number

    SetTmStpLsb( word_1 & 0x0FFFFFFF );             //Word 1, bits 0:27 - timestamp LSB

    //these are default values
    // SetBFillTree(false);
    //    SetBPushData(false);
  }

  //Sample trace: Waveform data
  else if(GetDataType() == 0){
    // 4 x 14 bits Waveform samples
    // NOT IMPLEMENTED PROPERLY, only first value saved now (need more variables in unpack_data_struct)

    //my_int=  (word_0 >> 16) & 0x00003FFF;    //sample n
    SetSampleLength( (word_0 >> 16) & 0x00003FFF );

    //these are default values
    //    SetBFillTree(false);
    //    SetBPushData(false); //not being processed for now.
  }

  //If 30:31 do not conform to any expected values, produce error message.
  else{
    //output error message!
    std::cerr << " **WARNING** Unpacker.cpp: DATA TYPE NOT RECOGNIZED! word 0: " << word_0 << ", data_type: "<<GetDataType()<<std::endl;
    SetBPushData(false);
  }

  //do following (settmstp, flags, etc..) only if BPushData
  SetTmStp(); //reconstruct full time stamp
  SetFlags();

  if(!GetSyncFlag()) SetBPushData(false); //do not forward data out of SYNC100 (need to improve test of SYNC status)

  // If an interesting data type, Fill TTree with this entry 
  if(GetBRootTree() && GetBFillTree()) out_root_tree->Fill();
    
  if(GetBHistograms()) FillHistograms();

  if(b_debug){
    if(my_src_data.GetItrData()<65){ //for first few entries in block output to screen
      printf("db      Unp:: DataType= %i, NNAIDA%i, tm-stp(lsb)= %lu, ",GetDataType(),GetFee64Id(),GetTmStpLsb());
      if(GetDataType()==3) printf("adc= %u\n",GetADCdata());
      else if(GetDataType()==2) printf("code= %i\n",GetInfoCode());
    }
    else if(GetFee64Id() == 3 && GetDataType() == 3 && GetChId()<2){
      printf("db      Unp:: NNAIDA3: ch= %i, adc= %u\n",GetChId(),GetADCdata());
    }
  }
}

void Unpacker::FillHistograms(){
  
  if(GetDataType()==3){ //ADC
    if(GetADCrange()==0){ //Low energy range
      hFEE64_ADClow->Fill(GetFee64Id());
      hCh_FEE64_ADClow->Fill(GetFee64Id(),GetChId());

    }
    else if(GetADCrange()==1){ //High energy range
      hFEE64_ADChigh->Fill(GetFee64Id());
      hCh_FEE64_ADChigh->Fill(GetFee64Id(),GetChId());
    }
  }

  else if(GetDataType()==2){ //InfoCode

    hFEE64_Info->Fill(GetFee64Id());

    hInfoCode_FEE64->Fill(GetFee64Id(),GetInfoCode());

    if(GetInfoCode()==2) hFEE64_PAUSE->Fill(GetFee64Id());
    else if(GetInfoCode()==3) hFEE64_RESUME->Fill(GetFee64Id());
    else if(GetInfoCode()==4) hFEE64_SYNC100->Fill(GetFee64Id());

    else if(GetInfoCode()==6) hCh_FEE64_DISC->Fill(GetFee64Id(), GetChId() );
  }

  else if(GetDataType()==1) hFEE64_Waveform->Fill(GetFee64Id());
}

bool Unpacker::IsValidFee64Id(){
  if(GetFee64Id() >0 && GetFee64Id() <= common::N_FEE64) return true; //common:n_fee64?
  else return false;
}

bool Unpacker::IsValidFee64Id(int mod_id){
  if(mod_id >0 && mod_id <= common::N_FEE64) return true; //common:n_fee64?
  else return false;
}

void Unpacker::InitUnpacker(int opt, char *file_name){

  //start_t = common::getRealTime();
  //start_clk = common::getCPUtime();

  // if first of 'opt' is 1, enable histogramming
  if( (opt & 0x01) == 1){

    std::cout << "     Unpacker::InitUnpacker(): defining histograms"<<std::endl;

    SetBHistograms(true);

    char hname[256];
    char htitle[256];

    cUnp1= new TCanvas("cUnp1","cUnpacker 1",20,20,900,900); cUnp1->Divide(4,3);

    sprintf(hname,"hFEE64_ADClow");
    sprintf(htitle,"FEE64 !ADClow;module ID");
    hFEE64_ADClow = new TH1I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1);

    sprintf(hname,"hFEE64_ADChigh");
    sprintf(htitle,"FEE64 !ADChigh;module ID");
    hFEE64_ADChigh = new TH1I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1);

    sprintf(hname,"hFEE64_Waveform");
    sprintf(htitle,"FEE64 !Waveform;module ID");
    hFEE64_Waveform = new TH1I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1);

    sprintf(hname,"hFEE64_Info");
    sprintf(htitle,"FEE64 !Info;module ID");
    hFEE64_Info = new TH1I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1);

    cUnp1->cd(1); hFEE64_ADClow->Draw();
    cUnp1->cd(2); hFEE64_ADChigh->Draw();
    cUnp1->cd(3); hFEE64_Waveform->Draw();
    cUnp1->cd(4); hFEE64_Info->Draw();

    sprintf(hname,"hInfoCode_FEE64");
    sprintf(htitle,"Info Code vs FEE64;module ID;info code");
    hInfoCode_FEE64 = new TH2I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1,16,0,16);
    cUnp1->cd(5); hInfoCode_FEE64->Draw("col");

    sprintf(hname,"hFEE64_PAUSE");
    sprintf(htitle,"FEE64 !PAUSE;module ID");
    hFEE64_PAUSE = new TH1I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1);
    cUnp1->cd(6); hFEE64_PAUSE->Draw("");

    sprintf(hname,"hFEE64_RESUME");
    sprintf(htitle,"FEE64 !RESUME;module ID");
    hFEE64_RESUME = new TH1I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1);
    //  hFEE64_RESUME->SetLineColor(2);
    cUnp1->cd(7); hFEE64_RESUME->Draw("");

    sprintf(hname,"hFEE64_SYNC100");
    sprintf(htitle,"FEE64 !SYNC100;module ID");
    hFEE64_SYNC100 = new TH1I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1);
    cUnp1->cd(8); hFEE64_SYNC100->Draw("");

    //  sprintf(hname,"hFEE64_Waveform");
    //  sprintf(htitle,"FEE64 !Waveform;module ID");
    //  hFEE64_Waveform = new TH1I(hname, htitle, common::N_FEE64,0,common::N_FEE64);
    //  cUnp1->cd(8); hFEE64_Waveform->Draw("");

    sprintf(hname,"hCh_FEE64_ADClow");
    sprintf(htitle,"channel vs FEE64 !ADClow;module ID;ch ID");
    hCh_FEE64_ADClow = new TH2I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1,64,0,64);
    cUnp1->cd(9); hCh_FEE64_ADClow->Draw("col");

    sprintf(hname,"hCh_FEE64_ADChigh");
    sprintf(htitle,"channel vs FEE64 !ADChigh;module ID;ch ID");
    hCh_FEE64_ADChigh = new TH2I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1,64,0,64);
    cUnp1->cd(10); hCh_FEE64_ADChigh->Draw("col");

    sprintf(hname,"hCh_FEE64_DISC");
    sprintf(htitle,"channel vs FEE64 !DISC;module ID;ch ID");
    hCh_FEE64_DISC = new TH2I(hname, htitle, common::N_FEE64,1,common::N_FEE64+1,64,0,64);
    cUnp1->cd(11); hCh_FEE64_DISC->Draw("col");


    if(b_debug) std::cout << "db    Unpacker.cpp: Initialized histograms and canvas for this step"<<std::endl;
  }

  //if second bit of 'opt' is 1, enable TTree for output
  //  if( (opt & 0x02) == 1){
  if( ( (opt & 0x02) >> 1)  == 1){

    std::cout << " ***     Unpacker::InitUnpacker(): initializing TTree" << std::endl;

    //Initialize TTree
    out_root_tree = new TTree("AIDA_unpack","AIDA_unpack");
    // 16/6: remove corr_flag
    out_root_tree->Branch("entry_unpack",&unp_data,"tm_stp/l:corr_scaler/l:tm_stp_lsb/l:info_field/l:adc_data/i:sample_length/i:data_type/b:fee64_id/b:ch_id/b:adc_range/b:info_code/b:sync_flag/O");

    out_root_tree->SetAutoSave(-500000000);
    out_root_tree->SetAutoFlush(-50000000);

    std::cout << "  !!!!! Unpacker:: TTrree->GetAutoSave(): "<< out_root_tree->GetAutoSave() <<" !!! "<<std::endl;

    SetBRootTree(true);

  }

  LoadParameters(file_name);

  ResetData(); //sets values of unp_data structure to default values

}

void Unpacker::UpdateHistograms(){

  if(GetBHistograms()){

    for(int i=0; i<12; i++){
      cUnp1->cd(i+1)->Modified(); 
    }
    cUnp1->Update();
  }
}

void Unpacker::ResetHistograms(){

  if(GetBHistograms()){
    hFEE64_ADClow->Reset();
    hFEE64_ADChigh->Reset();  
    hFEE64_Waveform->Reset();
    hFEE64_Info->Reset();

    hInfoCode_FEE64->Reset();

    hFEE64_PAUSE->Reset();
    hFEE64_RESUME->Reset();
    hFEE64_SYNC100->Reset();

    hCh_FEE64_ADClow->Reset();
    hCh_FEE64_ADChigh->Reset();
    hCh_FEE64_DISC->Reset();
  }

}

void Unpacker::ResetData(){

  b_push_data = false;
  b_fill_tree = false;

  unp_data.tm_stp      = 0;
  unp_data.corr_scaler = 0;
  unp_data.tm_stp_lsb  = 0;
  unp_data.info_field  = 0;

  unp_data.adc_data      = 0;
  unp_data.sample_length = 0;
  unp_data.data_type     = 0;
  unp_data.fee64_id      = 0;
  unp_data.ch_id         = 0;
  unp_data.adc_range     = 0;
  unp_data.info_code     = 0;

  unp_data.sync_flag = false;
  //  unp_data.corr_flag= false; 

}

void Unpacker::WriteHistograms(){

  //write histograms to TFile
  hFEE64_ADClow->Write();
  hFEE64_ADChigh->Write();  
  hFEE64_Waveform->Write();
  hFEE64_Info->Write();

  hInfoCode_FEE64->Write();

  hFEE64_PAUSE->Write();
  hFEE64_RESUME->Write();
  hFEE64_SYNC100->Write();

  hCh_FEE64_ADClow->Write();
  hCh_FEE64_ADChigh->Write();
  hCh_FEE64_DISC->Write();

  //and also TCanvas
  cUnp1->Write();

  std::cout << "      Unpacker::WriteHistograms(): done writing to file"<<std::endl;

}

void Unpacker::LoadParameters(char * file_name){

  std::ifstream param_file;

  param_file.open(file_name,std::ios::in); //need to specify as text?


  const int MAX_LINES= 40000; //4176 number of parameters used so far...
  int line_count=0;
  int par_count=0;

  std::string my_line;
  std::string par_name;
  int module;
  double data;
  double x, y, z;
  //  bool b_data;

  if(param_file.is_open()){

    for(;;){

      if(std::getline(param_file, my_line)){

	if(b_debug) std::cout << " * LoadParam() - line: " << my_line << std::endl;

	//skip empty lines and comments
	if(!my_line.empty() && '\n' != my_line[0] && '#' != my_line[0]){
	
	  std::istringstream iss(my_line);  
	  y = -1; //?? o a valid number
	  x = -1;
	  z= -1;
	  // par_name[0]= '\0';

	  iss >> par_name >> x >> y >> z;


	  if(b_debug) std::cout << " * LoadParameter() - values read: " << par_name << " "<< x << y << z << std::endl;
	  

	  if(par_name=="b_mod_enabled"){
	    module= x; data= y;
	    if(IsValidFee64Id(module)){
	      par_count++;
	      if(y==1){
		b_mod_enabled[module-1] = true;
	      }
	      else if(y==0) {
		b_mod_enabled[module-1] = false;
	      }
	    }
	  }
	  else{
	    if(b_debug) std::cout << " Unpacker.cpp:: parameter name not recognized: " << par_name << "." <<std::endl;
	  }
	
	}
	

      }
      else{
	std::cout << " Unpacker.cpp:: LoadParameters() reached end of parameter file. Loaded new value for "<< par_count << " parameters." << std::endl;
	break;
      }
      //check forever loop not running amok
      line_count++;
      if(line_count>MAX_LINES){

	std::cout << " Unpacker.cpp:: count of lines read by LoadParameters() exceeds maximum allowed value (" << MAX_LINES << "). Parameters loaded= " << par_count<<std::endl;
	std::cout << "                  Break out of loop and continue... *****************"  << std::endl << std::endl;
	break;
      }
    }

  }
  else{
    std::cout << " Unpacker.cpp:: Error opening parameter file: " << file_name << " *****************"  << std::endl << std::endl;
  }
  return;

}

void Unpacker::SetBSyncStatus(bool flag){
  b_sync_status= flag;
}

//void Unpacker::SetBCorrStatus(bool flag){
//  b_corr_status= flag;
//}

void Unpacker::SetBDebug(bool flag){
  b_debug= flag;
}

void Unpacker::SetBHistograms(bool flag){
  b_histograms= flag;
}

void Unpacker::SetBPushData(bool flag){
  b_push_data= flag;
}

void Unpacker::SetBFillTree(bool flag){
  b_fill_tree= flag;
}

void Unpacker::SetBRootTree(bool flag){
  if(flag){
    if(out_root_tree){
      b_root_tree= true;
    }
    else {
      std::cout << " ** WARNING ** Unpacker.cpp: attempted to set b_root_tree as valid when TTree is not initialized! *** " << std::endl;
      b_root_tree= false;
    }
  }
  else b_root_tree= false;
}

void Unpacker::SetASICid(unsigned char value) {
  asic_id = value;
}

void Unpacker::SetFastDiscHits(unsigned char value) {
  fast_disc_hits = value;
}

void Unpacker::SetTmStp(){
  if(GetBSyncStatus()){
    unp_data.tm_stp = ( tm_stp_msb << 28 ) | unp_data.tm_stp_lsb ; //or do I use a parameter in function call for tm_stp_lsb?
   }
  else unp_data.tm_stp = 0;
}

void Unpacker::SetFlags(){
  unp_data.sync_flag= GetBSyncStatus(); // b_sync_status;
  //  unp_data.corr_flag= b_corr_status;
}

void Unpacker::SetCorrScaler(unsigned long value){
  unp_data.corr_scaler= value;
  //  unp_data.corr_flag= true; (we'll have at some point a proper flag, for missing corr pulsers and such things)
}

void Unpacker::SetTmStpLsb(unsigned long value){
  unp_data.tm_stp_lsb= value;
}

void Unpacker::SetInfoField(unsigned long value){
  unp_data.info_field= value;
}

void Unpacker::SetADCdata(unsigned int value){
  unp_data.adc_data= value;
}

void Unpacker::SetSampleLength(unsigned int value){
  unp_data.sample_length= value;
}

//validity test done here or with specific functions?
void Unpacker::SetDataType(unsigned char value){
  //if(value>=0 && value<4) unp_data.data_type= value;
  //else unp_data.data_type= 99;
  unp_data.data_type= value;
}

void Unpacker::SetFee64Id(unsigned char value){
  //if(value>=0 && value<17) unp_data.fee64_id= value;
  //else{
  //  if(b_debug) std::cout << "db       Unpacker.cpp: Invalid FEE64 momdule id: " << value << std::endl;
  //  unp_data.fee64_id= 0;
  // }
  unp_data.fee64_id= value;
}

void Unpacker::SetChId(unsigned char value){
  unp_data.ch_id= value;
}

void Unpacker::SetADCrange(unsigned char value){
  unp_data.adc_range= value;
}

void Unpacker::SetInfoCode(unsigned char value){
  unp_data.info_code= value;
}

bool Unpacker::GetBSyncStatus(){ return b_sync_status; }
//bool Unpacker::GetBCorrStatus(){ return b_corr_status; }
bool Unpacker::GetBDebug(){ return b_debug; }

bool Unpacker::GetBPushData(){ return b_push_data; }
bool Unpacker::GetBFillTree(){ return b_fill_tree; }

bool Unpacker::GetBHistograms(){ return b_histograms; }
bool Unpacker::GetBRootTree(){ return b_root_tree; }

unsigned char Unpacker::GetASICid(){ return asic_id; }
unsigned char Unpacker::GetFastDiscHits(){ return fast_disc_hits; }

unsigned long  Unpacker::GetTmStp(){ return unp_data.tm_stp; }
unsigned long  Unpacker::GetCorrScaler(){ return unp_data.corr_scaler; }
unsigned long  Unpacker::GetTmStpLsb(){ return unp_data.tm_stp_lsb; }
unsigned long  Unpacker::GetInfoField(){ return unp_data.info_field; }
unsigned int  Unpacker::GetADCdata(){ return unp_data.adc_data; }
unsigned int  Unpacker::GetSampleLength(){ return unp_data.sample_length; }
unsigned char  Unpacker::GetDataType(){ return unp_data.data_type; }
unsigned char  Unpacker::GetFee64Id(){ return unp_data.fee64_id; }
unsigned char  Unpacker::GetChId(){ return unp_data.ch_id; }
unsigned char  Unpacker::GetADCrange(){ return unp_data.adc_range; }
unsigned char  Unpacker::GetInfoCode(){ return unp_data.info_code; }

bool  Unpacker::GetSyncFlag(){ return unp_data.sync_flag; }
//bool  Unpacker::GetCorrFlag(){ return unp_data.corr_flag; }

Unpacker::Unpacker(){

  b_debug = false;

  b_root_tree = false;
  b_histograms = false;

  b_push_data = false;
  b_fill_tree = false;


  tm_stp_msb = 0;
  for(int i=0;i<common::N_FEE64;i++){
    corr_scaler_data0[i]= 0;
    corr_scaler_data1[i]= 0;
  }
 
  b_sync_status = false;
  //  b_corr_status = false;


  //...Parameters....
  //by default, all channels enabled
  for(int i=0;i<common::N_FEE64;i++){
    b_mod_enabled[i]= true;
  }

  unp_data.tm_stp= 0;
  unp_data.corr_scaler= 0;
  unp_data.tm_stp_lsb= 0;

  unp_data.info_field= 0;
  unp_data.adc_data= 0;
  unp_data.sample_length= 0;

  unp_data.data_type= 0;
  unp_data.fee64_id= 0;
  unp_data.ch_id= 0;
  unp_data.adc_range= 0;
  unp_data.info_code= 0;

  unp_data.sync_flag= false;
  //  unp_data.corr_flag= false; 
}
