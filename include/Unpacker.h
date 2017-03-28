#ifndef _Unpacker_H //what is going on here?
#define _Unpacker_H

#include <fstream>
#include <iostream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>

#include "DataSource.h"
#include "Common.h"

class Unpacker{
 private:

  int n_unp;      //Number of data items processed.
  
  bool b_debug;      ///< Debug option
  bool b_root_tree;  ///< Write data to file option
  bool b_histograms; ///< Hitsogram output option
  bool b_push_data;  ///< Push data to calibrator stage
  bool b_fill_tree;  ///< Fill tree with unpacked data

  TTree * out_root_tree;  ///< Output root tree

  unsigned long tm_stp_msb; ///< Most significant bits of timestamp

  unsigned long corr_scaler_data0[common::N_FEE64]; //[16]; // [common::n_fee64]? need one for each module? how to check if it falls out of correlation?
  unsigned long corr_scaler_data1[common::N_FEE64]; //[16];

  unsigned char asic_id;
  unsigned char fast_disc_hits;

  double start_t, start_clk;
  int counter;

  bool b_sync_status;
  //  bool b_corr_status; 

  //-----------------------------------------------
  //   BEGIN PARAMETERS
  //-----------------------------------------------

  bool b_mod_enabled[common::N_FEE64]; 

  //-----------------------------------------------
  //   END PARAMETERS
  //-----------------------------------------------


  // !!!------------------------------------!!!
  //
  // Note: if new variables included, they must also be added to definition
  // of output_root_tree Branch in InitUnpacker()
  //
  // !!!-------------------------------------!!!
  struct unpack_data_struct{
    unsigned long tm_stp;       // Internal AIDA timestamp 
    unsigned long corr_scaler;  // External correlation timestamp
    unsigned long tm_stp_lsb;   // Least significant bits of AIDA timestamp
    unsigned long info_field;   // Info field of AIDA data word - contains some information. See GREAT data documentation

    unsigned int adc_data;      // ADC value
    unsigned int sample_length; // ?

    unsigned char data_type;    // Data type - identifies type of data: 0 = sample waveform, 1 = sample lenght, 2 = info data, 3 = ADC data
    unsigned char fee64_id;     // FEE64 ID number
    unsigned char ch_id;        // Channel number
    unsigned char adc_range;    // ADC range: 0 = low energy, 1 = high energy
    unsigned char info_code;    // Info code of AIDA 'information' word - identifies data content of info data words. See GREAT data documentation 

    bool sync_flag;             // Bool to track synchronisation.
  } unp_data;

  // *******************************************************
  // See http://npg.dl.ac.uk/documents/edoc504/edoc504.html
  // for full documentation of AIDA GREAT data format.
  // *******************************************************

 public:

  TCanvas *cUnp1;

  TH1I *hFEE64_ADClow;       // Tracks number of low energy ADC words in each FEE
  TH1I *hFEE64_ADChigh;      // Tracks number of high energy ADC words in each FEE
  TH1I *hFEE64_Waveform;     // Tracks number of waveform words in each FEE
  TH1I *hFEE64_Info;         // Tracks number of information words in each FEE
  TH1I *hFEE64_PAUSE;        // Tracks number of PAUSE items issued by each FEE
  TH1I *hFEE64_RESUME;       // Tracks number of RESUME items issued by each FEE
  TH1I *hFEE64_SYNC100;      // Tracks number of SYNC100 items issued by each FEE

  TH2I *hInfoCode_FEE64;     // 2D plot showing number of different info code words in each FEE
  TH2I * hCh_FEE64_ADClow;   // 2D plot showing number of low enegy ADC words in each strip of each FEE
  TH2I * hCh_FEE64_ADChigh;  // 2D plot showing number of high enegy ADC words in each strip of each FEE
  TH2I * hCh_FEE64_DISC;     // 2D plot showing number of discriminator words in each strip of each FEE

  Unpacker();
  ~Unpacker(){};

  void InitUnpacker(int opt, char *file_name); 
  void Process(DataSource & my_src_data);
  void Close();
  void LoadParameters(char *file_name);

  void FillHistograms();
  bool IsValidFee64Id();
  bool IsValidFee64Id(int mod_id);

  void UpdateHistograms();
  void WriteHistograms();
  void ResetHistograms();

  void ResetData();
  //  void Write();


  void SetBDebug(bool flag);
  void SetBHistograms(bool flag);
  void SetBPushData(bool flag);
  void SetBFillTree(bool flag);
  void SetBRootTree(bool flag);
  void SetBSyncStatus(bool flag);
  void SetASICid(unsigned char value);
  void SetFastDiscHits(unsigned char value);
  //  void SetBCorrStatus(bool flag);

  //  void SetTmStp(long long value);
  void SetTmStp();
  void SetFlags();

  void SetCorrScaler(unsigned long value);
  void SetTmStpLsb(unsigned long value);
  void SetInfoField(unsigned long value);
  void SetADCdata(unsigned int value);
  void SetSampleLength(unsigned int value);
  void SetDataType(unsigned char value);
  void SetFee64Id(unsigned char value);
  void SetChId(unsigned char value);
  void SetADCrange(unsigned char value);
  void SetInfoCode(unsigned char value);

  unsigned long GetTmStp();
  unsigned long GetCorrScaler();
  unsigned long GetTmStpLsb();
  unsigned long GetInfoField();
  unsigned int GetADCdata();
  unsigned int GetSampleLength();
  unsigned char GetDataType();
  unsigned char GetFee64Id();
  unsigned char GetChId();
  unsigned char GetADCrange();
  unsigned char GetInfoCode();

  bool GetSyncFlag();
  //  bool GetCorrFlag();

  bool GetBSyncStatus();
  //  bool GetBCorrStatus();

  bool GetBPushData();
  bool GetBFillTree();
  bool GetBHistograms();
  bool GetBDebug();
  bool GetBRootTree();
  
  unsigned char GetASICid();
  unsigned char GetFastDiscHits();

};


#endif
