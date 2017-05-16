#ifndef _Calibrator_H //what is going on here?
#define _Calibrator_H

#include <fstream>
#include <iostream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2I.h>

#include "Unpacker.h"
#include "Common.h"

class Calibrator{
 private:

  int n_calib;

  bool b_debug;

  bool b_root_tree; //do we write data to file
  bool b_histograms;

  bool b_push_data;
  bool b_fill_tree;

  //TFile * out_file_root;
  TTree * out_root_tree;

  // ----- correlations in time of different types -----
  unsigned long tm_stp_disc[common::N_FEE64][common::N_CHANNEL]; //latest disc tm stp!
  Long_t tm_stp_corr_offset;

  int hist_fill_count;
  static const int n_update = 10000000;

  double start_t, start_clk;

  //-----------------------------------------------
  //   BEGIN PARAMETERS
  //-----------------------------------------------
  //
  // ----- parameters: Set/(Get?) should be private functions and only performed once when initializing step -----

  // ++ mask modules & channels ++
  //  bool b_mod_enabled[common::N_FEE64]; 
  bool b_ch_enabled[common::N_FEE64][common::N_CHANNEL];
  
  // ++ geometry of stack ++
  int map_dssd[common::N_FEE64];       ///< Map to give module -> DSSD information.
  int map_side[common::N_FEE64];       ///< Map to give module -> vertical/horizontal strips, i.e is a FEE top/bottom/L/R half of a DSSD.
  int map_strip[common::N_FEE64];      ///< Map to give module -> ch(0:63)/ch(64:127) info.
  
  // ++ ADC energy calibration ++
  int adc_polarity[common::N_FEE64];                       ///< Polarity of ADC for given module -> n-/p-side strips. 
  double adc_offset[common::N_FEE64][common::N_CHANNEL];   ///< ADC offset for given module and channel
  int adc_thresh[common::N_FEE64][common::N_CHANNEL];      ///< Threshold of ADC channel for event to be considered 'real'/above noise.
  double adc_gain[common::N_FEE64][common::N_CHANNEL];     ///< Gain in low energy range for given module and channel. For now use simple calibration - 0.7keV/ch.
  double adc_gain_highE[common::N_FEE64];                  ///< Gain in high energy range for given module and channel. For now use simple calibration.

  int disc_time_window;                ///<
  double aida_time_calib;              ///<

  //-----------------------------------------------
  //   END PARAMETERS
  //-----------------------------------------------

  //ch enables, ch valid.....???
  bool b_enabled;                      //if c hannel is not masked 
  bool b_valid_ch;                     //if .... we have set the channels for this data bit?
  bool b_corr_status;                  //
  //  bool b_corr_offset;  //keep info if we have data of correlation scaler for synchronization

  // !!!------------------------------------!!!
  //
  // note: if new variables included, they must also be added to definition
  // of output_root_tree Branch in InitCalibrator()
  //
  // !!!-------------------------------------!!!
  
  /*! \struct
   * A structure containing the calibrated AIDA data events.
   *
   */
  struct calib_data_struct{

    double time_aida;           ///< Internal timestamp of AIDA.
    double time_disc;           ///< Timestamp of discriminator event associated with slow comparator event.
    double time_external;       ///< Timestamp of external correlation scalar.
    double adc_energy;          ///< ADC energy of data event.

    int adc_data;               ///< ADC value for data event.
    int dssd;                   ///< Number of the DSSD in which a data event occured. 
    int strip;                  ///< Strip number of data event.

    unsigned char adc_range;    ///< ADC range of data event. 1 = high (20GeV) and 0 = low (20MeV/1GeV) range.
    unsigned char side;         ///< Side data event occured in. 1 = p- and 0 = n-side.
    unsigned char module;       ///< Module in which data event occured.
    unsigned char channel;      ///< Channel in which data event occured.
    unsigned char data_type;    ///< Type of data event. ADC, discriminator, SYNC100, PAUSE/RESUME etc...
    unsigned char info_code;    ///< Info code of data event.

    //    bool sync_flag; //... should always be true... (?)... so remove; if false data is useless and will not be recorded?
    bool corr_flag;             ///< True if we have a good calculation of offset EXT->AIDA.
    bool disc_flag;             ///< True if we have good DISC value for ADC hit.

  } cal_data;


 public:

  TCanvas *cCal1;               ///<  
  TCanvas *cCal2;               ///<

  TH2I *hCalEhCh[common::N_FEE64];       ///< Calibrated energy distribution of highE in each FEE.
  TH2I *hCalElCh[common::N_FEE64];       ///< Calibrated energy distribution of lowE in each FEE.
  TH1I *hCalStripCount[common::N_FEE64]; ///< Number of times each strip in each FEE fires.        

  Calibrator();
  ~Calibrator(){};


  void InitCalibrator(int opt, char *file_name);
  void Process(Unpacker & my_unp_data);
  void Close();
  void LoadParameters(char *file_name);


  void FillHistograms();
  void UpdateHistograms();
  void ResetHistograms();
  bool IsValidChannel(int module, int channel);
  void Update();
  //  void SetCanvas2(TCanvas * canvas);

  void ResetData();
  void Write();

  //get and set: which ones!

  bool SetGeometry();
  void CalibrateADC();
  unsigned char OrderChannel( unsigned char ch );

  //Setters...
  void SetBDebug(bool flag);
  void SetBHistograms(bool flag);
  void SetBPushData(bool flag);
  void SetBFillTree(bool flag);
  void SetBRootTree(bool flag);

  void SetBEnabled(bool flag); //forget for now...
  void SetBValidCh(bool flag);
  void SetBCorrStatus(bool flag);
  //void SetBCorrOffset(bool flag);

  void SetTmStpDisc(unsigned long value);
  void SetTmStpOffset(Long_t value); //time difference AIDA->EXTERNAL


  //Setters... for cal_data structure
  void SetTimeAIDA(double value);
  bool SetTimeDisc();
  void SetADCTmStp();
  void SetTimeExternal();
  void SetADCenergy(double value);

  void SetADCdata(int value);
  void SetDSSD(int value);
  void SetStrip(int value);

  void SetADCrange(unsigned char value);
  void SetSide(unsigned char value);
  void SetModule(unsigned char value);
  void SetChannel(unsigned char value);
  void SetDataType(unsigned char value);
  void SetInfoCode(unsigned char value);

  //  void SetSyncFlag(bool value);
  //  void SetCorrFlag(bool value);
  void SetCorrFlag();
  void SetDiscFlag(bool value);


  //Getters...
  bool GetBDebug();
  bool GetBHistograms();
  bool GetBPushData();
  bool GetBFillTree();
  bool GetBRootTree();

  bool GetBEnabled();
  bool GetBValidCh();
  bool GetBCorrStatus();
  //bool GetBCorrOffset();

  unsigned long GetTmStpDisc();
  Long_t GetTmStpOffset();

  //Getters... for cal_data structure
  double GetTimeAIDA();
  double GetTimeDisc();
  double GetTimeExternal();
  double GetADCenergy();

  int GetADCdata();
  int GetDSSD();
  int GetStrip();

  unsigned char GetADCrange();
  unsigned char GetSide();
  unsigned char GetModule();
  unsigned char GetChannel();
  unsigned char GetDataType();
  unsigned char GetInfoCode();

  bool GetCorrFlag();
  //  bool GetSyncFlag();
  bool GetDiscFlag();

};


#endif
