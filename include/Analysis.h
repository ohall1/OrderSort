/*! \file Analysis.h
 * \class Analysis
 * 
 * \brief A brief description of the class.
 *
 * A class which takes unpakced and calibrated AIDA data, builds events and carries out some analysis.
 *
 * \author A. Estrade (additions by C. Griffin)
 * \date 03/02/17
 * \version 2.0
 */

#ifndef _Analysis_H //what is going on here?
#define _Analysis_H //Good question...

#include <iostream>
#include <vector>
#include <map>

//#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>

#include "DataSource.h"
#include "Calibrator.h"
#include "Common.h"

class Analysis{
 private:

  int n_analysed;

  bool b_debug;              ///< TRUE: output debug info to terminal. FALSE: no debug output  
  bool b_histograms;         ///< TRUE: Create histograms. FASLE: Don't create histograms.
  bool b_root_tree;          ///< TRUE: Write output tree. FALSE: Don't write a tree.

  TTree * out_root_tree;     ///< Output root tree to contain

  // ----- correlations in time of different types -----
  static const int n_update= 10000000;

  double start_t, start_clk;  ///< Start variables for system time measurement.

  //-----------------------------------------------
  //   BEGIN PARAMETERS
  //-----------------------------------------------
  //
  // ----- parameters: Set/(Get?) should be private functions and only performed once when initializing step -----
  bool b_mod_enabled[common::N_FEE64];     ///< Bool to represent whether a strip should be active in analysis. Range ?
  int geo_detector[common::N_FEE64];       ///< Int to determine which DSSD a FEE is connected to. (1->common::N_DSSD. -1=not used)
  int geo_side[common::N_FEE64];           ///< Int to represent whether a FEE is p- or n- side. (p-side = 1, n-side = 0. -1 = not used) 
  int geo_strip[common::N_FEE64];          ///< Int to define position on detector FEE deals with. (1 = left/bottom half, 2 = right/top half from beam direction)
  
  int64_t t0_aida;
  int64_t t0_ext;
  int64_t t_aida_prev;

  bool b_first_ts;
  int64_t first_ts, last_ts;

  int total_evt_mult[common::N_DSSD][2];   ///< Total multiplicity of event for high[1] and low[0] energy data items.
  int total_mult[2];//  = {0};
  int total_e_dep[2];// = {0};
  
  int numbers;// = 0;
  int closes, puls_num;
  int implant_words, imp_down, imp_up, imp_entry, imp_num;
  int decay_words, dec_hits, dec_num;

  double event_time_window;                ///< Size of time window for creating events (unit?)

  int64_t t_prev;                          ///< Time of previous data word.

  int dE_i_lim;                            ///< 
  int dX_i_lim;                            ///< Max strip deviation in implant particle passing between two DSSDs
  int dE_d_lim;                            ///<
  int dX_d_lim;                            ///< Max strip deviation in decay particle passing between....
  int E_i_min;                             ///< Minimum energy of implant event
  int E_d_min;                             ///< Minimum energy of decay event
  int E_d_max;                             ///< Maximum energy of decay event     

  // double t0_i, t0_d;

  int strip_min_i[common::N_DSSD][2];      ///< Lowest number strip fired during implant event, i.e. high gain range
  int strip_max_i[common::N_DSSD][2];      ///< Highest number strip fired during implant event, i.e. high gain range
  int strip_min_d[common::N_DSSD][2];      ///< Lowest number strip fired during decay event, i.e. low gain range
  int strip_max_d[common::N_DSSD][2];      ///< Highest number strip fired during decay event, i.e. low gain range

  double e_sum_d[common::N_DSSD][2];       ///< Total energy deposited during decay event
  double e_sum[common::N_DSSD][2];


  //-----------------------------------------------
  //   END PARAMETERS
  //-----------------------------------------------

  int event_count;                         //
  double t_low_prev;                       //
  double t_high_prev;                      //
  double t_disc_prev;                      //

  bool b_pulser, b_implant, b_decay;       //

  // !!!------------------------------------!!!
  //
  // note: if new variables included, they must also be added to definition
  // of output_root_tree Branch in InitCalibrator()
  //
  // !!!-------------------------------------!!!
  
  /*! \struct 
   * A structure to contain the more detailed event information.
   *
   */
  struct event_struct{

    int64_t e_i[common::N_DSSD][2];   ///< Energy from decay data per event.
    int64_t e_d[common::N_DSSD][2];   ///< Energy from implant data per event.

    int64_t t0;                       ///< Time of first data item in event.
    int64_t t0_ext;                   ///< Time of external scaler.
    int64_t dt;                       ///< 

    int multiplicity;                 ///< Multiplicity of event.

    int n_side_i[common::N_DSSD][2];  ///< Number of implant data per detector side.
    int n_side_d[common::N_DSSD][2];  ///< Number of decay data per detector side.

    int n_det_i[common::N_DSSD];      ///< Number of implant data per detector.
    int n_det_d[common::N_DSSD];      ///< Number of decay data per detector.

    int x_i[common::N_DSSD];          ///< Identified x-strip of implant event.
    int y_i[common::N_DSSD];          ///< Identified y-strip of implant event.
    int x_d[common::N_DSSD];          ///< Identified x-strip of decay event.
    int y_d[common::N_DSSD];          ///< Identified y-strip of decay event.

    int dx_d[common::N_DSSD];         ///<
    int dy_d[common::N_DSSD];         ///<

    unsigned char decay_flag;         ///< Flag to signify decay event.
    unsigned char implant_flag;       ///< Flag to signify implant event.

  } evt_data;
  
  struct dssd_evt {
  dssd_evt():t(999),t_ext(-999),det(-999),range(-999),side(-999),strip(-999),energy(-999){}
  dssd_evt(Calibrator & cal):t(cal.GetTimeAIDA()),
      t_ext(cal.GetTimeExternal()),
      det(cal.GetDSSD()),
      range(cal.GetADCrange()),
      side(cal.GetSide()),
      strip(cal.GetStrip()),
      energy(cal.GetADCenergy()){}
    Long64_t t;
    Long64_t t_ext;
    int det;
    int range;
    int side;
    int strip;
    int energy;
  };

  /*! \struct
   * Structure to hold identified cluster events.
   */
  struct cluster_evt{
  cluster_evt():t(-999),t_ext(-999),x(-999),y(-999),z(-999),energy(-999),mult(-999),flag(-999){}
    Long64_t t;
    Long64_t t_ext;
    int x;
    int y;
    int z;
    int energy;
    int mult;
    int flag;
  };
  
  /*! \struct
   * A structure containing the 'compact' event information.
   *
   */
  struct dssd_hit{
  dssd_hit():t(-999),t_ext(-999),x(-999),y(-999),z(-999),ex(-999),ey(-999),flag(-999){}
  dssd_hit(cluster_evt x, cluster_evt y):t(-999),t_ext(-999),x(x.x),y(y.y),z(x.z),ex(x.energy),ey(y.energy),flag(x.flag){}
  dssd_hit(const dssd_hit& a):t(a.t),t_ext(a.t_ext),x(a.x),y(a.y),z(a.z),ex(a.ex),ey(a.ey),flag(a.flag){}
    Long64_t t;
    Long64_t t_ext;
    int x;
    int y;
    int z;
    int ex;
    int ey;
    int flag;
  };
  
  dssd_hit hit;
  
  typedef std::multimap<int, std::multimap<int, dssd_evt> > Hit_array;  ///< Multimap typedef to hold dssd_hits. [DSSD]<side, <strip, evt_info> >
  Hit_array decay_hits[common::N_DSSD];            ///< Array to contain decay events in each detector
  Hit_array implant_hits[common::N_DSSD];          ///< Array to contain decay events in each detector

  typedef std::multimap<int, cluster_evt> Clust_array;  ///< Multimap typedef to hold event clusters. [DSSD]<side, cluster_info>
  Clust_array decay_clusts[common::N_DSSD];
  Clust_array implant_clusts[common::N_DSSD];

  Hit_array::iterator side_it;                     ///< Iterator for Hit_array. Loops over side.
  std::multimap<int, dssd_evt>::iterator strip_it; ///< Iterator for Hit_array internal array. Loops over strips.

  typedef std::multimap<Long64_t, dssd_hit> Event_array;     ///< Multimap typedef to hold identified events. [DSSD]<timestamp, cluster_evt>
  Event_array decay_evts;
  Event_array implant_evts;

 public:

  //*************************************
  //          ADC spectra
  //*************************************
  TH2I * hADClowCh[common::N_FEE64];        ///< 2D ADC spectra for low energy range (per FEE64).
  TH2I * hADChighCh[common::N_FEE64];       ///< 2D ADC spectra for high energy range (per FEE64).
  TH2I * hADCdiscCh[common::N_FEE64];       ///< 2D spectra fast discriminator spectra (per FEE64).
  TH1I * hCh_ADClow[common::N_FEE64];       ///< 1D ADC spectra for low energy range (per FEE64).
  TH1I * hCh_ADChigh[common::N_FEE64];      ///< 1D ADC spectra for high energy range (per FEE64).
  TH1I * hCh_ADCdisc[common::N_FEE64];      ///< 1D fast discriminator spectra (per FEE64).

  TH1I * hADClow_all[common::N_DSSD][3];       ///< Combined ADC spectra of calibrated strips for individual DSSDs [DSSD#][] and different decay conditions [][j].

  TH1I * hElow[common::N_FEE64];            ///< Energy spectra for the low energy range (per FEE64). 
  TH1I * hEhigh[common::N_FEE64];           ///< Energy spectra for the high energy range (per FEE64).
  TH1I * hEdisc[common::N_FEE64];           ///< Energy spectra for the fast discriminator (per FEE64).

  //TH1I * hElowFin[common::N_DSSD];          ///< Energy spectrum for low energy range of final 'identified' decay energy
  //TH2I * hElowFinExEy[common::N_DSSD];      ///< Ex vs Ey plot for final 'identified' decay energies

  //*************************************
  //      Performance monitoring
  //*************************************
  // . SYNC/PAUS/RESUME per FEE
  // . ADC(high) per FEE
  // . ADC(low) per FEE
  // . ADC(disc) per FEE

  //TH1I * hSyncFEE;                          ///< 
  //TH1I * hPauseFEE;                         ///< 
  //TH1I * hResumeFEE;                        ///< 
  //TH1I * hADClowFEE;                        ///< 
  //TH1I * hADChighFEE;                       ///< 
  //TH1I * hADCdiscFEE;                       ///< 

  //*************************************
  //         Time Distributions
  //*************************************
  TH1I * hTimeADClow[common::N_DSSD];       ///< Timestamp distribution of lowE range ADC words in each DSSD -> t_n - t_(n-1).
  TH1I * hTimeADCdisc[common::N_DSSD];      ///< Timestamp distribution of fast disc words in each DSSD -> t_n - t_(n-1).
  TH1I * hTimeADChigh[common::N_DSSD];      ///< Timestamp distribution of highE range ADC words in each DSSD -> t_n - t_(n-1).

  TH1I * hTimeStamp;                        ///< 
  TH1I * hTimeStampExt;                     ///< 
  TH1I * hTimeStampFlag;                    ///< 

  TH1I * hEvt_TmStpDist;                    ///< Timestamp distribution within an event. [0] = decays, [1] = implants.

  //*************************************
  //            Decay events
  //*************************************
  TH1I * hEvt_Eside_d[common::N_DSSD][2];     ///< Energy spectrum of decay events for DSSDi and n-/p-sdie [i][n=0,p=1]. 
  
  TH2I * hEvt_ExEy_d[common::N_DSSD][2];      ///< Ex vs Ey for all possible front/back decay cluster pairs [0] and matched front/back pairs [1] in each detector during an event.
  TH1I * hEvt_residualE_d;                    ///< Ex - Ey for all possible front/back decay cluster pairs.

  TH2I * hEvt_XY_d[common::N_DSSD];           ///< x-y distrution of decay events for each DSSD
  TH1I * hEvt_X_d[common::N_DSSD];            ///< X-strip distribution of decay events.
  TH1I * hEvt_Y_d[common::N_DSSD];            ///< Y-strip distribution of decay events.

  TH1I * hEvt_Mult_d[common::N_DSSD][2];      ///< Multiplicity of decay events.
  TH2I * hEvt_MultXY_d[common::N_DSSD][2];    ///< x vs y multiplicity of decay events in each DSSD.

  TH2I * hEvt_EPulser_d;                      ///< E_x vs E_y for pulser events. Found by calculating (E1n+...+EN_DSSDn)/N_DSSD vs (E1p+...+EN_DSSDp)/N_DSSD.
  TH2I * hEvt_pulserMult;

  //TH2I * hEvt_ExEy_sum_d[common::N_DSSD];   ///<
  
  //TH1I * hEvt_dX_d[common::N_DSSD];         ///< Change in x-position of decay events between DSSD planes.
  //TH1I * hEvt_dY_d[common::N_DSSD];         ///< Change in y-position of decay events between DSSD planes.

  //TH2I * hEvt_ExEy_df[common::N_DSSD];      ///< E_x vs E_y for decay events in DSSDi w/  condition decay_flag>0 && decay_flay\%10==i
  //TH2I * hEvt_XY_df[common::N_DSSD];        ///< x-y plot of decay events in DSSDi w/ condition decay_flag>0 && decay_flay\%10==i
  //TH1I * hEvt_Eside_df[common::N_DSSD][2];  ///< Energy spectrum of decay events each side in DSSDi w/ condition decay_flag>0 && decay_flay\%10==i

  //TH2I * hEvt_ExEy_df2[common::N_DSSD];     ///< E_x vs E_y for decay events in DSSDi w/  condition decay_flag>10 && decay_flay\%10==i
  //TH2I * hEvt_XY_df2[common::N_DSSD];       ///< x-y plot of decay events in DSSDi, w/ condition decay_flag>10 && decay_flay\%10==i
  //TH1I * hEvt_Eside_df2[common::N_DSSD][2]; ///< Energy spectrum of decay events each side in DSSDi w/ condition decay_flag>0 && decay_flay\%10==i

  //TH1I * hEvt_flag_d;

  //TH1I * hEvt_MultiDet_d;                      ///< Tally of decay events in DSSDs.
  //TH2I * hEvt_MultiSide_d;                     ///< Tally of DSSDs/sides with decay info.
  //TH1I * hEvt_MultiStrip_d[common::N_DSSD][2]; ///< Strip multiplicity of decay events (n- [][0] and p- [][1] sides)

  //TH2I * hEvt_MultidX_d[common::N_DSSD][2];    ///< Strip multiplicity vs strip_max-strip_min by DSSD side

  //TH2I * hEvt_MultiID;                       ///< Muliplicity of implant events vs decay events.
  //TH1I * hEvt_HitsFlag;                      ///< Number of events with implant/decay flags (0:implant, 1:decay).

  //*************************************
  //            Imaplants
  //*************************************

  TH1I * hEvt_Eside_i[common::N_DSSD][2];     ///< Energy spectrum of n- ([][0]) and p-side ([][1]) events

  TH1I * hEvt_X_i[common::N_DSSD];            ///< X-strip distribution of implant events.
  TH1I * hEvt_Y_i[common::N_DSSD];            ///< Y-strip distribution of implant events.
  TH2I * hEvt_XY_i[common::N_DSSD];           ///< x-y plot of implant events.
  TH2I * hEvt_ExEy_i[common::N_DSSD][2];      ///< Ex vs Ey for all combinations of front/back implant events [0] and matched front/back pairs [1].  
  TH1I * hEvt_residualE_i;                    ///< Ex - Ey for all front/back implants pairs.
  
  TH1I * hEvt_Mult_i[common::N_DSSD][2];      ///< Multiplicity of implant events.
  TH2I * hEvt_MultXY_i[common::N_DSSD][2];    ///< x vs y multiplicity of implant events.
  TH1I * hEvt_Mult_impdec;                    ///< LEC multiplicity during implant event.
 
  //TH1I * hEvt_dX[common::N_DSSD];           ///< Change in x-strip position between DSSD planes.
  //TH1I * hEvt_dY[common::N_DSSD];           ///< Change in y-strip position between DSSD planes.
  //TH2I * hEvt_dXdX[common::N_DSSD];         ///< Change in x-strip position between DSSDi/DSSDi+1 vs. DSSDi+1/DSSDi+2
  //TH2I * hEvt_dYdY[common::N_DSSD];         ///< Change in x-strip position between DSSDi/DSSDi+1 vs. DSSDi+1/DSSDi+2

  //TH2I * hEvt_ExEy_if[common::N_DSSD];      ///< E_x vs E_y for implant events in DSSDi
  //TH2I * hEvt_XY_if[common::N_DSSD];        ///< x-y plot of implant events in DSSDi
  //TH1I * hEvt_Eside_if[common::N_DSSD][2];  ///< Energy spectrum of decay events in each side of DSSDi

  //TH1I * hEvt_Eaida;                        ///< Total energy deposited by implant events over all DSSDs (per event).
  //TH1I * hEvt_Eaida_gE;                     ///<
  //TH1I * hEvt_Eaida_gX;                     ///<

  //TH2I * hEvt_ExEy[common::N_DSSD];         ///< E_x vs E_y for each DSSD
  //TH2I * hEvt_EdE;                          ///<
  //TH2I * HEvt_EdE_gE;
  //TH2I * HEvt_EdE_gX;
  //TH2I * HEvt_EdE_g***;                   //also advanced patterns... hits first det but not last dssd?

  //TH1I * hEvt_Multi[common::N_DSSD][2];     ///< Multiplicity of implant events in DSSDi and n-/p-side [i][n=0,p=1].
  //TH2I * hEvt_MultiSide[common::N_DSSD];    ///< 
  //TH1I * hEvt_Hits_gE[4][2];              ///<
  //TH1I * hEvt_Hits_gX[4][2];              ///<

  //TH1I * hEvt_HitsSide;                     ///< Tallies hits in each side of each DSSD
  //TH1I * hEvt_HitsDet;                      ///< Tallies hits in each of the DSSDs
  //TH1I * hEvt_HitsDet_gE;                   ///<
  //TH1I * hEvt_HitsDet_gX;                   ///<

  // *************************************
  //                 Canvaes
  // *************************************

  TCanvas *cADClow[2];                       ///< 2D and 1D low energy range ADC spectra for all active FEEs. Panel [0] 4xN_DSSD, panel [1] 8xN_DSSD (XxY) 
  TCanvas *cADCdisc[2];                      ///< 2D and 1D discriminator spectra for all active FEEs. Panel [0] 4xN_DSSD, panel [1] not currently used.
  TCanvas *cADChigh[2];                      ///< 2D and 1D high energy range ADC spectra for all active FEEs. Panel [0] 4xN_DSSD, panel [1] not currently used.
  TCanvas *cEall[2];                         ///< Overlay of lowE ADC and fast disc spectra [0] and highE energy spectra [1].
  TCanvas *cTimeDist[common::N_DSSD];        ///< Timestamp diganostics per DSSD. All panels 

  TCanvas *cEvtE1;                           ///< Energy plots for implant events: n- and p-side energy, E_x vs E_y and same w/&w/o flag constraint
  //TCanvas *cEvtE2;                           ///<
  TCanvas *cEvtXY;                           ///< XY plots for implant events.
  TCanvas *cEvtdXdY;                         ///< dX/dY plots for fast ions between planes.
  TCanvas *cEvtMulti;                        ///< Multiplicity and DSSD/side tally of implant events.
  //TCanvas *cEvtHits;                         ///<
  //TCanvas *cEvtTS;                           ///<

  TCanvas *cEvtE_d;                          ///< Energy plots for decay events.
  //TCanvas *cEvtE_df;                         ///< Energy plots for decay events, with flag constraints.
  TCanvas *cEvtXY_d;                         ///< XY plots for decay events (w + w/o constraints).
  TCanvas *cEvtXY2_d;                        ///< XY plots for decay events (n- and p-sides)
  TCanvas *cEvtMulti_d;                      ///< Multiplicity of deacy events by DSSD/side and some other things.
  //TCanvas *cEvtHit_d;                        ///<

  // *******************************************
  //                  METHODS
  // *******************************************

  Analysis();
  ~Analysis(){};

  /*! \fn InitAnalysis(int opt, cahr *file_name)
   * Initialises the analysis process. Calls: ResetEvents() and LoadParamters() 
   * and initialises all histograms and canvases used.
   * \param[in] opt Integer representing the level of output to produced during the analysis stage.
   * \param[in] file_name The name of the text file containing the input parameters.
   */
  void InitAnalysis(int opt, char *file_name); //
  
 /*! \fn
   * Processes data from my_source or my_cal_data
   * \param[in] my_source DataSource type which comes from somewhere...
   * \param[in] my_cal_data Calibrated data passed from Calibrator.cpp...
   * \return void
   */
  void Process(DataSource & my_source, Calibrator & my_cal_data);
 
  /*! \fn
   * Run at the end of the analysis process. Calls the WriteHistograms() method and writes the root tree to file.
   **/
  void Close();
  
  /*! \fn
   * Loads system parameters from an input text file.
   * \param[in] file_name File name from where to take the input parameters.
   */
  void LoadParameters(char *file_name);

  /*! \fn
   * Build events from teh calibrated data passed from Calibrator.cpp.
   * \param[in] my_cal_data Calibrated raw data from Calibrator.cpp.
   * \return Bool: true if ....., false if.....
   */
  bool BuildEvent(Calibrator & my_cal_data);
  
  void BuildImplant(Calibrator & my_cal_data);

  
  void BuildDecay(Calibrator & my_cal_data);


  /*! \fn
   * Examines event proerties and determines the type of a event (decay/implant/pulser), filling
   * a dssd_hits type event to the output root tree.
   */
  void CloseEvent();
  
  /*! \fn
   * Initialise an event, calling the ResetEvent() method to reset event 
   * variables and then calling BuildEvent().
   */
  void InitEvent(Calibrator & my_cal_data);

  /*! \fn
   * Fill 'diagnostic' histograms, i.e. histograms not specifically to do with an event.
   * \param[in] my_cal_data Description of it
   */
  void FillHistogramsSingles(Calibrator & my_cal_data);
  
  /*! \fn
   * Fill event specific histograms, e.g. xy plots, energy plots etc.
   */
  void FillHistogramsEvent();
  
 /*! \fn
   * If requested by the user during the run, this method updates all canvases shown on screen
   * to display their most recent contents.
   */
  void UpdateHistograms();

  /*! \fn
   * If requested by the user during the run, the bin contents of all histograms are 
   * reset to zero. Filling will continue from this point.
   *
   */
  void ResetHistograms();

  /*! \fn
   * Checks whether a given module/channel number is valid, i.e. does it fall within the
   * expected values set out in Common.h.
   */
  bool IsValidChannel(int module, int channel);

  /*! \fn
   *
   */
  void ResetEvent();
  
  /*! \fn
   * Write the histograms and canvases to file opened in main.cpp.
   */
  void WriteHistograms();

  /*! \fn
   * Does a thing. A very important thing.
   */
  void WriteOutBuffer(DataSource & my_source);

  /*! \fn
   * Checks whether a given channel is enabled.
   * \return true if enabled, false if not.
   */
  bool IsChEnabled(Calibrator & my_cal_data);

  /*! \fn
   * Prints out some useful information.
   */
  void PrintEvent();

  // ********************************
  //            Setters...
  // ********************************
  
  /*! \fn
   * Sets value of b_debug.
   * \param[in] flag Value to which b_debug will be set.
   */
  void SetBDebug(bool flag);

  /*! \fn
   * Sets value of b_histograms.
   * \param[in] flag Value to which b_histograms will be set.
   */
  void SetBHistograms(bool flag);

  /*! \fn
   * Sets value of b_push_data. Currently unused.
   * \param[in] flag Value to which b_push_data will be set.
   */
  void SetBPushData(bool flag);

  /*! \fn
   * Currently unused.
   */
  void SetBFillTree(bool flag);

  /*! \fn
   * Sets value of b_root_tree. Determined by initial run parameters.
   * \param[in] flag Value to which b_root_tree will be set.
   */
  void SetBRootTree(bool flag);

  /*! \fn
   * Set value of the time window over which events will be created.
   * \param[in] value Value of the time window in units of 10ns.
   */
  bool SetEventTimeWindow(double value);

  // *******************************
  //            Getters...
  // *******************************

  /*! \fn
   * Get the value of the time window over which events are created (units of 10ns).
   */
  double GetEventTimeWindow();

  /*! \fn
   * Return the current value of b_debug.
   */
  bool GetBDebug();

  /*! \fn
   * Return the current value of b_histograms.
   */
  bool GetBHistograms();

  /*! \fn
   * Return the value of b_push_data. Currently unused.
   */
  bool GetBPushData();

  /*! \fn
   * Currently unused.
   */
  bool GetBFillTree();

  /*! \fn
   * Return the current value of b_root_tree.
   */
  bool GetBRootTree();

  /*! \fn
   * Return the mutliplicity (int) of the current event.
   */
  int GetMultiplicity();

  // ******** Unused *********
  //  void Update();
  //  void SetCanvas2(TCanvas * canvas);
  //  void CalibrateAdc();

};


#endif
