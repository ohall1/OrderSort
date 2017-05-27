/*! 
 * \file Analysis.cpp
 * \class Analysis
 *
 * Semi-commented C++ file responsible for the analysis section of the AIDA data processing.
 * A class which takes unpakced and calibrated AIDA data, builds events and carries out some analysis.
 *
 * \author A. Estrade (additions from C. Giffin)
 * \date
 */

// ************* To-do ***************
// + Look into gradual slow down
//     > memory hogging?
//     > maps not cleared?
// + Optimise use of maps
// + Writing incremental output files
// + General efficiency
//     > multiple front/back matched hits per event
// ***********************************

#include "Analysis.h"

#include <iomanip>
#include <stdio.h>
#include <sstream>
#include <string>

//the heart of the matter

/* Main process building events, closing them when finished and writing to relevant histograms.
 * 
 */
void Analysis::Process(DataSource & my_source, Calibrator & my_cal_data){

  //ResetData(); //set values of unp_data to defaults (zero)
  //channel will be set to default values for some bits of data where they don't apply (SYNC pulse)

  if(b_first_ts) {first_ts = my_cal_data.GetTimeAIDA(); b_first_ts = false;}
  else last_ts = my_cal_data.GetTimeAIDA();

  if(GetBHistograms()) FillHistogramsSingles(my_cal_data);
  
  if( BuildEvent(my_cal_data) ){  //If true, finish and build event and start new. If false, continues adding to current event.

    for(int i=0; i<common::N_DSSD; ++i) {
      for(int j=0; j<2; ++j) {
	total_e_dep[j] += e_sum_d[i][j];
	total_mult[j]  += total_evt_mult[i][j];
      }
    }
    
    if( total_mult[0] > 0 || total_mult[1] > 0 ){
       CloseEvent();
              
       if(b_debug){
	PrintEvent();
      }

      ++n_analysed;
      if(GetBHistograms()) FillHistogramsEvent();
      if(my_source.GetBSendData()) WriteOutBuffer(my_source);
    }
  
    InitEvent(my_cal_data);
    
  }
}

//the real heart of the matter

/* Build events for the calibrated data.
 * Takes timestamp of current event and if inside time window continues to add data to ongoing event.
 * If outside time window, closes current event and initialised new one.
 * Assigns basic decay/implant/pulser flag to data events.
 */
bool Analysis::BuildEvent(Calibrator & my_cal_data){

    if(tmStpCorrOffeset != my_cal_data.GetTmStpOffset()){
      tmStpCorrOffeset = my_cal_data.GetTmStpOffset();
    }


  ++numbers;

  if(!IsChEnabled(my_cal_data)) return false; //skip channels not enabled

  //Using discrim time for events - currently not enabled
  if(0){
    if( (my_cal_data.GetTimeDisc() > ( t0_aida + event_time_window) || my_cal_data.GetTimeDisc() < ( t0_aida - event_time_window)) 
	/*&& my_cal_data.GetTimeDisc() - t_disc_prev > 250 */){

      return true;  //start new event
    }
  }
  else  if( my_cal_data.GetTimeAIDA() > ( t0_aida + event_time_window) || // If above or below time window, or more than 2.5us since last data word,
	    my_cal_data.GetTimeAIDA() < ( t0_aida - event_time_window) || // return true, i.e. start new event.
	    ((my_cal_data.GetTimeAIDA() - t_aida_prev) > 250 && t_aida_prev != -999) ){
    return true; //time to start a new event!
  }

  //If good data to add to current event....

  t_aida_prev = my_cal_data.GetTimeAIDA();
  evt_data.dt = my_cal_data.GetTimeAIDA() - evt_data.t0; //assume monotonically increasing tm-stps
  if( (int)my_cal_data.GetADCrange() == 1 )        {BuildImplant(my_cal_data); ++implant_words;}
  else if( (int)my_cal_data.GetADCrange() == 0 )   {BuildDecay(my_cal_data); ++decay_words;}


  // -----------------------------------------------------------
  //   Fill some histograms with internal event data
  // ----------------------------------------------------------
  if(GetBHistograms()) hEvt_TmStpDist->Fill( my_cal_data.GetTimeAIDA() - t0_aida );
  /*if(GetBHistograms()){
    
    hEvt_TmStpDist[0]->Fill( evt_data.dt );
    hEvt_TmStpDist[3]->Fill( my_cal_data.GetTimeAIDA() - evt_data.t0 );
    
    if(my_cal_data.GetDiscFlag()) hEvt_TmStpDist[2]->Fill( my_cal_data.GetTimeDisc() - evt_data.t0 );
    
    }*/
  return false;
}

 void Analysis::BuildImplant(Calibrator & cal){
  
   dssd_evt imp_evt;
   //if(cal.GetDiscFlag())        dec_evt.t = cal.GetTimeDisc();
   //else if(!cal.GetDiscFlag())  dec_evt.t = cal.GetTimeAIDA();
   imp_evt.t      = cal.GetTimeAIDA();
   imp_evt.t_ext  = cal.GetTimeExternal();
   imp_evt.det    = cal.GetDSSD();
   imp_evt.range  = cal.GetADCrange();
   imp_evt.side   = cal.GetSide();
   imp_evt.strip  = cal.GetStrip();
   imp_evt.energy = cal.GetADCenergy();
  
   if(imp_evt.det < 1) {
   std::cout << "***ANA IMP DSSD# < 1***" << std::endl;
   std::cout << "Evt data > t: " << std::fixed << cal.GetTimeAIDA() << "\t det: " << (int)cal.GetDSSD() << "\t side: " << (int)cal.GetSide() << "\t strip: " << (int)cal.GetStrip() << "\t e: " << (int)cal.GetADCenergy() << std::endl;
   }
   
   if(imp_evt.det > 0) {
     implant_hits[cal.GetDSSD()-1][cal.GetSide()].insert({cal.GetStrip(),imp_evt});
    // std::cout << implant_hits[cal.GetDSSD()-1][cal.GetSide()].size() << " hello" << std::endl;

     e_sum_d[cal.GetDSSD()-1][cal.GetADCrange()] += cal.GetADCenergy();
     ++total_evt_mult[cal.GetDSSD()-1][cal.GetADCrange()];
   }
 }

void Analysis::BuildDecay(Calibrator & cal){
  
  dssd_evt dec_evt;
  //if(cal.GetDiscFlag())        imp_evt.t      = cal.GetTimeDisc();
  //else if(!cal.GetDiscFlag())  imp_evt.t_ext  = cal.GetTimeAIDA();
  dec_evt.t      = cal.GetTimeAIDA();
  dec_evt.t_ext  = cal.GetTimeExternal();
  dec_evt.det    = cal.GetDSSD();
  dec_evt.range  = cal.GetADCrange();
  dec_evt.side   = cal.GetSide();
  dec_evt.strip  = cal.GetStrip();
  dec_evt.energy = cal.GetADCenergy();
  
  if(dec_evt.det < 1) {
    std::cout << "***ANA DEC DSSD# < 1***" << std::endl;
    std::cout << "Evt data > t: " << std::fixed << cal.GetTimeAIDA() << "\t det: " << (int)cal.GetDSSD() << "\t side: " << (int)cal.GetSide() << "\t strip: " << (int)cal.GetStrip() << "\t e: " << (int)cal.GetADCenergy() << std::endl;
  }
    
  if(dec_evt.det > 0) {
    decay_hits[cal.GetDSSD()-1][cal.GetSide()].insert({cal.GetStrip(),dec_evt});
    
    e_sum_d[cal.GetDSSD()-1][cal.GetADCrange()] += cal.GetADCenergy();
    ++total_evt_mult[cal.GetDSSD()-1][cal.GetADCrange()];
  }

}

/* Determines type of event and writes event info to output tree
 *
 */
void Analysis::CloseEvent(){

  ++closes;

  #ifdef DEBUG_ANA
  std::cout << "++CLOSE EVENT++ >> " << std::endl << 
    "\t total implant E: " << total_e_dep[1] << std::endl << 
    "\t total implant m: " << total_mult[1] << std::endl <<
    "\t total decay E: " << total_e_dep[0] << std::endl << 
    "\t total decay m: " << total_mult[0] << std::endl;
  #endif

  bool b_evt_down;  //Booleans for up/down stream impant events
  bool b_evt_up;

  b_implant = false; //Boolean for implant/decay and pulser events
  b_decay = false;
  b_pulser = false;

  int NmaxI=4;        //Max number of strips to be triggered in implant event
  int EminI=500;      //Min Edep for event to be considered real (keV)
  int N_max_decay=8;  //Man number of strips to be triggered in decay event

  // *********************
  //   IMPLANTS IMPLANTS
  // *********************
  for(int det=0; det < common::N_DSSD; det++) {  
    b_evt_down = false; //Initialises upstream and downstream conditions for positive implant
    b_evt_up = true;
  
    //Check detector has events and no downstream events
    if( implant_hits[det][0].size() > 0 || implant_hits[det][1].size() > 0 ) {
      //std::cout << "    **** " << implant_hits[det].size() << " new implant events in DSSD " << det+1 << " ****" << std::endl;
      for (int down_det = det+1; down_det < common::N_DSSD; ++down_det) {
	      if ( implant_hits[down_det][0].size() > 0 || implant_hits[down_det][1].size() > 0 ) {
	          b_evt_down = true;
	          break;
        }
      }
    
      //If no events downstream check for events in all upstream detectors. 
      if (!b_evt_down) {
	      ++imp_down;
	      for (int up_det = 0; up_det < det; up_det++) {
	        if( implant_hits[up_det][0].size() < 1 && implant_hits[up_det][1].size() < 1) {
	          b_evt_up = false;
	          break;
	        }
	      }

		
	//If no events downstream and events in all upstream dets, make clusters in final DSSD.
	if(b_evt_up || det == 0) {
	  ++imp_up;


	  	    
	  //std::cout << "size of DSSD map: " << implant_hits[det].size() << std::endl;

	  //for( side_it = implant_hits[det].begin(); side_it != implant_hits[det].end(); ++side_it ) { //Loop over side
    for(int hitSide = 0 ; hitSide <2 ; hitSide++){ //Loop over side

      //Initialise variables for start of cluster checks
      int multi      = 0;        //Total number of strips in each cluster
      int cluster_e  = 0;        //Total energy deposited in each cluster
      int max_e      = -999;     //Maximum energy deposited in one strip in cluster
      int max_strip  = -999;     //Strip containing maximum energy deposit
      int strip_prev = -999;     //Strip previoulsy processed
      int side_prev  = -999;     //Side previously processed
      int n_hits     = 0;        //
      Long64_t t_min = -999;     //Lowest timestamp within cluster
      Long64_t t_ext_min = -999; //Lowest t_ext within cluster

      for( strip_it = implant_hits[det][hitSide].begin() ; strip_it != implant_hits[det][hitSide].end() ; ++strip_it){
        ++n_hits;
        //std::cout << "Implant n_hit " << n_hits << " >> det: " << (int)(strip_it->second).det << "\t side: " << (int)(side_it->first) << "\t strip: " << (int)((strip_it->second).strip) << "\t e: " << (int)((strip_it->second).energy) << std::endl;
        /*for(int count = 0; count <6 ;count++){
          std::cout << count << " " << det << " " << implant_hits[count][0].size() <<std::endl;
        }*/

        if( (abs(strip_it->first  - strip_prev) == 1 || strip_prev == -999) ){ //If first strip or neighbouring strips on same side, add to cluster
          cluster_e += (strip_it->second).energy;
          ++multi;
    
          if( (strip_it->second).energy > max_e)                     { max_e = (strip_it->second).energy; max_strip = (strip_it->second).strip;} //Max_e pos
          if( (strip_it->second).t < t_min || t_min < 0)             { t_min = (strip_it->second).t;}                                            //Earliest tm_stmp in cluster
          if( (strip_it->second).t_ext < t_ext_min || t_ext_min < 0) { t_ext_min = (strip_it->second).t_ext;}                                    //Earliest external tm_stmp in cluster
          strip_prev = strip_it->first;
          side_prev = hitSide;
        } //End cluster building

        if( (abs((strip_it->first - strip_prev)) > 1) ) { //If end of cluser write to cluster map
        
          cluster_evt clust;
          clust.t      = t_min;
          clust.t_ext  = t_ext_min;
          if(side_prev == 0)      {clust.x = -1; clust.y = max_strip;}
          else if(side_prev == 1) {clust.x = max_strip; clust.y = -1;}
          clust.z      = det+1;
          clust.energy = cluster_e;
          clust.mult   = multi;
          clust.flag   = 10;

          if (clust.z > implantMaxZ){
            implantMaxZ = clust.z;
          }
    

          ++imp_entry;
    
          implant_clusts[det].insert({side_prev,clust});
        
          //std::cout << "Side implant cluster identified in DSSD " << clust.z << " >> side: " << side_prev << "\t strip: " << max_strip << "\t mult: " << multi << "\t e: " << cluster_e  << "\t t: " << t_min << std::endl;

          multi      = 1;    //
          cluster_e  = (strip_it->second).energy;    //
          max_e      = (strip_it->second).energy;    //
          max_strip  = strip_it->first;              // Reset cluster variables to
          strip_prev = strip_it->first;              // those of current hit event.
          side_prev  = hitSide;               //
          t_min      = (strip_it->second).t;         //
          t_ext_min  = (strip_it->second).t_ext;     //   
        } //End of cluster writing for side change or discontinuous strips firing

        if( n_hits == implant_hits[det][hitSide].size() ) {
          cluster_evt clust;
          clust.t      = t_min;
          clust.t_ext  = t_ext_min;
          if(side_prev == 0)      {clust.x = -1; clust.y = max_strip;}
          else if(side_prev == 1) {clust.x = max_strip; clust.y = -1;}
          clust.z      = det+1;
          clust.energy = cluster_e;
          clust.mult   = multi;
          clust.flag   = 10;

          if (clust.z > implantMaxZ){
            implantMaxZ = clust.z;
          }
    
          ++imp_entry;
    
          implant_clusts[det].insert({side_prev,clust});
        
          //std::cout << "Side implant cluster identified in DSSD " << clust.z << " >> side: " << side_prev << "\t strip: " << max_strip << "\t mult: " << multi << "\t e: " << cluster_e <<"\t t: " << t_min << std::endl;
          //if(n_hits == implant_hits[det].size()) std::cout << "End of hits in this event!!! " << std::endl;
        } //End clustering on last event
      }//End loop over side
	  } //End clustering over this DSSD
	} //End if events upstream
      } //End if events downstream
    } //End if implants in DSSDx
  }//End cluster finding over all DSSDs

  //Loop through front and back cluster events to find Ex ~ Ey events
  for(int det=0; det<common::N_DSSD;++det) {
    std::pair<Clust_array::iterator, Clust_array::iterator> y_it = implant_clusts[det].equal_range(0); //Select y (back) strip clusters
    
    for(Clust_array::iterator y_clusts = y_it.first; y_clusts != y_it.second; ++y_clusts) {
      std::pair<Clust_array::iterator, Clust_array::iterator> x_it = implant_clusts[det].equal_range(1); //select x (front) strip clusters
      
      for(Clust_array::iterator x_clusts = x_it.first; x_clusts != x_it.second; ++x_clusts) {
	if(b_histograms) {
	  hEvt_ExEy_i[det][0]->Fill( (x_clusts->second).energy, (y_clusts->second).energy );
	  hEvt_MultXY_i[det][0]->Fill( (x_clusts->second).mult, (y_clusts->second).mult );
	  hEvt_residualE_i->Fill( (x_clusts->second).energy - (y_clusts->second).energy );
	}
	if( (y_clusts->second).energy > ((x_clusts->second).energy - 500) && (y_clusts->second).energy < ((x_clusts->second).energy + 500) ) {
	  if(b_histograms) {
	    hEvt_Eside_i[det][0]->Fill((y_clusts->second).energy);
	    hEvt_Eside_i[det][1]->Fill((x_clusts->second).energy);
	    hEvt_ExEy_i[det][1]->Fill( (x_clusts->second).energy, (y_clusts->second).energy );
	    hEvt_XY_i[det]->Fill( x_clusts->first, y_clusts->first );
	    hEvt_X_i[det]->Fill( x_clusts->first );
	    hEvt_Y_i[det]->Fill( y_clusts->first );
	    hEvt_Mult_i[det][0]->Fill( (y_clusts->second).mult );
	    hEvt_Mult_i[det][1]->Fill( (x_clusts->second).mult );
	    hEvt_MultXY_i[det][1]->Fill( (x_clusts->second).mult, (y_clusts->second).mult );
	  }
	  
	  dssd_hit evt(x_clusts->second,y_clusts->second);
	  if( (x_clusts->second).t < (y_clusts->second).t && ((x_clusts->second).t > 0) && ((y_clusts->second).t > 0)){
      evt.t = (x_clusts->second).t; evt.t_ext = (x_clusts->second).t_ext;
    }
	  else if( (y_clusts->second).t < (x_clusts->second).t && ((x_clusts->second).t > 0) && ((y_clusts->second).t > 0)){
      evt.t = (y_clusts->second).t; evt.t_ext = (y_clusts->second).t_ext;
    }
    else if( (x_clusts->second).t == (y_clusts->second).t){
      evt.t = (x_clusts->second).t; evt.t_ext = (x_clusts->second).t_ext;
    }
    else{
      evt.t = 0; evt.t_ext = 0;
    }
	  implant_evts.insert({evt.t,evt});
	  ++imp_num;
	  b_implant = true;
	  //std::cout << "Implant pair -> front... det: " << det+1 << "\t x/y: " << (x_clusts->second).x << "/" << (y_clusts->second).y << "\t ex/ey: " << (x_clusts->second).energy << "/" << (y_clusts->second).energy << "\t t: " << evt.t << std::endl;
	}
      }
    }
  }
  
  // ******************************************
  //             PULSER PULSER
  // ******************************************
  
  if(!b_implant){ //If not an implant event
    if( (total_evt_mult[0][0] > 100 && /*total_evt_mult[1][0] > 100  &&*/ total_evt_mult[2][0] > 100  && total_evt_mult[3][0] > 100 && total_evt_mult[4][0] > 100 && total_evt_mult[5][0] > 100) 
	|| 
	(e_sum_d[0][0] > 5000 && /* e_sum_d[1][0] > 5000 && */e_sum_d[2][0] > 5000 && e_sum_d[3][0] > 5000 && e_sum_d[3][0] > 5000 && e_sum_d[4][0] > 5000 && e_sum_d[5][0] > 5000) ) {
      if(GetBHistograms()){hEvt_pulserMult->Fill(total_evt_mult[0][0], total_evt_mult[0][1]);}
      ++puls_num;
      b_pulser = true;
    }
  }
  
  // *******************************************
  //              DECAYS DECAYS
  // *******************************************

  if(!b_pulser && !b_implant){    // If not a pulser or implant
    
    for(int det=0; det < common::N_DSSD; ++det) {  // Loop over detectors ( 0 -> common::N_DSSD )
      if ( decay_hits[det][0].size() > 0 && decay_hits[det][1].size() > 0){
      for(int hitSide = 0 ; hitSide < 2 ; ++hitSide){
        int multi      = 0;        //Total number of strips in each cluster
        int cluster_e  = 0;        //Total energy of each cluster
        int max_e      = -999;     //Max individual strip energy
        int max_strip  = -999;     //Stip with highest energy deposit
        int strip_prev = -999;     //Previously processed strip
        int side_prev  = -999;     //Previously processed side
        Long64_t t_min = -999;     //Earliest timestamp for cluster strips
        Long64_t t_ext_min = -999; //Earliest t_ext for cluster strips
        int n_hits = 0;
        //std::cout << "Size of DSSD" << det+1 << " map: " << decay_hits[det].size() << std::endl;
      
//      for( side_it = decay_hits[det].begin(); side_it != decay_hits[det].end(); ++side_it ) {
          for( strip_it = decay_hits[det][hitSide].begin() ; strip_it != decay_hits[det][hitSide].end() ; ++strip_it){  

            ++n_hits;
            //std::cout << "Decay n_hit " << n_hits << " >> det: " << (int)(strip_it->second).det << "\t side: " << (int)(side_it->first) << "\t strip: " << (int)((strip_it->second).strip) << "\t e: " << (int)((strip_it->second).energy) << std::endl;
            //std::cout << det <<"det nhits"<<n_hits << " Hello side " << side_prev << " Side now " << hitSide << " strip prev " << strip_prev << " new strip " << (strip_it->first) <<std::endl; 

            if( (abs(strip_it->first  - strip_prev) == 1 || strip_prev == -999) ){ //If first strip or neighbouring strips on same side, add to cluster
              cluster_e += (strip_it->second).energy;
              ++multi;
              
              if( (strip_it->second).energy > max_e)                     { max_e = (strip_it->second).energy; max_strip = (strip_it->second).strip;} //Max_e pos
              if( (strip_it->second).t < t_min || t_min < 0)             { t_min = (strip_it->second).t;}                                            //Earliest tm_stmp in cluster
              if( (strip_it->second).t_ext < t_ext_min || t_ext_min < 0) { t_ext_min = (strip_it->second).t_ext;}                                    //Earliest external tm_stmp in cluster
              strip_prev = strip_it->first;
              side_prev = hitSide;            
            } //End cluster building  
            else if( (abs((strip_it->first - strip_prev)) > 1) ) { //If end of cluser write to cluster map
              //if(multi<7) {  //If cluster has mulitplicity below 7 -> write it
              cluster_evt clust;
              clust.t      = t_min;
              clust.t_ext  = t_ext_min;
              if(side_prev == 0)      {clust.x = -1; clust.y = max_strip;}
              else if(side_prev == 1) {clust.x = max_strip; clust.y = -1;}
              clust.z      = det+1;
              clust.energy = cluster_e;
              clust.mult   = decay_hits[det][hitSide].size();//multi;
              clust.flag   = 1;           

                
              ++dec_hits;
                
              decay_clusts[det].insert({side_prev,clust});
              if(GetBHistograms()){hEvt_Mult_d[det][side_prev]->Fill(multi);}
              //}   //Otherwise reset variables for next cluster build.
              //std::cout << "Side decay cluster identified in DSSD " << clust.z << " >> side: " << side_prev << "\t strip: " << max_strip << "\t mult: " << multi << "\t e: " << cluster_e  << "\t t: " << t_min << std::endl;
              
              multi      = 1;    //
              cluster_e  = (strip_it->second).energy;    //
              max_e      = (strip_it->second).energy;    //
              max_strip  = strip_it->first;              // Reset cluster variables to
              strip_prev = strip_it->first;              // those of current hit event.
              side_prev  = hitSide;               //
              t_min      = (strip_it->second).t;         //
              t_ext_min  = (strip_it->second).t_ext;     //                  
            } //End of cluster writing for side change or discontinuous strips firing 

            if( n_hits == decay_hits[det][hitSide].size() /*&& multi < 7*/ ) {
              cluster_evt clust;
              clust.t      = t_min;
              clust.t_ext  = t_ext_min;
              if(side_prev == 0)      {clust.x = -1; clust.y = max_strip;}
              else if(side_prev == 1) {clust.x = max_strip; clust.y = -1;}
              clust.z      = det+1;
              clust.energy = cluster_e;
              clust.mult   = decay_hits[det][hitSide].size();//multi;
              clust.flag   = 1;
                
              ++dec_hits;
                
              decay_clusts[det].insert({side_prev,clust});
              //hEvt_Mult_d[det][side_prev]->Fill(multi);   
              //std::cout << "End of size decay cluster identified in DSSD " << clust.z << " >> side: " << side_prev << "\t strip: " << max_strip << "\t mult: " << multi << "\t e: " << cluster_e <<"\t t: " << t_min << std::endl;
                
            } //End clustering on last event
          }//End looping these strips
        }//End clustering this side
      } //If events in this det
    } //End clustering over this DSSD
  } //End if not pulser or implant event

  //Loop through front and back cluster events to find Ex ~ Ey events
  for(int det=0; det<common::N_DSSD;++det) {
    std::pair<Clust_array::iterator, Clust_array::iterator> y_it = decay_clusts[det].equal_range(0); //Select y (back) strip clusters
    
    for(Clust_array::iterator y_clusts = y_it.first; y_clusts != y_it.second; ++y_clusts) {
      std::pair<Clust_array::iterator, Clust_array::iterator> x_it = decay_clusts[det].equal_range(1); //select x (front) strip clusters
      
      for(Clust_array::iterator x_clusts = x_it.first; x_clusts != x_it.second; ++x_clusts) {
	if(b_histograms) {
	  hEvt_ExEy_d[det][0]->Fill( (x_clusts->second).energy, (y_clusts->second).energy );
	  hEvt_MultXY_d[det][0]->Fill((x_clusts->second).mult, (y_clusts->second).mult );
	  hEvt_residualE_d->Fill( (x_clusts->second).energy - (y_clusts->second).energy );
	}
	if( (y_clusts->second).energy > ((x_clusts->second).energy - 140) && (y_clusts->second).energy < ((x_clusts->second).energy + 140)) {
	  if(b_histograms) {
	    hEvt_Eside_d[det][0]->Fill((y_clusts->second).energy);
	    hEvt_Eside_d[det][1]->Fill((x_clusts->second).energy);
	    hEvt_ExEy_d[det][1]->Fill( (x_clusts->second).energy, (y_clusts->second).energy );
	    hEvt_XY_d[det]->Fill( x_clusts->first, y_clusts->first );
	    hEvt_X_d[det]->Fill( x_clusts->first );
	    hEvt_Y_d[det]->Fill( y_clusts->first );
	    hEvt_Mult_d[det][0]->Fill( (y_clusts->second).mult );
	    hEvt_Mult_d[det][1]->Fill( (x_clusts->second).mult );
	    hEvt_MultXY_d[det][1]->Fill( (x_clusts->second).mult, (y_clusts->second).mult );
	  }

	  dssd_hit evt(x_clusts->second,y_clusts->second);
	  if( (x_clusts->second).t < (y_clusts->second).t && ((x_clusts->second).t > 0) && ((y_clusts->second).t > 0)){
      evt.t = (x_clusts->second).t; evt.t_ext = (x_clusts->second).t_ext;
    }
	  else if( (y_clusts->second).t < (x_clusts->second).t && ((x_clusts->second).t > 0) && ((y_clusts->second).t > 0)){
      evt.t = (y_clusts->second).t; evt.t_ext = (y_clusts->second).t_ext;
    }
    else if ((x_clusts->second).t == (y_clusts->second).t){
      evt.t = (x_clusts->second).t; evt.t_ext = (x_clusts->second).t_ext;
    }
	  decay_evts.insert({evt.t,evt});
	  ++dec_num;
	  b_decay = true;
	  //std::cout << "Decay pair -> det: " << det+1 << "\t x/y: " << (x_clusts->second).x << "/" << (y_clusts->second).y << "\t ex/ey: " << (x_clusts->second).energy << "/" << (y_clusts->second).energy << "\t t: " << evt.t << std::endl;
	     }
      }
    }
  }
  
  if(GetBRootTree()){
    
    if( (b_decay || b_implant) &&  !b_pulser ){     // If implant event present... Write implant event and ignore decay

      if(b_implant) {
	for(Event_array::iterator it=implant_evts.begin(); it!=implant_evts.end(); ++it) {
	  /*//Set info
	  hit.t     = (it->second).t;
	  hit.t_ext = (it->second).t_ext;
	  hit.x     = (it->second).x;
	  hit.y     = (it->second).y;
	  hit.z     = (it->second).z;
	  hit.ex    = (it->second).ex;           // p-side strips
	  hit.ey    = (it->second).ey;           // n-side strips
    hit.nx    = (it->second).nx;           // p-side strips
    hit.ny    = (it->second).ny;           // n-side stripa
    hit.flag  = 4;//(it->second).flag;       //1, 2, ..., 6. Implant flag from Calibrator.cpp 10+det# (i.e. DSSD 1 = det0) (?)*/
    root_evt root_imp((it->second));
    root_imp.T = ((it->second).t + tmStpCorrOffeset)*10;
    root_imp.ID = 4; // ID 4 defines implant 5 defines decay
    //std::cout << root_imp.T <<std::endl;
    hit = root_imp;

    if (tmStpCorrOffeset != 0){
	   out_root_tree->Fill();  // Write to tree
    }
	  hEvt_Mult_impdec->Fill(total_evt_mult[(it->second).z][0]);
	  //std::cout << "Implant event written to tree." << std::endl;
	}
      }
      else if( b_decay && !b_implant && !b_pulser ){    // Else if decay event present... (do the same)
	for(Event_array::iterator it=decay_evts.begin(); it!=decay_evts.end(); ++it){
    /*
	  hit.t     = (it->second).t;
	  hit.t_ext = (it->second).t_ext;
	  hit.x     = (it->second).x;
	  hit.y     = (it->second).y;
	  hit.z     = (it->second).z;
	  hit.ex    = (it->second).ex;           // p-side strips
	  hit.ey    = (it->second).ey;           // n-side strips
    hit.nx    = (it->second).nx;           // p-side strips
    hit.ny    = (it->second).ny;           // n-side stripa
	  hit.flag  = 5;//(it->second).flag;       //1, 2, ..., 6. Implant flag from Calibrator.cpp 10+det# (i.e. DSSD 1 = det0) (?)*/
    root_evt root_dec((it->second));
    root_dec.T = ((it->second).t + tmStpCorrOffeset)*10;
    root_dec.ID = 5; // ID 4 defines implant 5 defines decay
    hit = root_dec;
	  if(tmStpCorrOffeset != 0 && root_dec.ny < 10 && root_dec.nz <10 && root_dec.z > implantMaxZ){
	   out_root_tree->Fill();
    }
	}
      }
    }
  }
}

/* Initialise event
 * 
 */
void Analysis::InitEvent(Calibrator & my_cal_data){
    
  ResetEvent();

  implantMaxZ = 0;

  //check again just in case we're trying to initialize with wrong data
  if(/*!my_cal_data.GetADCrange() || */!IsChEnabled(my_cal_data)){
    t0_aida = -9999999; //bogus timestamp... should take care of things
  }
  /*
  if(my_cal_data.GetDiscFlag()){
    t0_aida = my_cal_data.GetTimeDisc();
    t0_ext  = my_cal_data.GetTimeExternal();

    }*/
  else{
    t0_aida = my_cal_data.GetTimeAIDA();
    t0_ext  = my_cal_data.GetTimeExternal();
    
    t_aida_prev = my_cal_data.GetTimeAIDA();
  }

  #ifdef DEBUG_ANA
  std::cout << std::endl << std::endl << "***********INIT_EVENT**********" << std::endl << "\t t0_aida: " << t0_aida << "\t t: " << my_cal_data.GetTimeAIDA() << std::endl;
  //if(my_cal_data.GetADCrange() || IsChEnabled(my_cal_data)) std::cout << "ENABLED" << std::endl;
  #endif

  //evt_data.t0= my_cal_data.GetTimeAIDA();

  BuildEvent(my_cal_data);


}

/* Not sure about this one....
 *
 */
void Analysis::WriteOutBuffer(DataSource & my_source){

  //int s_double= sizeof(double);
  //int s_int= sizeof(int);
  //int s_char= sizeof(char);

  int offset=my_source.GetBuffOffset(); 
  //int 4
  //char 1

  //  int my_int[4][2];
  //  for(int i=0;i<4;i++){
  //    my_int[i][0]= evt_data.e_i[1][0];
  //  }

  int32_t aida_id= 0xA1DA;

  memcpy(my_source.BufferOut+offset, (char*) &aida_id, sizeof(int32_t) );
  memcpy(my_source.BufferOut+offset+sizeof(int32_t), (char*) &evt_data, sizeof(evt_data) );


  #ifdef DEBUG_ANA2
    std::cout << "\n size of evt_data: "<< sizeof(evt_data) << std::endl ;
    
    int j=0;
    for(int i= offset; i< offset+sizeof(evt_data)+sizeof(int32_t); ++i){
      if((j%16)==0) std::cout<< std::endl <<" -- "; //printf("\n");
      if( (j%4)==0) std::cout << " 0x";
      
      printf("%02hhx",my_source.BufferOut[i]);
      
      ++j;
    }
  #endif

  offset= offset + sizeof(evt_data)+sizeof(int32_t);

  my_source.SetBuffOffset(offset);
  my_source.TransferBuffer(evt_data.t0);

  return;

}

/* Fills "diagnostic" histograms, e.g. not related to a specific event.
 *
 */
void Analysis::FillHistogramsSingles(Calibrator & my_cal_data){
  
  int mod   = my_cal_data.GetModule();
  int det   = my_cal_data.GetDSSD();
  int ch    = my_cal_data.GetChannel();
  int side  = my_cal_data.GetSide();
  int range = my_cal_data.GetADCrange();

  if(b_debug) std::cout << " about to fill histograms-singles....mod= "<< mod << " range= "<<range  << std::endl;

  if(b_mod_enabled[mod-1]){
    /*if(det<0 || det>common::N_DSSD){
      std::cout << "\n******************     det" << det<< "->0  ***********************"<<std::endl;
      det=0;
      }*/
    
    if(ch<0 || ch>63){
      std::cout << "\n******************     ch" << ch<< "->0  ***********************"<<std::endl;
      ch=0;
    }
    
    if(mod<1 || mod>common::N_FEE64){
      std::cout << "\n******************     mod" << mod<< "->0  ***********************"<<std::endl;
      mod=0;
    }
    
    if(side<0 || side>2){
      std::cout << "\n******************     side" << side<< "->0  ***********************"<<std::endl;
      side=0;
      }
    
    //DECAY DECAY DECAY
    if(range==0){
      hADClowCh[mod-1]->Fill(ch,my_cal_data.GetADCdata());
      hCh_ADClow[mod-1]->Fill(ch);
      hElow[mod-1]->Fill(my_cal_data.GetADCenergy());

      if(my_cal_data.GetDiscFlag()) hEdisc[mod-1]->Fill(my_cal_data.GetADCdata());
      
      if(my_cal_data.GetDiscFlag()){
	hADCdiscCh[mod-1]->Fill(ch,my_cal_data.GetADCdata());
	hCh_ADCdisc[mod-1]->Fill(ch);
      }
      
      hTimeADClow[0]->Fill(my_cal_data.GetTimeAIDA()-t_low_prev);
      t_low_prev= my_cal_data.GetTimeAIDA();
      
      if(my_cal_data.GetDiscFlag()){
	hTimeADCdisc[0]->Fill(my_cal_data.GetTimeDisc()-t_disc_prev);
	t_disc_prev= my_cal_data.GetTimeDisc();
      }
    }

    //IMPLANT IMPLANT IMPLANT
    else if(range==1){
      hADChighCh[mod-1]->Fill(ch,my_cal_data.GetADCdata());
      hCh_ADChigh[mod-1]->Fill(ch);
      hEhigh[mod-1]->Fill(my_cal_data.GetADCdata());
      
      hTimeADChigh[0]->Fill(my_cal_data.GetTimeAIDA()-t_high_prev);
      t_high_prev= my_cal_data.GetTimeAIDA();
    }
    
    if(0) std::cout << " TS(aida), TS(ext): " << my_cal_data.GetTimeAIDA()<<",  "<< my_cal_data.GetCorrFlag() << ":  "<< my_cal_data.GetTimeExternal()<<std::endl;
    
    hTimeStamp->Fill( my_cal_data.GetTimeAIDA() );
    
    if(my_cal_data.GetCorrFlag()){
      hTimeStampExt->Fill( my_cal_data.GetTimeExternal() );
      
      hTimeStampFlag->Fill( 1 );
    }
    else  hTimeStampFlag->Fill( 0 );
    
  }
}

/* Fills histograms related to specific events, e.g. xp or energy plots etc.
 *
 */
/*void Analysis::FillHistogramsEvent(){

  double e_det[common::N_DSSD] = {0}; // Highest single E_dep per DSSD in FillHistogramsEvents()
  double e_aida = 0;                  // Total E_dep from implants over all DSSDs in FillHistogramsEvents()
  //bool b_gE = true;                   //
  //bool b_gX = true;                   //

  int multi_d = 0;                    // Total number of decay events across all detectors in FillHistogramsEvents()
  int multi_i = 0;                    // Total number of implant events across all detectors in FillHistogramsEvents()

  if(b_debug) std::cout<<" ...evt histogrms..." << std::endl;

  for(int i=0;i<common::N_DSSD; ++i) {

    multi_i += evt_data.n_det_i[i];    //Add to cumulative number of events
    multi_d += evt_data.n_det_d[i];    //

    if(evt_data.e_i[i][0]>evt_data.e_i[i][1]) e_det[i]=evt_data.e_i[i][0];  //Find highest E_dep in each DSSD
    else e_det[i]=evt_data.e_i[i][1];
    
    if(evt_data.n_det_i[i]>0) e_aida += e_det[i];    //Total E_dep over all DSSDs
  }  //End DSSD loop - summation calcs
  
  hEvt_MultiID->Fill(multi_i,multi_d);
  if(evt_data.implant_flag>0) hEvt_HitsFlag->Fill(0);
  if(evt_data.decay_flag>0) hEvt_HitsFlag->Fill(1);

  if(b_debug) {
    // PrintEvent();
    printf("        eaida %f,    edet[i]  %f  %f  %f  %f \n", e_aida, e_det[0], e_det[1], e_det[2], e_det[3]);
  }

  if(multi_i>0) {   //If implant events....

    for(int i=0;i<common::N_DSSD; ++i){
      
      for(int j=0; j<2; ++j){
	if(evt_data.n_side_i[i][j]>0){
	  hEvt_Eside[i][j]->Fill( evt_data.e_i[i][j] );
	  hEvt_Multi[i][j]->Fill(evt_data.n_side_i[i][j]);
	  
	  if( (evt_data.implant_flag%10)==i+1 && evt_data.implant_flag>0) hEvt_Eside_if[i][j]->Fill( evt_data.e_i[i][j] );
	  
	}  // End if implant events in side
      }  // End sides
      
      if(evt_data.n_det_i[i]>0){    //If implant events in DSSDi+1... Fill hists.
	
	hEvt_ExEy[i]->Fill(evt_data.e_i[i][1],evt_data.e_i[i][0]);
	
	hEvt_X[i]->Fill(evt_data.x_i[i]);
	hEvt_Y[i]->Fill(evt_data.y_i[i]);
	hEvt_XY[i]->Fill(evt_data.x_i[i],evt_data.y_i[i]);
	
	if((evt_data.implant_flag%10)==i+1 && evt_data.implant_flag>0){
	  hEvt_ExEy_if[i]->Fill(evt_data.e_i[i][1],evt_data.e_i[i][0]);
	  hEvt_XY_if[i]->Fill(evt_data.x_i[i],evt_data.y_i[i]);
	}
      }

      if(evt_data.n_det_i[i]>0){
	hEvt_HitsDet->Fill(i+1,1);
      }
    } // End DSSDs
        
    hEvt_Eaida->Fill(e_aida);         //Fill total energy deposited histo
    //hEvt_EdE->Fill(e_det[0],e_aida);  //
    
    // **** Track heavy particles ****
    for(int i=0; i<common::N_DSSD; ++i) {
      
      if(evt_data.n_det_i[i] > 0 && evt_data.n_det_i[i+1] > 0 && i<common::N_DSSD-1) {   // Position change between DSSDi and DSSDi+1
	hEvt_dX[i]->Fill(evt_data.x_i[i+1] - evt_data.x_i[i]);
	hEvt_dY[i]->Fill(evt_data.y_i[i+1] - evt_data.y_i[i]);
      }

      if(evt_data.n_det_i[i] > 0 && evt_data.n_det_i[i+1] > 0 && evt_data.n_det_i[i+2] > 0 && i<common::N_DSSD-2) {  //Pos change between DSSDi+1/DSSDi vs DSSDi+2/DSSDi+1
	hEvt_dXdX[i]->Fill(evt_data.x_i[i+1] - evt_data.x_i[i], evt_data.x_i[i+2] - evt_data.x_i[i+1]);
	hEvt_dYdY[i]->Fill(evt_data.y_i[i+1] - evt_data.y_i[i], evt_data.y_i[i+2] - evt_data.y_i[i+1]);
      }
    } // End loop over DSSDs

    if(evt_data.n_det_i[0] > 0 && evt_data.n_det_i[common::N_DSSD-1] > 0) {     //Position change between first and last DSSD
      hEvt_dX[common::N_DSSD-1]->Fill(evt_data.x_i[common::N_DSSD-1] - evt_data.x_i[0]);
      hEvt_dY[common::N_DSSD-1]->Fill(evt_data.y_i[common::N_DSSD-1] - evt_data.y_i[0]);

      if(evt_data.n_det_i[common::N_DSSD/2] > 0) {     //Position change between first/middle vs middle/last
	hEvt_dXdX[common::N_DSSD-2]->Fill(evt_data.x_i[common::N_DSSD/2] - evt_data.x_i[0], evt_data.x_i[common::N_DSSD-1] - evt_data.x_i[common::N_DSSD/2]);
	hEvt_dYdY[common::N_DSSD-2]->Fill(evt_data.y_i[common::N_DSSD/2] - evt_data.y_i[0], evt_data.y_i[common::N_DSSD-1] - evt_data.y_i[common::N_DSSD/2]);
      }
    }
    // **** End particle tracking ****
    
    for(int i=0; i<common::N_DSSD; ++i){     //Fill number of hits per side
      if(evt_data.n_side_i[i][0]>0) hEvt_HitsSide->Fill(i+1,1); 
      if(evt_data.n_side_i[i][1]>0) hEvt_HitsSide->Fill(i+1+0.5,1); 
    }
        
    hEvt_TmStpDist[1]->Fill( evt_data.dt );
    
    if(b_debug) std::cout << "db     Analysis::FillHistrgramsEvent():   done with evt_implant (multi="<<multi_i<<")"<< std::endl;
    
  } //End implant events
  
  if(b_pulser) {    //If pulser event...
    hEvt_EPulser_d->Fill( (evt_data.e_d[0][0] + evt_data.e_d[1][0] + evt_data.e_d[2][0] + evt_data.e_d[3][0] + evt_data.e_d[4][0] + evt_data.e_d[5][0])/6.,
			  (evt_data.e_d[0][1] + evt_data.e_d[1][1] + evt_data.e_d[2][1] + evt_data.e_d[3][1] + evt_data.e_d[4][1] + evt_data.e_d[5][1])/6. );
  } //End pulser events
  
  else if(multi_d>0){
    
    // EVENT : DECAY : HISTOGRAMS
    int multi_det_d=0;
    
    for(int i=0;i<common::N_DSSD; ++i) {   //Loop over all DSSDs 
      
      if(evt_data.n_det_d[i]>0) ++multi_det_d;   //Counter for decay events
      
      for(int j=0; j<2; ++j) {               //Loop over sides
	if(evt_data.n_side_d[i][j]>0) {   //If side has decays events...
	  
	  hEvt_Eside_d[i][j]->Fill( evt_data.e_d[i][j] );  //Fill hist
	  
	  if(evt_data.decay_flag>0 && (evt_data.decay_flag%10) == i+1){
	    hEvt_Eside_df[i][j]->Fill( evt_data.e_d[i][j] );
	    if(evt_data.decay_flag>10) hEvt_Eside_df2[i][j]->Fill( evt_data.e_d[i][j] );
	  }
	}
      }
      
      if(evt_data.n_det_d[i]>0) {   //If DSSD has decay events...
	
	hEvt_ExEy_d[i]->Fill(evt_data.e_d[i][1],evt_data.e_d[i][0]);
	hEvt_XY_d[i]->Fill(evt_data.x_d[i],evt_data.y_d[i]);
	hEvt_flag_d->Fill(0);
	hADClow_all[i][0]->Fill(evt_data.e_d[i][1]);
	
	if(evt_data.decay_flag>0 && (evt_data.decay_flag%10) == i+1){
	  hEvt_ExEy_df[i]->Fill(evt_data.e_d[i][1],evt_data.e_d[i][0]);
	  hEvt_XY_df[i]->Fill(evt_data.x_d[i],evt_data.y_d[i]);
	  hEvt_flag_d->Fill(1);
	  hADClow_all[i][1]->Fill(evt_data.e_d[i][1]);

	  if(evt_data.decay_flag>10){
	    hEvt_ExEy_df2[i]->Fill(evt_data.e_d[i][1],evt_data.e_d[i][0]);
	    hEvt_XY_df2[i]->Fill(evt_data.x_d[i],evt_data.y_d[i]);
	    hEvt_flag_d->Fill(2);
	    hADClow_all[i][2]->Fill(evt_data.e_d[i][1]);
	  }
	}
      }
      
      hEvt_MultiDet_d->Fill(multi_det_d);
      
      if(evt_data.n_side_d[i][0]>0 && evt_data.n_side_d[i][1]>0)       hEvt_MultiSide_d->Fill(2,i+1); //If decay events in both sides
      else if(evt_data.n_side_d[i][1]>0 && evt_data.n_side_d[i][0]==0) hEvt_MultiSide_d->Fill(1,i+1); //Else if only decay events in p-side
      else if(evt_data.n_side_d[i][0]>0 && evt_data.n_side_d[i][1]==0) hEvt_MultiSide_d->Fill(0.,i+1.); //Else if only decay events in n-side
      
      if(evt_data.n_side_d[i][0]>0){
	hEvt_MultiStrip_d[i][0]->Fill(evt_data.n_side_d[i][0]);
	hEvt_MultidX_d[i][0]->Fill(evt_data.n_side_d[i][0],strip_max_d[i][0]-strip_min_d[i][0]);
		
	if(evt_data.n_side_d[i][1]>0){
	  hEvt_MultiStrip_d[i][1]->Fill(evt_data.n_side_d[i][0]);
	  hEvt_MultidX_d[i][1]->Fill(evt_data.n_side_d[i][1],strip_max_d[i][1]-strip_min_d[i][1]);
	}
      }
    }
    if(b_debug) std::cout << "db     Analysis::FillHistrgramsEvent():   done with evt_decay (multi="<<multi_d<<")"<< std::endl;
    
  }
  
  }*/

/* Update histograms and canvases to display current information during run.
 *
 */

void Analysis::FillHistogramsEvent() {

}

void Analysis::UpdateHistograms(){
  
  if(GetBHistograms()){
    std::cout << "  Analysis::UpdateHistograms()... updating"<<std::endl;
    
    for(int i=0; i<2; ++i){
      for(int j=0; j<4*common::N_DSSD; ++j){
	if(i==0){
	  cADClow[i] ->cd(j+1)->Modified();
	  cADCdisc[i]->cd(j+1)->Modified();
	  cADChigh[i]->cd(j+1)->Modified();
	} //End if(i==0)
	
	cEall[i]->cd(j+1)->Modified();
	
      } //End j loop
      
      for(int j=0; j<(8*4); ++j){
	cADClow[1]->cd(j+1)->Modified();
      }
      
      cEall[i]->Update();
    }
    
    cADClow[0]  ->Update();
    cADClow[1]  ->Update();
    cADCdisc[0] ->Update();
    cADChigh[0] ->Update();
        
    for(int i=1; i<8; ++i){
      cTimeDist[0]->cd(i)->Modified();
    }
    cTimeDist[0]->Update();
    
    
    for(int i=0; i<common::N_DSSD*6; ++i) {
      if(i<common::N_DSSD*2) {
	//cEvtXY2_d  ->cd(i+1)->Modified();
      }
      if(i<common::N_DSSD*3) {
	std::cout << "Change to panel " << i+1 << std::endl;
	cEvtXY     ->cd(i+1)->Modified();
	std::cout << "cEvtXY" << std::endl;
	cEvtMulti  ->cd(i+1)->Modified();
	std::cout << "cEvtMulti" << std::endl;
	cEvtXY_d   ->cd(i+1)->Modified();
	std::cout << "cEvtXY_d" << std::endl;
	//cEvtMulti_d->cd(i+1)->Modified();
	std::cout << "cEvtMulti_d" << std::endl;
      }
      if(i<common::N_DSSD*4) {
	cEvtdXdY   ->cd(i+1)->Modified();
      }
      if(i<common::N_DSSD*5) {
	cEvtE_d    ->cd(i+1)->Modified();
      }
      if(i<common::N_DSSD*6) {
	cEvtE1     ->cd(i+1)->Modified();
      }
    }
    cEvtE1      ->Update();
    cEvtXY      ->Update();
    cEvtdXdY    ->Update();
    cEvtMulti   ->Update();
    cEvtE_d     ->Update();
    cEvtXY_d    ->Update();
    //cEvtXY2_d   ->Update();
    //cEvtMulti_d ->Update();
  }
}

/* Initialise analysis - load parameters from file, initialise histograms etc.
 *
 */
void Analysis::InitAnalysis(int opt, char *file_name){
  
  event_count = 0;
  t0_aida     = 0;
  t_low_prev  = 0;
  t_high_prev = 0;
  t_disc_prev = 0;

  b_pulser= false;

  ResetEvent();

  LoadParameters(file_name);

  std::string nombre[32];
  nombre[0]="NNAIDA1 (Det6, Pside)";
  nombre[1]="NNAIDA2 (Det5, Pside)";
  nombre[2]="NNAIDA3 (Det6, Nside)";
  nombre[3]="NNAIDA4 (Det5, Nside)";
  nombre[4]="NNAIDA5 (Det6, Pside)";
  nombre[5]="NNAIDA6 (Det5, Pside)";
  nombre[6]="NNAIDA7 (Det5, Nside)";
  nombre[7]="NNAIDA8 (Det6, Nside";
  nombre[8]="NNAIDA9 (Det4, Pside)";
  nombre[9]="NNAIDA10 (Det3, Pside)";
  nombre[10]="NNAIDA11 (Det4, Nside)";
  nombre[11]="NNAIDA12 (Det3, Nside)";
  nombre[12]="NNAIDA13 (Det4, Pside)";
  nombre[13]="NNAIDA14 (Det3, Pside)";
  nombre[14]="NNAIDA15 (Det3, Nside)";
  nombre[15]="NNAIDA16 (Det4, Nside)";
  nombre[16]="NNAIDA17 (Det1, Pside)";
  nombre[17]="NNAIDA18 (Det2, Pside)";
  nombre[18]="NNAIDA19 (Det2, Nside)";
  nombre[19]="NNAIDA20 (Det1, Nside)";
  nombre[20]="NNAIDA21 (Det2, Pside)";
  nombre[21]="NNAIDA22 (Det1, Pside)";
  nombre[22]="NNAIDA23 (Det1, Nside)";
  nombre[23]="NNAIDA24 (Det2, Nside)";
  nombre[24]="NNAIDA25 (DetX, Xside)";
  nombre[25]="NNAIDA26 (DetX, Xside)";
  nombre[26]="NNAIDA27 (DetX, Xside)";
  nombre[27]="NNAIDA28 (DetX, Xside)";
  nombre[28]="NNAIDA29 (DetX, Xside)";
  nombre[29]="NNAIDA30 (DetX, Xside)";
  nombre[30]="NNAIDA31 (DetX, Xside)";
  nombre[31]="NNAIDA32 (DetX, Xside)";
  
  // if first of 'opt' is 1, enable histogramming
  if( (opt & 0x01) == 1){
    char hname[256];
    char htitle[256];
    std::string stitle;
    std::string full_title;

    SetBHistograms(true);
    
    sprintf(hname,"cADCspec_lowE_2D");
    cADClow[0]= new TCanvas(hname, hname, 10,10,1200,900);  cADClow[0]->Divide(4,common::N_DSSD);
    sprintf(hname,"cADCspec_lowE_1D");
    cADClow[1]= new TCanvas(hname, hname, 10,10,1800,1000); cADClow[1]->Divide(8,common::N_DSSD);
    sprintf(hname,"cADCspec_highE_0");
    cADChigh[0]= new TCanvas(hname, hname, 20,20,1200,900); cADChigh[0]->Divide(4,common::N_DSSD);
    sprintf(hname,"cADCdisc_0");
    cADCdisc[0]= new TCanvas(hname, hname, 30,30,1200,900); cADCdisc[0]->Divide(4,common::N_DSSD);
    sprintf(hname,"cEall_0");
    cEall[0]= new TCanvas(hname, hname, 40,40,1200,900);    cEall[0]->Divide(4,common::N_DSSD);
    sprintf(hname,"cEall_1");
    cEall[1]= new TCanvas(hname, hname, 40,40,1200,900);    cEall[1]->Divide(4,common::N_DSSD);
    
    for(int i=0; i<common::N_FEE64; ++i){
      if(b_mod_enabled[i]){
	
	int k;
	k = 4*(geo_detector[i]-1)+geo_side[i]*2+geo_strip[i];    //Unique i.d. number for all FEE64s 1->common::N_FEE64(inc). Groups DSSDs together.
	int l;
	l = 8*(geo_detector[i]-1)+geo_side[i]*2+geo_strip[i];    //Unique i.d. number for all FEE64s 1->common::N_FEE64(inc). Groups high/low ADC range together w/ DSSD.
	int m;
	m = 8*(geo_detector[i]-1)+geo_side[i]*2+geo_strip[i]+4;  //Unique i.d. number for all FEE64s 1->common::N_FEE64(inc). Groups high/low ADC range together w/ DSSD.
	
	if(k<1 || k>common::N_FEE64 || l<1 || l>((common::N_DSSD*8)-4) || m<5 || m>(common::N_DSSD*8)) {
	  std::cout << " **** Analysis::InitAnalysis(): ERROR IN GEOMETRY"<<
	    "\n        mod, det, side, strip: " << i << " " << geo_detector[i]<< " " << geo_side[i]<< " " << "  " <<geo_strip[i]<< std::endl;
	}
	
	// *******************************************
	//                ADC singles
	// *******************************************
	
	//Low E range 2D ADC plots
	sprintf(hname,"hADClowFEE%i",i+1);
	sprintf(htitle,"2D per-strip ADC spectrum for lowE range;ch;Number of ADC data");
	stitle= htitle;
	full_title= nombre[i]+stitle;
	hADClowCh[i]= new TH2I(hname, full_title.data(), 64, 0, 64, 1024, 0, 65536);
	cADClow[0]->cd(k); hADClowCh[i]->Draw("colz"); gPad->SetLogz(1);
	
	//Low E range 1D ADC spectrum
	sprintf(hname,"hCh_ADClowFEE%i",i+1);
	sprintf(htitle,"Number of lowE ADC events per strip in FEE%i;ch;", i);
	stitle= htitle;
	full_title= nombre[i]+stitle;
	hCh_ADClow[i]= new TH1I(hname, full_title.data(), 64, 0, 64);
	cADClow[1]->cd(l); hCh_ADClow[i]->Draw("");

	//High E range 2D ADC plots
	sprintf(hname,"hADChighFEE%i",i+1);
	sprintf(htitle,"2D ADC plot for highE range;ch;Number of ADC data");
	stitle= htitle;
	full_title= nombre[i]+stitle;
	hADChighCh[i]= new TH2I(hname, full_title.data(), 64, 0, 64, 1024, 0, 65536);
	cADChigh[0]->cd(k); hADChighCh[i]->Draw("colz"); gPad->SetLogz(1);

	//High E range 1D ADC spectrum
	sprintf(hname,"hCh_ADChighFEE%i",i+1);
	sprintf(htitle,"Number of highE ADC events per strip in FEE%i;ch;", i);
	stitle= htitle;
	full_title= nombre[i]+stitle;
	hCh_ADChigh[i]= new TH1I(hname, full_title.data(), 64, 0, 64);
	cADClow[1]->cd(m); hCh_ADChigh[i]->Draw("");

	// ADC(decay range) ! discriminator hit
	sprintf(hname,"hADCdiscFEE%i",i+1);
	sprintf(htitle," ADC(low) !DISC hits;ch;ADC data");
	stitle= htitle;
	full_title= nombre[i]+stitle;
	hADCdiscCh[i]= new TH2I(hname, full_title.data(), 64, 0, 64, 1024, 0, 65536);
	cADCdisc[0]->cd(k); hADCdiscCh[i]->Draw("colz"); gPad->SetLogz(1);

	//Fast disc hits
	sprintf(hname,"hCh_ADCdiscFEE%i",i+1);
	sprintf(htitle,"Number of fast disc events per strip in FEE%i;ch;", i);
	stitle= htitle;
	full_title= nombre[i]+stitle;
	hCh_ADCdisc[i]= new TH1I(hname, full_title.data(), 64, 0, 64);
	
	//Low E range energy spectrum
	sprintf(hname,"hElowFEE%i",i+1);
	sprintf(htitle," E(low); E [MeV]");
	stitle= htitle;
	full_title= nombre[i]+stitle;
	if(geo_side[i]==0)       hElow[i]= new TH1I(hname, full_title.data(), 1000, 30000, 65000);
	else if(geo_side[i]==1)  hElow[i]= new TH1I(hname, full_title.data(), 1000, 0, 34000);
	cEall[0]->cd(k); hElow[i]->Draw(""); gPad->SetLogy(1);

	//Fast disc energy spectrum
	sprintf(hname,"hEdiscFEE%i",i+1);
	sprintf(htitle," E(disc); E [MeV]");
	stitle= htitle;
	full_title= nombre[i]+stitle;
	if(geo_side[i]==0)       hEdisc[i]= new TH1I(hname, full_title.data(), 1000, 30000, 65000);
	else if(geo_side[i]==1)  hEdisc[i]= new TH1I(hname, full_title.data(), 1000, 0, 34000);
	hEdisc[i]->SetLineColor(kRed);
	cEall[0]->cd(k); hEdisc[i]->Draw("same"); gPad->SetLogy(1);

	//High energy range enrgy spectrum
	sprintf(hname,"hEhighFEE%i",i+1);
	sprintf(htitle," E(high); E [GeV]");
	stitle= htitle;
	full_title= nombre[i]+stitle;
	if(geo_side[i]==0)       hEhigh[i]= new TH1I(hname, full_title.data(), 1000, 30000, 50000);
	else if(geo_side[i]==1)  hEhigh[i]= new TH1I(hname, full_title.data(), 1000, 15000, 34000);
	cEall[1]->cd(k); hEhigh[i]->Draw(""); gPad->SetLogy(1);


      } //End if(b_mod_enabled)
    } //End loop over common::N_FEE64
    
    //*************************************
    //        Timestamp monitoring
    //*************************************

    for(int i=0; i<common::N_DSSD; ++i){

      sprintf(hname,"cTimeDist_DSSD%i",i+1);
      cTimeDist[i] = new TCanvas(hname, hname, 40,40,1100,700); cTimeDist[i]->Divide(4,2);
     
      //time stamp distribution for ADClow
      sprintf(hname,"hTimeADClow_DSSD%i",i+1);
      sprintf(htitle,"Difference between timestamps in LEC (DSSD%i);t_n - t_(n-1) [clk ticks]",i+1);
      hTimeADClow[i] = new TH1I(hname, htitle, 1000, -500, 9500);

      cTimeDist[i]->cd(1); hTimeADClow[i]->Draw(""); gPad->SetLogy(1);

      //time stamp distribution for ADCdisc
      sprintf(hname,"hTimeADCdisc_DSSD%i",i+1);
      sprintf(htitle,"Difference between timestamps of fast discriminators (DSSD%i);t_n - t_(n-1) [clk ticks]",i+1);
      hTimeADCdisc[i] = new TH1I(hname, htitle, 1000, -500, 9500);
      cTimeDist[i]->cd(2); hTimeADCdisc[i]->Draw(""); gPad->SetLogy(1);
      
      //Time stamp distribution for ADChigh
      sprintf(hname,"hTimeADChigh_DSSD%i",i+1);
      sprintf(htitle,"#Difference between timestamps in HEC (DSSD%i);t_n - t_(n-1) [clk ticks]",i+1);
      hTimeADChigh[i] = new TH1I(hname, htitle, 1000, -500, 9500);
      cTimeDist[i]->cd(3); hTimeADChigh[i]->Draw(""); gPad->SetLogy(1);
    
    } //End loop over common::N_DSSD

    hTimeStamp = new TH1I("hTimeStamp","Time stamp;tmstp - t_first [1/1e6]",1000,0,1e5);
    hTimeStampExt = new TH1I("hTimeStampExt","Time stamp External;tmstp [1/1e6]",1000,0,1e5);
    hTimeStampFlag = new TH1I("hTimeStampFlag","Time stamp Ext Flag;corr_flag",2,0,2);
    cTimeDist[0]->cd(4);  hTimeStampFlag->Draw(""); gPad->SetLogy(1);

    hEvt_TmStpDist = new TH1I("hEvt_TmStpDist0","Time Stamp Dist (within event);ts - ts_{0} [arb. units]",500,-2500,4500);

    //sprintf(hname,"hEvt_Eaida");
    //sprintf(htitle,"Energy AIDA (DSSD sum);Energy [ch]");
    //hEvt_Eaida= new TH1I(hname, htitle, 1024, 0, 49152);

    // *********************************************************
    //       Histograms to monitor event reconstruction
    // *********************************************************

    std::cout << "     InitAnalysis(): initialize histograms for event clustering"<<std::endl;

    //Create canvases for implant event plots
    cEvtE1 = new TCanvas("cEvtE1","cEvt ImplantEnergyPlots", 100,100,1200,1200); cEvtE1->Divide(common::N_DSSD,6);
    cEvtXY = new TCanvas("cEvtXY","cEvt ImplantXY", 140,140,900,1000); cEvtXY->Divide(common::N_DSSD,3);
    cEvtdXdY = new TCanvas("cEvtdXdY","cEvt Implant_dX:dY", 140,140,1200,800); cEvtdXdY->Divide(common::N_DSSD,4);
    cEvtMulti = new TCanvas("cEvtMulti","cEvt Multiplicity", 140,140,1200,1000); cEvtMulti->Divide(common::N_DSSD,3);

    //Create canvases for decay event plots
    cEvtE_d = new TCanvas("cEvtE_d","cEvt Energy1 (DECAY)", 100,100,1200,1000); cEvtE_d->Divide(common::N_DSSD,5);
    cEvtXY_d = new TCanvas("cEvtXY_d","cEvt DecayXY", 140,140,1200,900); cEvtXY_d->Divide(common::N_DSSD,3);
    cEvtXY2_d = new TCanvas("cEvtXY2_d","Multiplicity vs Cluster size (DECAY)", 140,140,1000,800); cEvtXY2_d->Divide(common::N_DSSD,2);
    cEvtMulti_d = new TCanvas("cEvtMulti_d","cEvt Multiplicity (DECAY)", 140,140,1200,800); cEvtMulti_d->Divide(common::N_DSSD,3);

    for(int i=0; i<common::N_DSSD; ++i) {

      //*************************************
      //        IMPLANTS IMPLANTS
      //*************************************

      //x-y plot of implant event distribution
      sprintf(hname,"hEvt_XYi_DSSD%i",i+1);
      sprintf(htitle,"x-y plot of implant event distribution in DSSD%i;X [ch];Y [ch]",i+1);
      hEvt_XY_i[i] = new TH2I(hname, htitle, 32, 0, 128, 32, 0, 128);
      cEvtXY->cd(i+(2*common::N_DSSD)+1); hEvt_XY_i[i]->Draw("colz"); gPad->SetLogz(0);
      
      //x-strip distribution of implant events
      sprintf(hname,"hEvt_Xi_DSSD%i",i+1);
      sprintf(htitle,"x-strip distribution of implant events in DSSD%i;X [ch]",i+1);
      hEvt_X_i[i] = new TH1I(hname, htitle, 128, 0, 128);
      cEvtXY->cd(i+1); hEvt_X_i[i]->Draw(""); gPad->SetLogy(1);
      
      //y-strip distribution of implant events
      sprintf(hname,"hEvt_Yi_DSSD%i",i+1);
      sprintf(htitle,"y-strip distribution of implant events in DSSD%i;Y [ch]",i+1);
      hEvt_Y_i[i] = new TH1I(hname, htitle, 128, 0, 128);
      cEvtXY->cd(i+common::N_DSSD+1); hEvt_Y_i[i]->Draw(""); gPad->SetLogy(1);
      
      //E_x vs E_y for implant events
      sprintf(hname,"hEvt_ExEy_i_all%i",i+1);
      sprintf(htitle,"Ex vs Ey for unpaired implant cluster events (DSSD%i);Ex [ch];Ey [ch]",i+1);
      hEvt_ExEy_i[i][0] = new TH2I(hname, htitle, 1024, 0, 16384, 1024, 0, 16384);
 
      sprintf(hname,"hEvt_ExEy_i_clust%i",i+1);
      sprintf(htitle,"Ex vs Ey for ExEy matched implant cluster events (DSSD%i);Ex [ch];Ey [ch]",i+1);
      hEvt_ExEy_i[i][1] = new TH2I(hname, htitle, 1024, 0, 16384, 1024, 0, 16384);
      cEvtE1->cd(i+(2*common::N_DSSD)+1); hEvt_ExEy_i[i][1]->Draw("colz"); gPad->SetLogz(0);
 
      //Energy distribution of n-side implant events
      sprintf(hname,"hEvt_Eside_i%i_n",i+1);
      sprintf(htitle,"Energy of n-side implant events (DSSD%i);Energy [ch]",i+1);
      hEvt_Eside_i[i][0] = new TH1I(hname, htitle, 512, 0, 32768);
      cEvtE1->cd(i+1); hEvt_Eside_i[i][0]->Draw(""); gPad->SetLogy(1);

      //Energy distribution of p-side implant events
      sprintf(hname,"hEvt_Eside_i%i_p",i+1);
      sprintf(htitle,"Energy of p-side implant events (DSSD%i);Energy [ch]",i+1);
      hEvt_Eside_i[i][1] = new TH1I(hname, htitle, 512, 0, 32768);
      cEvtE1->cd(i+common::N_DSSD+1); hEvt_Eside_i[i][1]->Draw(""); gPad->SetLogy(1);
      
      //Multiplicity of implant events
      sprintf(hname,"hEvt_Mult_i_n%i",i+1);
      sprintf(htitle,"Multiplicty of implant events (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_Mult_i[i][0] = new TH1I(hname, htitle, 128, 0, 128);
      sprintf(hname,"hEvt_Mult_i_p%i",i+1);
      sprintf(htitle,"Multiplicty of implant events (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_Mult_i[i][1] = new TH1I(hname, htitle, 128, 0, 128);

      //x_mult vs y_mult for implant events
      sprintf(hname,"hEvt_Mult_i_all%i",i+1);
      sprintf(htitle,"Multiplicty of all possible implant events (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_MultXY_i[i][0] = new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);
      sprintf(hname,"hEvt_Mult_i_pair %i",i+1);
      sprintf(htitle,"Multiplicty of n/p paired implant events (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_MultXY_i[i][1] = new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);
 
      //*************************************
      //     Heavy ion implant tracking
      //*************************************
     
      /*if(i<common::N_DSSD-1) {
	
	//Change in x/y-position between planes
	//sprintf(hname,"hEvt_dXi_%i-%i",i+2,i+1);
	//sprintf(htitle,"dX of fast ions between DSSD%i and DSSD%i;X%i - X%i [ch]",i+2,i+1,i+2,i+1);
	//hEvt_dX[i] = new TH1I(hname, htitle, 64, -128, 128);
	//cEvtdXdY->cd(i+1); hEvt_dX[i]->Draw(""); gPad->SetLogy(1);
	
	sprintf(hname,"hEvt_dYi_%i-%i",i+2,i+1);
	sprintf(htitle,"dY of fast ions between DSSD%i and DSSD%i;Y%i - Y%i [ch]",i+2,i+1,i+2,i+1);
	hEvt_dY[i] = new TH1I(hname, htitle, 64, -128, 128);
	cEvtdXdY->cd(i+1+common::N_DSSD); hEvt_dY[i]->Draw(""); gPad->SetLogy(1);
	}
      
      else if(i==common::N_DSSD-1) {
	sprintf(hname,"hEvt_dXi_FirstLast");
	sprintf(htitle,"Fast ion dX between the first and last DSSDs;X_first - X_last [ch]");
	hEvt_dX[i] = new TH1I(hname, htitle, 64, -128, 128);
	cEvtdXdY->cd(i+1); hEvt_dX[i]->Draw(""); gPad->SetLogy(1);

	sprintf(hname,"hEvt_dYi_FirstLast");
	sprintf(htitle,"Fast ion dY between the first and last DSSDs;Y_first - Y_last [ch]");
	hEvt_dY[i] = new TH1I(hname, htitle, 64, -128, 128);
	cEvtdXdY->cd(i+1+common::N_DSSD); hEvt_dX[i]->Draw(""); gPad->SetLogy(1);
	}
      
      if(i<common::N_DSSD-2) {

	//Change in x/y-position between DSSDi+1/DSSDi vs change between DSSDi+2/DSSDi+1
	sprintf(hname,"hEvt_dXdXi_%i-%i-%i",i+3,i+2,i+1);
	sprintf(htitle,"Fast ion dX between DSSD%i/DSSD%i vs DSSD%i/DSSD%i;X%i - X%i [ch];X%i - X%i [ch]",i+2,i+1,i+3,i+2,i+2,i+1,i+3,i+2);
	hEvt_dXdX[i] = new TH2I(hname, htitle, 32, -128, 128, 32, -128, 128);
	cEvtdXdY->cd(i+1+(2*common::N_DSSD)); hEvt_dXdX[i]->Draw("colz"); gPad->SetLogz(1);

	sprintf(hname,"hEvt_dYdYi_%i-%i-%i",i+3,i+2,i+1);
	sprintf(htitle,"Fast ion dY between DSSD%i/DSSD%i vs DSSD%i/DSSD%i;X%i - X%i [ch];X%i - X%i [ch]",i+2,i+1,i+3,i+2,i+2,i+1,i+3,i+2);
	hEvt_dYdY[i] = new TH2I(hname, htitle, 32, -128, 128, 32, -128, 128);
	cEvtdXdY->cd(i+1+(3*common::N_DSSD)); hEvt_dYdY[i]->Draw("colz"); gPad->SetLogz(1);
      }

      else if(i==common::N_DSSD-2) {
	
	//Change in x-y position between first/middle/last planes
	sprintf(hname,"hEvt_dXdXi_%i-%i-%i",common::N_DSSD,(common::N_DSSD/2)+1,1);
	sprintf(htitle,"Fast ion dX between DSSD%i/DSSD%i vs DSSD%i/DSSD%i;X%i - X%i [ch];X%i - X%i [ch]",
		(common::N_DSSD/2)+1,1,common::N_DSSD,(common::N_DSSD/2)+1,(common::N_DSSD/2)+1,1,common::N_DSSD,(common::N_DSSD/2)+1);
	hEvt_dXdX[i] = new TH2I(hname, htitle, 32, -128, 128, 32, -128, 128);
	cEvtdXdY->cd(i+1+(2*common::N_DSSD)); hEvt_dXdX[i]->Draw("colz"); gPad->SetLogz(1);
	
	sprintf(hname,"hEvt_dYdYi_%i-%i-%i",common::N_DSSD,(common::N_DSSD/2)+1,1);
	sprintf(htitle,"Fast ion dY between DSSD%i/DSSD%i vs DSSD%i/DSSD%i;X%i - X%i [ch];X%i - X%i [ch]",
		(common::N_DSSD/2)+1,1,common::N_DSSD,(common::N_DSSD/2)+1,(common::N_DSSD/2)+1,1,common::N_DSSD,(common::N_DSSD/2)+1);
	hEvt_dYdY[i] = new TH2I(hname, htitle, 32, -128, 128, 32, -128, 128);
	cEvtdXdY->cd(i+1+(3*common::N_DSSD)); hEvt_dYdY[i]->Draw("colz"); gPad->SetLogz(1);
	
	}

      //Strip multiplicity of n- and p- side implant events
      sprintf(hname,"hEvt_Multi_DSSD%i_n",i+1);
      sprintf(htitle,"Multiplicity of n-side implant events (DSSD%i);N hits (n-side)",i+1);
      hEvt_Multi[i][0] = new TH1I(hname, htitle, 15, 0, 15);
      cEvtMulti->cd(i+1); hEvt_Multi[i][0]->Draw(""); gPad->SetLogy(1);
      
      sprintf(hname,"hEvt_Multi_DSSD%i_p",i+1);
      sprintf(htitle,"Multiplicity of p-side events (DSSD%i);N hits (p-side)",i+1);
      hEvt_Multi[i][1] = new TH1I(hname, htitle, 15, 0, 15);
      cEvtMulti->cd(i+1+common::N_DSSD); hEvt_Multi[i][1]->Draw(""); gPad->SetLogy(1);
      */

      //************************************
      //          DECAYS DECAYS
      //************************************
      
      //x-y distribution of decay evnets
      sprintf(hname,"hEvt_XY_d%i",i+1);
      sprintf(htitle,"XY distribution of decay events (DSSD%i);x;y",i+1);
      hEvt_XY_d[i]= new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);
      cEvtXY_d->cd(i+1);  hEvt_XY_d[i]->Draw("colz"); gPad->SetLogz(1);

      //x distribution of decay events
      sprintf(hname,"hEvt_X_d%i",i+1);
      sprintf(htitle,"X distribution of decay events (DSSD%i);x;Counts",i+1);
      hEvt_X_d[i]= new TH1I(hname, htitle, 128, 0, 128);

      //y distribution of decay events
      sprintf(hname,"hEvt_Y_d%i",i+1);
      sprintf(htitle,"Y distribution of decay events (DSSD%i);x;y",i+1);
      hEvt_Y_d[i]= new TH1I(hname, htitle, 128, 0, 128);

      //E_x vs E_y for decay events
      sprintf(hname,"hEvt_ExEy_d_all%i",i+1);
      sprintf(htitle,"Ex vs Ey for unpaired decay cluster events (DSSD%i);Ex [keV];Ey [keV]",i+1);
      hEvt_ExEy_d[i][0] = new TH2I(hname, htitle, 1024, 0, 16384, 1024, 0, 16384);

      sprintf(hname,"hEvt_ExEy_d_clust%i",i+1);
      sprintf(htitle,"Ex vs Ey for ExEy matched decay cluster events (DSSD%i);Ex [keV];Ey [kEV]",i+1);
      hEvt_ExEy_d[i][1] = new TH2I(hname, htitle, 1024, 0, 16384, 1024, 0, 16384);
      cEvtE_d->cd(i+1+(2*common::N_DSSD));  hEvt_ExEy_d[i][1]->Draw("colz"); gPad->SetLogz(1); 
     
      // Energy distribution of decay events
      sprintf(hname,"hEvt_Eside_d%i_n",i+1);
      sprintf(htitle,"Energy distribution of n-side decay events (DSSD%i);E [***]",i+1);
      hEvt_Eside_d[i][0] = new TH1I(hname, htitle, 512, 0, 32768);

      sprintf(hname,"hEvt_Eside_d%i_p",i+1);
      sprintf(htitle,"Energy distribution of p-side decay events (DSSD%i);E [***]",i+1);
      hEvt_Eside_d[i][1] = new TH1I(hname, htitle, 512, 0, 32768);

      cEvtE_d->cd(i+1);  hEvt_Eside_d[i][0]->Draw(); gPad->SetLogy(1);
      cEvtE_d->cd(i+1+common::N_DSSD);  hEvt_Eside_d[i][1]->Draw(); gPad->SetLogy(1);

      //Multiplicity of decay events
      sprintf(hname,"hEvt_Mult_d_n%i",i+1);
      sprintf(htitle,"Multiplicty of n-side decay events (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_Mult_d[i][0] = new TH1I(hname, htitle, 128, 0, 128);
      sprintf(hname,"hEvt_Mult_d_p%i",i+1);
      sprintf(htitle,"Multiplicty of p-side decay events (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_Mult_d[i][1] = new TH1I(hname, htitle, 128, 0, 128);

      //x_mult vs y_mult for decay events
      sprintf(hname,"hEvt_Mult_d_all%i",i+1);
      sprintf(htitle,"Multiplicty of all possible decay event combos (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_MultXY_d[i][0] = new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);
      sprintf(hname,"hEvt_Mult_d_pair %i",i+1);
      sprintf(htitle,"Multiplicty of n/p paired decay event combos (DSSD%i); Multilpicity [ch];Counts",i+1);
      hEvt_MultXY_d[i][1] = new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);      
      
    } //End loop over common::N_DSSD

    /*for(int i=0; i<common::N_DSSD; ++i){

      // ************************************************
      //               Decay event historgams
      // ************************************************

      //E_x vs E_y for decay events
      sprintf(hname,"hEvt_ExEy_d%i",i+1);
      sprintf(htitle,"Ex vs Ey for decay events (DSSD%i);Ex [***];Ey [***]",i+1);
      hEvt_ExEy_d[i]= new TH2I(hname, htitle, 256, 0, 16384, 256, 0, 16384);
      cEvtE_d->cd(i+1+(2*common::N_DSSD));  hEvt_ExEy_d[i]->Draw("colz"); gPad->SetLogz(1);
      
      //Multiplicity of decay events by n- or p-side
      sprintf(hname,"hEvt_MultiStrip_d%i_n",i+1);
      sprintf(htitle,"Strip multiplicity of n-side decay events (DSSD%i);N strips",i+1);
      hEvt_MultiStrip_d[i][0] = new TH1I(hname,htitle,40,0,40);
      cEvtMulti_d->cd(i+1+common::N_DSSD);  hEvt_MultiStrip_d[i][0]->Draw(""); gPad->SetLogy(1);

      sprintf(hname,"hEvt_MultiStrip_d%i_p",i+1);
      sprintf(htitle,"Strip multiplicity of p-side decay events (DSSD%i);N strips",i+1);
      hEvt_MultiStrip_d[i][1] = new TH1I(hname,htitle,40,0,40);
      cEvtMulti_d->cd(i+1+(2*common::N_DSSD));  hEvt_MultiStrip_d[i][1]->Draw(""); gPad->SetLogy(1);

      //Energy spectrum of decay events
      sprintf(hname,"hEvt_Eside_df%i_n",i+1);
      sprintf(htitle,"Energy spectrum of n-side decay events (DSSD%i);E [***]",i+1);
      hEvt_Eside_df[i][0]= new TH1I(hname, htitle, 512, 0, 16384);
      hEvt_Eside_df[i][0]->SetLineColor(kRed);

      sprintf(hname,"hEvt_Eside_df%i_p",i+1);
      sprintf(htitle,"Energy spectrum of p-side decay events (DSSD%i);E [***]",i+1);
      hEvt_Eside_df[i][1]= new TH1I(hname, htitle, 512, 0, 16384);
      hEvt_Eside_df[i][1]->SetLineColor(kRed);

      //Energy spectrum of decay events by side
      sprintf(hname,"hEvt_Eside_df2%i_n",i+1);
      sprintf(htitle,"E(decay) Det%i n-side !decay2;E [***]",i+1);
      hEvt_Eside_df2[i][0]= new TH1I(hname, htitle, 512, 0, 16384);
      hEvt_Eside_df2[i][0]->SetLineColor(kGreen);

      sprintf(hname,"hEvt_Eside_df2%i_p",i+1);
      sprintf(htitle,"E(decay) Det%i p-side !decay2;E [***]",i+1);
      hEvt_Eside_df2[i][1]= new TH1I(hname, htitle, 512, 0, 16384);
      hEvt_Eside_df2[i][1]->SetLineColor(kGreen);

      //E_x vs E_y for decay, w/ extra conditions
      sprintf(hname,"hEvt_ExEy_df%i",i);
      sprintf(htitle,"Ex vs Ey for decay events (DSSD%i, decay_flag>0 && decay_flag%%10==%i);Ex [***];Ey [***]",i+1,i+1);
      hEvt_ExEy_df[i]= new TH2I(hname, htitle, 256, 0, 16384, 256, 0, 16384);

      sprintf(hname,"hEvt_ExEy_df2%i",i);
      sprintf(htitle,"Ex vs Ey for decay events (DSSD%i, decay_flag>10 && decay_flag%%10==%i);Ex [***];Ey [***]",i+1,i+1);
      hEvt_ExEy_df2[i]= new TH2I(hname, htitle, 256, 0, 16384, 256, 0, 16384);

      cEvtE_d->cd(i+1);                 hEvt_Eside_df[i][0]->SetLineColor(kRed); hEvt_Eside_df[i][0]->Draw("same");// gPad->SetLogy(1);
      cEvtE_d->cd(i+1+common::N_DSSD);  hEvt_Eside_df[i][1]->SetLineColor(kRed); hEvt_Eside_df[i][1]->Draw("same");// gPad->SetLogy(1);

      cEvtE_d->cd(i+1);                 hEvt_Eside_df2[i][0]->SetLineColor(kGreen); hEvt_Eside_df2[i][0]->Draw("same");// gPad->SetLogy(1);
      cEvtE_d->cd(i+1+common::N_DSSD);  hEvt_Eside_df2[i][1]->SetLineColor(kGreen); hEvt_Eside_df2[i][1]->Draw("same");// gPad->SetLogy(1);

      cEvtE_d->cd(i+1+(3*common::N_DSSD));   hEvt_ExEy_df[i]->Draw("colz"); gPad->SetLogz(1);
      cEvtE_d->cd(i+1+(4*common::N_DSSD));   hEvt_ExEy_df2[i]->Draw("colz"); gPad->SetLogz(1);
	
      //XY plots of decay events with conditions
      sprintf(hname,"hEvt_XY_df%i",i);
      sprintf(htitle,"X vs Y (decay) Det%i !decay;Ex [***];Ey [***]",i);
      hEvt_XY_df[i]= new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);

      sprintf(hname,"hEvt_XY_df2%i",i);
      sprintf(htitle,"X vs Y (decay) Det%i !decay2;Ex [***];Ey [***]",i);
      hEvt_XY_df2[i]= new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);

      cEvtXY_d->cd(i+1+common::N_DSSD);      hEvt_XY_df[i]->Draw("colz"); gPad->SetLogz(0);
      cEvtXY_d->cd(i+1+(2*common::N_DSSD));  hEvt_XY_df2[i]->Draw("colz"); gPad->SetLogz(0);
      
      //Strip multiplicty vs strip_max - strip_min by side (decay events)
      sprintf(hname,"hEvt_MultidX_d%i_n",i+1);
      sprintf(htitle,"Multiplicity vs strip_max-strip_min for n-side decays (DSSD%i);hit multiplicity;x_max - x_min [ch]",i+1);
      hEvt_MultidX_d[i][0]= new TH2I(hname, htitle, 32, 0, 32, 64, 0, 64);
      cEvtXY2_d->cd(i+1);  hEvt_MultidX_d[i][0]->Draw("colz"); gPad->SetLogz(1);

      sprintf(hname,"hEvt_MultidX_d%i_p",i+1);
      sprintf(htitle,"Multiplicity vs strip_max-strip_min for p-side decays (DSSD%i), p-side;hit multiplicity;x_max - x_min [ch]",i+1);
      hEvt_MultidX_d[i][1]= new TH2I(hname, htitle, 32, 0, 32, 64, 0, 64);
      cEvtXY2_d->cd(i+1+common::N_DSSD);  hEvt_MultidX_d[i][1]->Draw("colz"); gPad->SetLogz(1);


      // ******************************************************
      //                Implant event histograms
      // ******************************************************
      
      //Energy spectrum of implant events by side
      sprintf(hname,"hEvt_Eside_if%i_n",i+1);
      sprintf(htitle,"Energy spectrum of n-side implant events (DSSD%i);Energy [ch]",i+1);
      hEvt_Eside_if[i][0]= new TH1I(hname, htitle, 512, 0, 16384);
      hEvt_Eside_if[i][0]->SetLineColor(kRed);
      cEvtE1->cd(i+1); hEvt_Eside_if[i][0]->Draw("same"); gPad->SetLogy(1);
      
      sprintf(hname,"hEvt_Eside_if%i_p",i+1);
      sprintf(htitle,"Energy spectrum of p-side implant events (DSSD%i);Energy [ch]",i+1);
      hEvt_Eside_if[i][1]= new TH1I(hname, htitle, 512, 0, 16384);
      hEvt_Eside_if[i][1]->SetLineColor(kRed);
      cEvtE1->cd(i+5); hEvt_Eside_if[i][1]->Draw("same"); gPad->SetLogy(1);
      
      //E_x vs E_y for implant events
      sprintf(hname,"hEvt_ExEy_if%i",i+1);
      sprintf(htitle,"E_x vs E_y for implant events (DSSD%i);Energy n-side [ch];Energy p-side [ch]",i+1);
      hEvt_ExEy_if[i]= new TH2I(hname, htitle, 256, 0, 16384, 256, 0, 16384);
      cEvtE1->cd(i+13); hEvt_ExEy_if[i]->Draw("colz"); gPad->SetLogz(0);

      //XY plot of implant events
      sprintf(hname,"hEvt_XY_if%i",i+1);
      sprintf(htitle,"X vs Y distribution of implant events (DSSD%i);Ex [***];Ey [***]",i+1);
      hEvt_XY_if[i]= new TH2I(hname, htitle, 128, 0, 128, 128, 0, 128);

      //Composite ADC spectrum combinining all calibrated strips
      sprintf(hname,"hADClow_all_DSSD%i",i+1);
      sprintf(htitle,"Composite ADC spectrum for whole DSSD%i - no cuts;Channel;Counts",i+1);
      hADClow_all[i][0] = new TH1I(hname, htitle, 4096, 0, 65536);
      
      sprintf(hname,"hADClow_all_DSSD%if",i+1);
      sprintf(htitle,"Composite ADC spectrum for whole DSSD%i - 1st cut;Channel;Counts",i+1);
      hADClow_all[i][1] = new TH1I(hname, htitle, 4096, 0, 65536);
      
      sprintf(hname,"hADClow_all_DSSD%if2",i+1);
      sprintf(htitle,"Composite ADC spectrum for whole DSSD%i - 2nd cut;Channel;Counts",i+1);
      hADClow_all[i][2] = new TH1I(hname, htitle, 4096, 0, 65536);
      
    }*/

    //Low energy multiplicity during implant events
    hEvt_Mult_impdec = new TH1I("hEvt_Mult_impdec", "Low energy multiplicity during implant events (DSSD%i); Multilpicity [ch];Counts", 128, 0, 128);

    //Energy difference Ex - Ey plots for implant and decay pairs
    hEvt_residualE_i = new TH1I("hEvt_residualE_i","Ex-Ey for unpaired implant clusters; Ex - Ey [MeV];Counts",1024,-2000,2000);
    hEvt_residualE_d = new TH1I("hEvt_residualE_d","Ex-Ey for unpaired decay clusters; Ex - Ey [keV];Counts",1024,-2000,2000);
    
    hEvt_pulserMult = new TH2I("hEvt_pulserMult","Multiplicity of pulser events; x_mult [ch]; y_mult [ch]", 128, 0, 128, 128, 0 ,128);

    //cEvtHits_d= new TCanvas("cEvtHits_d","cEvt Hit Patterns (DECAY)", 150,150,1200,800); cEvtHits_d->Divide(2,2);
    //cEvtE_d->cd(3);  hEvt_EPulser_d->Draw("colz"); gPad->SetLogz(0);
    /*    
    //Tally of implant events per DSSD side
    hEvt_HitsSide = new TH1I("hEvt_HitsSide","Detector implant events (each side);DSSD# + side(n=0, p=1)*0.5;n_hits",2*common::N_DSSD,1,common::N_DSSD+1);
    cEvtMulti->cd((2*common::N_DSSD)+1); hEvt_HitsSide->Draw(""); gPad->SetLogy(1);
    
    //Tally of implant per DSSD
    hEvt_HitsDet = new TH1I("hEvt_HitsDet","Detector implant events; DSSD#; n_hits",common::N_DSSD,1,common::N_DSSD+1);
    cEvtMulti->cd((2*common::N_DSSD)+2); hEvt_HitsDet->Draw(""); gPad->SetLogy(1);

    //E_x vs E_y for pulser events
    hEvt_EPulser_d = new TH2I("hEvt_EPulser_E","Pulser spectra;<E n-side> [ch];<E p-side> [ch]",100,0,32000, 100, 0, 32000); 
    
    hEvt_TmStpDist[0] = new TH1I("hEvt_TmStpDist0","Time Stamp Dist (within event);ts - ts_{0} [arb. units]",500,-2500,4500); 
    cTimeDist[0]->cd(5);

    //cEvtTS->cd(1); 
    hEvt_TmStpDist[0]->Draw(); gPad->SetLogy(1);

    hEvt_TmStpDist[1] = new TH1I("hEvt_TmStpDist1","Time Stamp Dist (last hit);#Delta ts [arb. units]",500,-2500,4500); 
    cTimeDist[0]->cd(6);
    //    cEvtTS->cd(2); 
    hEvt_TmStpDist[1]->Draw(); gPad->SetLogy(1);

    hEvt_TmStpDist[2] = new TH1I("hEvt_TmStpDist2","Time Stamp Dist (within event !DISC);#Delta ts [arb. units]",500,-2500,4500); 
    cTimeDist[0]->cd(7);
    //    cEvtTS->cd(3); 
    hEvt_TmStpDist[2]->Draw(); gPad->SetLogy(1);

    hEvt_TmStpDist[3] = new TH1I("hEvt_TmStpDist3","Time Stamp Dist (within event !ADC);#Delta ts [arb. units]",500,-2500,4500); 
    cTimeDist[0]->cd(8);
    //    cEvtTS->cd(4); 
    hEvt_TmStpDist[3]->Draw(); gPad->SetLogy(1);          
    //x_mult vs y_mult for implant events
    hEvt_Mult_d = new TH1I("hEvt_Mult_d_all", "Multiplicty of decay events; Multilpicity [ch];Counts", 128, 0, 128);

    //Histogram to show how many decay events satisfy the various conditions
    //0 = n_decay>0, 1 = 0 && decay_flag>0 && flag%10==det# 1, 2 = 0 && 1 && decay+flag>10
    hEvt_flag_d = new TH1I("hEvt_flag_d","Number of decays satifying each condition", 5, 0, 5);

    //Number of decay events by detector and side
    hEvt_MultiDet_d = new TH1I("hEvt_MultiDet_d","Tally of decay events per DSSD;DSSD#;Events#",common::N_DSSD,1,common::N_DSSD+1);
    hEvt_MultiSide_d = new TH2I("hEvt_MultiSide_d","Tally of decay events by DSSD/side;Side (0 = n, 1 = p, 2 = both);DSSD#",3,0,3,common::N_DSSD,1,common::N_DSSD+1);

    //Multiplicity of implant events vs decay events
    hEvt_MultiID= new TH2I("hEvt_MultiID","Multiplicity of Implant vs Decay events;N implant hits;N decay hits",40,0,40,40,0,40);

    //Implant and decay flags
    hEvt_HitsFlag= new TH1I("hEvt_HitsFlag","Number of Implant & Decay Flags;Flag (0: implant, 1: decay)",2,0,2);

    cEvtMulti_d->cd(1);  hEvt_MultiID->Draw("colz"); gPad->SetLogz(1);
    cEvtMulti_d->cd(2);  hEvt_MultiDet_d->Draw(); gPad->SetLogy(1);
    cEvtMulti_d->cd(3);  hEvt_MultiSide_d->Draw("colz"); gPad->SetLogz(1);
    cEvtMulti_d->cd(4);  hEvt_HitsFlag->Draw(""); gPad->SetLogy(1);
    */  
    if(b_debug) std::cout << "db    Analysis.cpp: Initialized histograms and canvas for this step"<<std::endl;
    
  }

  //if second bit of 'opt' is 1, enable TTree for output
  if( ( (opt & 0x02) >> 1) == 1){
  
    std::cout << " ***     Analysis::InitAnalysis(): initializing TTree" << std::endl;
    //Initialize TTree
    out_root_tree = new TTree("AIDA_hits","AIDA_hits");
    out_root_tree->Branch("aida_hit",&hit,"T/l:Tfast/l:E/D:EX/D:EY/D:x/D:y/D:z/D:nx/I:ny/I:nz/I:ID/b");
    SetBRootTree(true);
  }

}

/* Load system parameters from text file.
 *
 */
void Analysis::LoadParameters(char *file_name) {

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

	++line_count;    //Add to line count to keep forever loop in check

	if(GetBDebug()) std::cout << " * LoadParam() - line: " << my_line << std::endl;

	if(line_count < MAX_LINES) {

	  //skip empty lines and comments
	  if(!my_line.empty() && '\n' != my_line[0] && '#' != my_line[0]){
	
	    std::istringstream iss(my_line);  
	    y = -1;
	    x = -1;
	    z= -1;

	    iss >> par_name >> x >> y >> z;

	    if(GetBDebug()) std::cout << " * LoadParameter() - values read: " << par_name << " "<< x << y << z << std::endl;

	    if(par_name=="b_mod_enabled"){
	      module= x; data = y;
	      if(IsValidChannel(module, 0)){
		if(y == 0) {
		  b_mod_enabled[module-1] = false;
		  ++par_count;
		}
		else if (y == 1) {
		  b_mod_enabled[module-1]= true;
		  ++par_count;
		}
		else if(y != 0 && y != 1) {
		  std::cout << " *** Incorrect input parameter for b_mod_enabled[" << module << "] = " << y << ". y should be = 1 or 0." << std::endl;
		}
	      }
	      else if(!IsValidChannel(module, 0)) {
		std::cout << "Invalid channel for b_mod_enabled: mod " << module << " does not exist." << std::endl;
		std::cout << "Module enabled by default!!" << std::endl;	      
	      }
	    }
	    else if(par_name=="map_dssd"){
	      module= x; data= y;
	      if(IsValidChannel(module, 0)){
		++par_count;
		geo_detector[module-1]= int(data);
	      }
	      else if(!IsValidChannel(module, 0)) {
		std::cout << "Invalid channel for map_dssd: mod " << module << " does not exist." << std::endl;
		std::cout << "Module enabled by default!!" << std::endl;	      
	      }
	    }
	    else if(par_name=="map_side"){
	      module= x; data= y;
	      if(IsValidChannel(module, 0)){
		++par_count;
		geo_side[module-1]= int(data);
	      }
	      else if(!IsValidChannel(module, 0)) {
		std::cout << "Invalid channel for map_side: mod " << module << " does not exist." << std::endl;
		std::cout << "Module enabled by default!!" << std::endl;
	      }
	    }
	    else if(par_name=="map_strip"){
	      module= x; data= y;
	      if(IsValidChannel(module, 0)){
		++par_count;
		geo_strip[module-1]= int(data);
	      }
	      else if(!IsValidChannel(module, 0)) {
		std::cout << "Invalid channel for map_strip: mod " << module << " does not exist." << std::endl;
		std::cout << "Module enabled by default!!" << std::endl;
	      }
	    }

	    else{
	      if(b_debug) std::cout << " Analysis.cpp:: parameter name not recognized: " << par_name << "." <<std::endl;
	    }

	  } //End if empty line
	} //End if line_count< MAX_LINES
	else {
	  std::cout << " ***** Analysis.cpp::LoadParameters() ***** " << std::endl 
		    << "That is an awful lot of prameters you're trying to load there buddy... This param file is over " << MAX_LINES << " lines long..." << std::endl
		    << par_count << " parameters loaded, now breaking out of the loop and continuing..." << std::endl;
	  break;
	}
      } //End if myLine
      
      else{
	std::cout << " Analysis.cpp:: LoadParameters() reached end of parameter file. Loaded new value for "<< par_count << " parameters." << std::endl;
	break;
      }
      
    } //End for loop
    
  } //End if param_file.open()
  
  else{
    std::cout << " Analysis.cpp:: Error opening parameter file: " << file_name << " *****************"  << std::endl << std::endl;
  }
  return;
}

/* Reset internal variables.
 *
 */
void Analysis::ResetEvent(){

  b_decay   = false;
  b_pulser  = false;
  b_implant = false;

  total_mult[0] = 0;
  total_mult[1] = 0;
  total_e_dep[0] = 0;
  total_e_dep[1] = 0;

  decay_evts.clear();
  implant_evts.clear();

  evt_data.multiplicity= 0;
  
  t_aida_prev = -999;
  
  /*hit.T     = -1;
  hit.t_ext = -1;
  hit.x     = -1;
  hit.y     = -1;
  hit.z     = -1;
  hit.ex    = -1;
  hit.ey    = -1;
  hit.flag  = -1;*/

  for(int i=0; i<common::N_DSSD; i++){
    
    total_evt_mult[i][0] = 0;
    total_evt_mult[i][1] = 0;

    //std::cout << implant_hits[i][0].size() << " first side det " << i << " backside " << implant_hits[i][1].size() <<std::endl;
    
    decay_hits[i][0].clear();
    implant_hits[i][0].clear();
    decay_hits[i][1].clear();
    implant_hits[i][1].clear();
    decay_clusts[i].clear();
    implant_clusts[i].clear();

    //std::cout << implant_hits[i][0].size() << " first side det " << i << " backside " << implant_hits[i][1].size() <<std::endl;

    
    e_sum_d[i][0] = 0;
    e_sum_d[i][1] = 0;
    
    e_sum[i][0] = 0;
    e_sum[i][1] = 0;
    
    strip_max_d[i][0] = 0;
    strip_max_d[i][1] = 0;
    strip_min_d[i][0] = 128;
    strip_min_d[i][1] = 128;
  }
  
  
  /*for(int i=0; i<common::N_DSSD; ++i){
    
    total_evt_mult[i][0] = 0;
    total_evt_mult[i][1] = 0;
    
    decay_hits[i].clear();
    implant_hits[i].clear();
    decay_clusts[i].clear();
    implant_clusts[i].clear();

    e_sum_d[i][0] = 0;
    e_sum_d[i][1] = 0;

    e_sum[i][0] = 0;
    e_sum[i][1] = 0;
    
    strip_max_d[i][0] = 0;
    strip_max_d[i][1] = 0;
    strip_min_d[i][0] = 128;
    strip_min_d[i][1] = 128;

    evt_data.n_det_i[i]     = 0;  
    evt_data.n_side_i[i][0] = 0;  
    evt_data.n_side_i[i][1] = 0;  
    evt_data.x_i[i]         = -999;  
    evt_data.y_i[i]         = -999;  
    evt_data.e_i[i][0]      = -99999;  
    evt_data.e_i[i][1]      = -99999;  

    evt_data.n_det_d[i]     = 0;  
    evt_data.n_side_d[i][0] = 0;  
    evt_data.n_side_d[i][1] = 0;  
    evt_data.x_d[i]         = -999;  
    evt_data.y_d[i]         = -999;  
    evt_data.e_d[i][0]      = -99999;  
    evt_data.e_d[i][1]      = -99999;  

    evt_data.dx_d[i] = -1;
    evt_data.dy_d[i] = -1;  

    evt_data.t0     = -99999;
    evt_data.t0_ext = -99999;
    evt_data.dt     = 0; 

    }*/

  //evt_data.decay_flag   = 0;
  //evt_data.implant_flag = 0;
  
}

/* Write historgrams to file.
 */
void Analysis::WriteHistograms(){
    
  hTimeStamp->Write();
  hTimeStampExt->Write();
  hTimeStampFlag->Write();

  hEvt_TmStpDist->Write();  

  /*for(int i=0; i<4; ++i) {
    hEvt_TmStpDist[i]->Write();
    }*/

  // ******* ADC spectra *******
  for(int i=0; i<common::N_FEE64; ++i){
    if(b_mod_enabled[i]){
      hADClowCh[i]->Write();
      hADCdiscCh[i]->Write();
      hADChighCh[i]->Write();
      hCh_ADClow[i]->Write();
      hCh_ADChigh[i]->Write();
      hCh_ADCdisc[i]->Write();
            
      hElow[i]->Write();
      hEhigh[i]->Write();
      hEdisc[i]->Write();
    }
  }
  
  //******** Event monitoring ********
  for(int i=0; i<common::N_DSSD; ++i){
    
    //******* Decay events *******
    hEvt_Eside_d[i][0]->Write();
    hEvt_Eside_d[i][1]->Write();
    hEvt_ExEy_d[i][0]->Write();
    hEvt_ExEy_d[i][1]->Write();
    
    hEvt_XY_d[i]->Write();
    hEvt_X_d[i]->Write();
    hEvt_Y_d[i]->Write();

    hEvt_Mult_d[i][0]->Write();
    hEvt_Mult_d[i][1]->Write();
    hEvt_MultXY_d[i][0]->Write();
    hEvt_MultXY_d[i][1]->Write();

    //****** Implant events *******
    hEvt_Eside_i[i][0]->Write();
    hEvt_Eside_i[i][1]->Write();
    hEvt_ExEy_i[i][0]->Write();
    hEvt_ExEy_i[i][1]->Write();
    
    hEvt_XY_i[i]->Write();
    hEvt_X_i[i]->Write();
    hEvt_Y_i[i]->Write();

    hEvt_Mult_i[i][0]->Write();
    hEvt_Mult_i[i][1]->Write();
    hEvt_MultXY_i[i][0]->Write();
    hEvt_MultXY_i[i][1]->Write();
  }
  
  hEvt_Mult_impdec->Write();
  hEvt_residualE_i->Write();
  hEvt_residualE_d->Write();

  //hEvt_EPulser_d->Write();
  hEvt_pulserMult->Write();

  /*for(int i=0; i<common::N_DSSD; ++i){
    
    hEvt_X[i]->Write();
    hEvt_Y[i]->Write();
    hEvt_XY[i]->Write();

    hEvt_dX[i]->Write(); // 1:2, 1:3, 2:3
    hEvt_dY[i]->Write();
    if(i<common::N_DSSD-1) {
      hEvt_dXdX[i]->Write();
      hEvt_dYdY[i]->Write();
      }
  
    hEvt_ExEy_if[i]->Write();
    hEvt_XY_if[i]->Write();

    //hEvt_ExEy[i]->Write();
 
    hEvt_ExEy_d[i]->Write();
    hEvt_XY_d[i]->Write();   

    hEvt_ExEy_df[i]->Write();
    hEvt_XY_df[i]->Write();
    hEvt_ExEy_df2[i]->Write();
    hEvt_XY_df2[i]->Write();
    
    hTimeADClow[i]->Write();
    hTimeADCdisc[i]->Write();
    hTimeADChigh[i]->Write();
    
    hADClow_all[i][0]->Write();
    hADClow_all[i][1]->Write();
    hADClow_all[i][2]->Write();

    for(int j=0; j<2; ++j){
      
      hEvt_Eside_i[i][j]->Write();
      
      hEvt_Eside_d[i][j]->Write();
      hEvt_Eside_df[i][j]->Write();
      hEvt_Eside_df2[i][j]->Write();
      
      hEvt_Eside_if[i][j]->Write();
      
      hEvt_Multi[i][j]->Write();
      hEvt_MultidX_d[i][j]->Write();
      
      }
    }

  hEvt_Mult_d->Write();
  hEvt_MultXY_i[0]->Write();
  hEvt_MultXY_i[1]->Write();
  hEvt_MultXY_d[0]->Write();
  hEvt_MultXY_d[1]->Write();
 
  hEvt_Eaida->Write();
  
  hEvt_HitsSide->Write();
  hEvt_HitsDet->Write();
  
  hEvt_flag_d->Write();

  hEvt_MultiID->Write();
  hEvt_HitsFlag->Write();
  */

  //Write canvases
  for(int i=0; i<2; ++i) {
    cADClow[i]->Write();
    if(i==0) {
      cADChigh[i]->Write();
      cADCdisc[i]->Write();
    }
    cEall[i]->Write();
  }

  for(int i=0; i<common::N_DSSD; ++i) {
    cTimeDist[i]->Write();
  }

  cEvtE1->Write(); 
  cEvtXY->Write();
  cEvtdXdY->Write();
  cEvtMulti->Write();
  
  cEvtE_d->Write();
  cEvtXY_d->Write();
  //cEvtXY2_d->Write();
  //cEvtMulti_d->Write();

  std::cout << "\n done writing Analysis histograms to file..." << std::endl;

}

/* Reset histograms so all bin content starts back from zero.
 *
 */
void Analysis::ResetHistograms(){

  if(GetBHistograms()){
    
  hTimeStamp->Reset();
  hTimeStampExt->Reset();
  hTimeStampFlag->Reset();

  hEvt_TmStpDist->Reset();  

  // ******* ADC spectra *******
  for(int i=0; i<common::N_FEE64; ++i){
    if(b_mod_enabled[i]){
      hADClowCh[i]->Reset();
      hADCdiscCh[i]->Reset();
      hADChighCh[i]->Reset();
      hCh_ADClow[i]->Reset();
      hCh_ADChigh[i]->Reset();
      hCh_ADCdisc[i]->Reset();
            
      hElow[i]->Reset();
      hEhigh[i]->Reset();
      hEdisc[i]->Reset();
    }
  }
  
  //******** Event monitoring ********
  for(int i=0; i<common::N_DSSD; ++i){
    
    //******* Decay events *******
    hEvt_Eside_d[i][0]->Reset();
    hEvt_Eside_d[i][1]->Reset();
    hEvt_ExEy_d[i][0]->Reset();
    hEvt_ExEy_d[i][1]->Reset();
    
    hEvt_XY_d[i]->Reset();
    hEvt_X_d[i]->Reset();
    hEvt_Y_d[i]->Reset();

    hEvt_Mult_d[i][0]->Reset();
    hEvt_Mult_d[i][1]->Reset();
    hEvt_MultXY_d[i][0]->Reset();
    hEvt_MultXY_d[i][1]->Reset();

    //****** Implant events *******
    hEvt_Eside_i[i][0]->Reset();
    hEvt_Eside_i[i][1]->Reset();
    hEvt_ExEy_i[i][0]->Reset();
    hEvt_ExEy_i[i][1]->Reset();
    
    hEvt_XY_i[i]->Reset();
    hEvt_X_i[i]->Reset();
    hEvt_Y_i[i]->Reset();

    hEvt_Mult_i[i][0]->Reset();
    hEvt_Mult_i[i][1]->Reset();
    hEvt_MultXY_i[i][0]->Reset();
    hEvt_MultXY_i[i][1]->Reset();
  }
  
  hEvt_Mult_impdec->Reset();
  hEvt_residualE_i->Reset();
  hEvt_residualE_d->Reset();

  //hEvt_EPulser_d->Reset();
  hEvt_pulserMult->Reset();
  
  /*
  for(int i=0; i<common::N_FEE64; ++i){
    if(b_mod_enabled[i]){
      hADClowCh[i]->Reset();
      hADCdiscCh[i]->Reset();
      hADChighCh[i]->Reset();
      hCh_ADClow[i]->Reset();
      hCh_ADChigh[i]->Reset();
      hCh_ADCdisc[i]->Reset();
	
      hElow[i]->Reset();
      hEhigh[i]->Reset();
      hEdisc[i]->Reset();
	
    }
  }
    
  // *** Timestamp distributions ***
  hTimeADClow[0]->Reset();
  hTimeADCdisc[0]->Reset();
  hTimeADChigh[0]->Reset();
    
  hTimeStamp->Reset();
  hTimeStampExt->Reset();
  hTimeStampFlag->Reset();
        
  // ******* Event Monitoring *******
  for(int i=0; i<common::N_DSSD; ++i){

    /*hEvt_ExEy[i]->Reset();
      hEvt_X[i]->Reset();
      hEvt_Y[i]->Reset(); 
      hEvt_XY[i]->Reset();     
      hEvt_dX[i]->Reset();
      hEvt_dY[i]->Reset();      
      if(i<common::N_DSSD-1) {
      hEvt_dXdX[i]->Reset();
      hEvt_dYdY[i]->Reset();
      }

      hEvt_ExEy_if[i]->Reset();
      hEvt_XY_if[i]->Reset();

      hEvt_ExEy_d[i]->Reset();
      hEvt_ExEy_df[i]->Reset();
      hEvt_ExEy_df2[i]->Reset();
      hEvt_XY_df[i]->Reset();
      hEvt_XY_df2[i]->Reset();
      hEvt_XY_d[i]->Reset();

      hADClow_all[i][0]->Reset();
      hADClow_all[i][1]->Reset();
      hADClow_all[i][2]->Reset();
      
      for(int j=0; j<2; ++j){
      hEvt_Eside_i[i][j]->Reset();
      hEvt_Eside_if[i][j]->Reset();
      hEvt_Eside_df[i][j]->Reset();
      hEvt_Eside_df2[i][j]->Reset();
      hEvt_Multi[i][j]->Reset();
      hEvt_Eside_d[i][j]->Reset();
      hEvt_MultidX_d[i][j]->Reset();
      }
      }
    
    hEvt_TmStpDist[1]->Reset();
    hEvt_TmStpDist[0]->Reset();
    hEvt_TmStpDist[2]->Reset();
    hEvt_TmStpDist[3]->Reset();
    
    hEvt_Eaida->Reset();

    hEvt_flag_d->Reset();
    
    hEvt_HitsSide->Reset();
    hEvt_HitsDet->Reset();
    
    hEvt_MultiID->Reset();
    hEvt_HitsFlag->Reset();*/

    std::cout << "        Analysis::ResetHistograms(): all histograms have been reset..." << std::endl;
  }

}

bool Analysis::IsValidChannel(int module, int channel){
  if(module > 0 && module <= common::N_FEE64){
    if(channel >= 0 && channel < common::N_CHANNEL){
      return true;
    }
  }
  return false;
}

bool Analysis::IsChEnabled(Calibrator & my_cal_data){

  if(b_mod_enabled[my_cal_data.GetModule()-1]){
    return true;
  }
  else return false;
}

bool Analysis::SetEventTimeWindow(double value){

  if(value>0){
    event_time_window= value;
    return true;
  }
  else{
    event_time_window =0;
    return false;
  }
}

void Analysis::PrintEvent(){

  printf("\n  *EVT*   Multiplicity : ");
  printf("\n  *EVT*        N total= %i",evt_data.multiplicity);
  for(int i=0; i<common::N_DSSD; ++i) {
    printf("\n  *EVT*i       N det%i =  %i  (%i, %i)", i, evt_data.n_det_i[i], evt_data.n_side_i[i][0], evt_data.n_side_i[i][1]);
  }
  printf("\n  *EVT*i  Position     : ");
  for(int i=0; i<common::N_DSSD; ++i) {
    printf("\n  *EVT*i       X det%i =  %i   Y det%i = %i",i, i, evt_data.x_i[0], evt_data.y_i[0]);
  }
  printf("\n  *EVT*i  Energy       : ");
  for(int i=0; i<common::N_DSSD; ++i) {
    printf("\n  *EVT*i       En det%i =  %li   Ep det0= %li", i, evt_data.e_i[0][0], evt_data.e_i[0][1]);
  }

  printf("\n  *EVT*   Time         : ");
  printf("\n  *EVT*        t0=  %li      dt=  %li",evt_data.t0, evt_data.dt);
  printf("\n  *EVT*    t0_ext=  %lX",evt_data.t0_ext);

  if(evt_data.decay_flag)  printf("\n  *EVT*   Type:  decay \n");
  else printf("\n  *EVT*   Type:  implant \n");

}


void Analysis::Close(){
  if(GetBHistograms()) WriteHistograms();

  if(GetBRootTree()) out_root_tree->Write(0,TObject::kOverwrite); 

  std::cout << "!!!!!!!!! END OF ANALYSIS: " << implant_words << "/" << decay_words << " implant/decay words processed. " << closes << " EVENTS CLOSED!!!!!!!!!" << std::endl;
  std::cout << imp_down << " implants with no downstream events. " << "Of these " << imp_up << " events also have events in all upstream DSSDs." << std::endl;
  std::cout << imp_entry << " clusters with mult >= 1 found and added to clusters. " << imp_num << " front/back implant pairs found. " << std::endl;
  
  std::cout << "\n!!!!!!!!! Decays: " << dec_hits << " decay events found. " << dec_num << " front/back decay pairs add to events. " << std::endl << "\t Pulsers: " << puls_num << std::endl;

  std::cout << "Time covered by file: " << (last_ts - first_ts) * 1.0e-8 << "s, giving pulser rate of " << puls_num / ((last_ts - first_ts) * 1.0e-8) << "Hz." << std::endl;
}


void Analysis::SetBDebug(bool flag){
  b_debug= flag;
}

void Analysis::SetBHistograms(bool flag){
  b_histograms= flag;
}

/*void Analysis::SetBPushData(bool flag){
 *  b_push_data= flag;
 }*/

void Analysis::SetBRootTree(bool flag){
  b_root_tree= flag;
}

// ***********************************************************
//                       Getters...
// ***********************************************************

double Analysis::GetEventTimeWindow(){ return event_time_window; }
bool Analysis::GetBDebug(){ return b_debug; }
bool Analysis::GetBHistograms(){ return b_histograms; }
bool Analysis::GetBRootTree(){ return b_root_tree; }
int Analysis::GetMultiplicity(){ return evt_data.multiplicity; }
//bool Analysis::GetBPushData(){ return b_push_data; }

// ***********************************************************
//                         Methods
// ***********************************************************

Analysis::Analysis(){

  b_debug = false;
  b_histograms= false;
  b_root_tree= false;
  b_pulser= false;

  b_first_ts = true;

  event_count = 0;
  t_low_prev  = 0;
  t_high_prev = 0;
  t_disc_prev = 0;
  t0_aida     = 0;

  for(int i=0; i<common::N_DSSD; ++i) {
    total_evt_mult[i][0] = 0;
    total_evt_mult[i][1] = 0;
  }
  event_time_window = 2500;
  dE_i_lim= 2000;
  dX_i_lim= 15;
  E_i_min= 300;

  dE_d_lim= 3000;
  dX_d_lim= 5;
  E_d_min= 150;
  E_d_max= 3000;

  evt_data.multiplicity= 0;
  for(int i=0; i<common::N_DSSD; ++i){
    //
    e_sum_d[i][0]=0;
    e_sum_d[i][1]=0;
    //


    evt_data.n_det_i[i]=0;  
    evt_data.n_side_i[i][0]=0;  
    evt_data.n_side_i[i][1]=0;  
    evt_data.x_i[i]=-999;  
    evt_data.y_i[i]=-999;  
    evt_data.e_i[i][0]=-99999;  
    evt_data.e_i[i][1]=-99999;  

    evt_data.n_det_d[i]=0;  
    evt_data.n_side_d[i][0]=0;  
    evt_data.n_side_d[i][1]=0;  
    evt_data.x_d[i]=-999;  
    evt_data.y_d[i]=-999;  
    evt_data.e_d[i][0]=-99999;  
    evt_data.e_d[i][1]=-99999;  

    evt_data.dx_d[i]=-1;  
    evt_data.dy_d[i]=-1;  


    evt_data.t0= -99999;
    evt_data.t0_ext= -99999;
    evt_data.dt= 0; 
  }
  evt_data.decay_flag= 0;
  evt_data.implant_flag= 0;

  /*hit.t= -1;
  hit.t_ext= -1;
  hit.x= -1;
  hit.y= -1;
  hit.z= -1;
  hit.ex= -1;
  hit.ey= -1;
  hit.flag= -1;*/


}


