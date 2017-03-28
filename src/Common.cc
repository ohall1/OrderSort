#include "Common.h"
#include "DutClass.h"
#include <ctime>
#include <sys/time.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TObject.h>

#include <sstream>

void common::SaveObject(TObject* object, 
			const std::string& name, 
			DutClass* dut){
  
  std::stringstream object_name;
  object_name << dut->GetDutName()
	      << name;
    
  TFile* file_output = dut->GetFileOutput();
  file_output->cd();
  object->Write(object_name.str().c_str());
  
}

void common::SaveCanvas(TCanvas* canvas, const std::string& name, DutClass* dut){
  
  std::stringstream canvas_name;
  canvas_name << dut->GetDutName()
	      << name;
  canvas->SaveAs(canvas_name.str().c_str());
  
}

double common::getRealTime() {
  struct timeval time;
  gettimeofday(&time, NULL);
  return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}

double common::getCPUtime() {
  return (double)clock() / CLOCKS_PER_SEC;
}
