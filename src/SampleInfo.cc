#include "interface/SampleInfo.h"

SampleInfo::SampleInfo()
{
  sampleName="";
  //mcEvtType = 0;
  gen_weight = 0.0;
  pileup.clear();
  truepileup.clear();
  ///puWeights.clear();
  ///truepuWeights.clear();
  hepNUP = -99;
  topPtWeights.clear();
  hepIDUP.clear();
  hepISTUP.clear();
}

SampleInfo::~SampleInfo()
{
}
