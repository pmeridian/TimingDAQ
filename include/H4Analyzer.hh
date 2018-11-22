#ifndef H4Analyzer_HH
#define H4Analyzer_HH
#define VME_CHANNELS 36
#define VME_TIMES 4
#define VME_SAMPLES 1024

#include "DatAnalyzer.hh"
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <math.h>
#include <assert.h>

// This is the class that should be used for parsing and analyzing
// VME data files in .dat format.

class H4Analyzer : public DatAnalyzer {
public:

#define WORDSIZE 4
  typedef unsigned long long dataTypeSize_t;
  typedef uint32_t WORD;
  
#define MAX_ADC_CHANNELS 100000
#define MAX_DIGI_SAMPLES 100000
#define MAX_TDC_CHANNELS 200
#define MAX_SCALER_WORDS 16
#define MAX_PATTERNS 16
#define MAX_PATTERNS_SHODO 16
#define SMALL_HODO_X_NFIBERS 8
#define SMALL_HODO_Y_NFIBERS 8
#define MAX_TRIG_WORDS 32
#define MAX_RO 100

  // data format used as bridge between the high level structures and the root tree
  struct treeStructData
  {
    unsigned int runNumber ;
    unsigned int spillNumber ;
    unsigned int evtNumber ;
    unsigned int evtTimeDist ;
    unsigned int evtTimeStart ;

    unsigned int 	nEvtTimes ;
    ULong64_t 	evtTime [MAX_RO] ;
    unsigned int 	evtTimeBoard [MAX_RO] ;

    //  unsigned int triggerBits ;

    unsigned int nAdcChannels ;
    //PG FIXME tranform these into vectors
    unsigned int adcBoard[MAX_ADC_CHANNELS] ;
    unsigned int adcChannel[MAX_ADC_CHANNELS] ;
    unsigned int adcData[MAX_ADC_CHANNELS] ;

    unsigned int nDigiSamples ;
    unsigned int digiFrequency[MAX_DIGI_SAMPLES] ;
    unsigned int digiStartIndexCell[MAX_DIGI_SAMPLES] ;
    unsigned int digiGroup[MAX_DIGI_SAMPLES] ;
    unsigned int digiChannel[MAX_DIGI_SAMPLES] ;
    unsigned int digiSampleIndex[MAX_DIGI_SAMPLES] ;
    unsigned int digiBoard[MAX_DIGI_SAMPLES] ;
    float digiSampleValue[MAX_DIGI_SAMPLES] ;
    float digiSampleTime[MAX_DIGI_SAMPLES] ;

    unsigned int nTdcChannels ;
    unsigned int tdcBoard[MAX_TDC_CHANNELS] ;
    unsigned int tdcChannel[MAX_TDC_CHANNELS] ;
    unsigned int tdcData[MAX_TDC_CHANNELS] ;

    unsigned int nScalerWords ;
    WORD scalerWord[MAX_SCALER_WORDS] ;
    unsigned int scalerBoard[MAX_SCALER_WORDS] ;

    unsigned int nPatterns ;
    WORD pattern[MAX_PATTERNS] ;
    WORD patternBoard[MAX_PATTERNS] ;
    WORD patternChannel[MAX_PATTERNS] ;

    unsigned int nTriggerWords;
    unsigned int triggerWords[MAX_TRIG_WORDS] ;
    unsigned int triggerWordsBoard[MAX_TRIG_WORDS] ;

    float         xIntercept;
    float         yIntercept;
    float         xSlope;
    float         ySlope;
    float         chi2;
    int           ntracks;
  } ;

  struct adcData
  {
    unsigned int board ;
    unsigned int channel ;
    unsigned int adcReadout ;
  } ;

  struct tdcData
  {
    unsigned int board ;
    unsigned int channel ;
    unsigned int tdcReadout ;
  } ;

  struct digiData
  {
    unsigned int board ;
    unsigned int group ;
    unsigned int frequency ;
    unsigned int startIndexCell ;
    unsigned int channel ;
    unsigned int sampleIndex ;
    float sampleValue;
    float sampleTime;
  } ;

  struct patternData
  {
    unsigned int board ;
    unsigned int channel ;
    unsigned int patternValue;
  };

  struct triggerWordData
  {
    unsigned int board;
    unsigned int triggerWord;
  };

  struct scalerData
  {
    unsigned int board ;
    unsigned int scalerValue;
  };

  struct timeData
  {
    unsigned int board;
    uint64_t time;
  };

  struct eventId
  {
    unsigned int runNumber;
    unsigned int spillNumber;
    unsigned int evtNumber;
  };

  struct fnal_track_reco
  {
    float         xIntercept;
    float         yIntercept;
    float         xSlope;
    float         ySlope;
    float         chi2;
    int           ntracks;
  };

  struct H4Event
  {
    H4Event (TFile * outFile, TTree * outTree) :
      outFile_ (outFile), 
      outTree_ (outTree) 
    { 
      createOutBranches (outTree_, thisTreeEvent_) ;   
    }

    ~H4Event () { }

    eventId id;
    std::vector<triggerWordData> 		triggerWords ;
    std::vector<bool> 	triggerBits ;
    std::vector<adcData> 	adcValues ; 
    std::vector<tdcData> 	tdcValues ; 
    std::vector<digiData> digiValues ; 
    std::vector<scalerData> 	scalerWords ; 
    std::vector<patternData> 	patterns ; 

    unsigned int 		evtTimeDist ;
    unsigned int 		evtTimeStart ;

    vector<timeData> 	evtTimes ;

    fnal_track_reco       track;

    void clear () ;
    void Fill () ;

  private :
  
    TFile * outFile_ ;
    TTree * outTree_ ;
    treeStructData thisTreeEvent_ ;

    // move events from the structures to the variables to be put into the root tree
    void fillTreeData (treeStructData & treeData) ;
    void createOutBranches (TTree* tree,treeStructData& treeData) ;

  } ;

  struct FTBFPixelEvent {
    double xSlope;
    double ySlope;
    double xIntercept;
    double yIntercept;
    double chi2;
    int trigger;
    int runNumber;
    Long64_t timestamp;
  };

  H4Analyzer() : DatAnalyzer(VME_CHANNELS, VME_TIMES, VME_SAMPLES, 4096, 1.) {}

  void GetCommandLineArgs(int argc, char **argv);

  void LoadCalibration();

  void InitLoop();

  int FixCorruption(int);

  int GetChannelsMeasurement();

  unsigned int GetTimeIndex(unsigned int n_ch) { return n_ch/9; }

  void Analyze();
protected:
  // Set by command line arguments or default
  TString pixel_input_file_path;
  TString calibration_file_path = "";

  // Calibration vars
  double off_mean[4][9][1024];
  double tcal[VME_TIMES][1024];

  //Raw values
  unsigned short tc[VME_TIMES]; // trigger counter bin
  unsigned short raw[VME_CHANNELS][VME_SAMPLES]; // ADC counts

  //VME binary
  unsigned int triggerNumber = -1;
  unsigned short N_corr = 0;
  unsigned long Max_corruption = 0;
  unsigned int event_time_tag = 0;
  unsigned int group_time_tag = 0;

  unsigned int ref_event_size = 0;
  unsigned int N_false = 0;

  vector<int> manual_skip = {0};

  // Pixel events variables
  FTBFPixelEvent* pixel_event;
  TFile *pixel_file = nullptr;
  TTree *pixel_tree = nullptr;

  unsigned long int idx_px_tree = 0;
  unsigned long int entries_px_tree = 0;

  float xIntercept;
  float yIntercept;
  float xSlope;
  float ySlope;
  vector<float> x_DUT;
  vector<float> y_DUT;
  float chi2;
  int ntracks;

  H4Event * event_ ;

};

#endif
