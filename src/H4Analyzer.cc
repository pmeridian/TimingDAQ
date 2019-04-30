#include "H4Analyzer.hh"
#include <assert.h>

#define DEBUG_UNPACKER 0

using namespace std;

static_assert(sizeof(unsigned int) == sizeof(float));

unsigned int castFloatToUInt(float f) {
  union { float f; unsigned int i; } u;
  u.f = f;
  return u.i;
}

float castUIntToFloat(unsigned int i) {
  union { float f; unsigned int i; } u;
  u.i = i;
  return u.f;
}

void H4Analyzer::H4Event::createOutBranches (TTree* tree,treeStructData& treeData)
{
  //Instantiate the tree branches
  tree->Branch("runNumber"	,&treeData.runNumber,"runNumber/i");
  tree->Branch("spillNumber"	,&treeData.spillNumber,"spillNumber/i");
  tree->Branch("evtNumber"	,&treeData.evtNumber,"evtNumber/i");
  tree->Branch("evtTimeDist"	,&treeData.evtTimeDist,	"evtTimeDist/i");
  tree->Branch("evtTimeStart"	,&treeData.evtTimeStart,"evtTimeStart/i");

  tree->Branch("nEvtTimes"	,&treeData.nEvtTimes,	"nEvtTimes/i"); // l-> ULong64_t
  tree->Branch("evtTime"	,treeData.evtTime,	"evtTime[nEvtTimes]/l"); // l-> ULong64_t
  tree->Branch("evtTimeBoard"	,treeData.evtTimeBoard,	"evtTimeBoard[nEvtTimes]/i"); // l-> ULong64_t

//  tree->Branch("triggerBits",&treeData.triggerBits,"triggerBits/i");

  tree->Branch("nAdcChannels"	,&treeData.nAdcChannels,"nAdcChannels/i");
  tree->Branch("adcBoard"	,treeData.adcBoard,"adcBoard[nAdcChannels]/i");
  tree->Branch("adcChannel"	,treeData.adcChannel,"adcChannel[nAdcChannels]/i");
  tree->Branch("adcData"	,treeData.adcData,"adcData[nAdcChannels]/i");

  tree->Branch("nTdcChannels"	,&treeData.nTdcChannels,"nTdcChannels/i");
  tree->Branch("tdcBoard"	,treeData.tdcBoard,"tdcBoard[nTdcChannels]/i");
  tree->Branch("tdcChannel"	,treeData.tdcChannel,"tdcChannel[nTdcChannels]/i");
  tree->Branch("tdcData"	,treeData.tdcData,"tdcData[nTdcChannels]/i");

  tree->Branch("nDigiSamples"	,&treeData.nDigiSamples,"nDigiSamples/i");
  tree->Branch("digiFrequency"	,treeData.digiFrequency,"digiFrequency[nDigiSamples]/i");
  tree->Branch("digiStartIndexCell"	,treeData.digiStartIndexCell,"digiStartIndexCell[nDigiSamples]/i");
  tree->Branch("digiGroup"	,treeData.digiGroup,"digiGroup[nDigiSamples]/i");
  tree->Branch("digiChannel"	,treeData.digiChannel,"digiChannel[nDigiSamples]/i");
  tree->Branch("digiSampleIndex",treeData.digiSampleIndex,"digiSampleIndex[nDigiSamples]/i");
  tree->Branch("digiSampleValue",treeData.digiSampleValue,"digiSample[nDigiSamples]/F");
  tree->Branch("digiSampleTime",treeData.digiSampleTime,"digiSampleTime[nDigiSamples]/F");

  tree->Branch("digiBoard"	,treeData.digiBoard,"digiBoard[nDigiSamples]/i");

  tree->Branch("nScalerWords"	,&treeData.nScalerWords,"nScalerWords/i");
  tree->Branch("scalerWord"	,treeData.scalerWord,	"scalerWord[nScalerWords]/i");
  tree->Branch("scalerBoard"	,treeData.scalerBoard,	"scalerBoard[nScalerWords]/i"); // size is Words

  tree->Branch("nPatterns"	,&treeData.nPatterns,"nPatterns/i");
  tree->Branch("pattern"	,treeData.pattern,"pattern[nPatterns]/i");
  tree->Branch("patternBoard"	,treeData.patternBoard,"patternBoard[nPatterns]/i");
  tree->Branch("patternChannel"	,treeData.patternChannel,"patternChannel[nPatterns]/i");

  tree->Branch("nTriggerWords"	,&treeData.nTriggerWords,"nTriggerWords/i");
  tree->Branch("triggerWords"	,treeData.triggerWords,"triggerWords[nTriggerWords]/i");
  tree->Branch("triggerWordsBoard",treeData.triggerWordsBoard,"triggerWordsBoard[nTriggerWords]/i");

  return ;
} 


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void H4Analyzer::H4Event::fillTreeData (treeStructData & treeData)
{
  treeData.runNumber = id.runNumber ;
  treeData.spillNumber = id.spillNumber ;
  treeData.evtNumber = id.evtNumber ;
  //  if (DEBUG_UNPACKER) printf (" =  =  =  =  =  =  FILLING EVENT %d  =  =  =  =  = \n", treeData.evtNumber) ;
  treeData.nEvtTimes = evtTimes.size () ;
  if (DEBUG_UNPACKER)
    {
      cout << "[H4Event][fillTreeData]      | FILLING " << treeData.nEvtTimes << " Times Values\n" ;
    }
  for (unsigned int i = 0 ;i<fmin (MAX_RO, evtTimes.size ()) ;++i)
    {
      treeData.evtTimeBoard[i] =  evtTimes[i].board;
      treeData.evtTime[i] =  evtTimes[i].time;
    }


  treeData.evtTimeStart = evtTimeStart ;
  treeData.evtTimeDist = evtTimeDist ;

  //treeData.triggerWord = triggerWord ;
  //cout << "[H4Event][fillTreeData]      | Trigger word: " << treeData.triggerWord << "\n" ;

//  treeData.triggerBits = 0 ;
//  for (unsigned int i = 0 ; i<fmin (32, triggerBits.size ()) ;++i)
//    treeData.triggerBits += triggerBits[i]>>i ;

  treeData.nAdcChannels = adcValues.size () ;
  if (DEBUG_UNPACKER)
    {
      cout << "[H4Event][fillTreeData]      | FILLING " << treeData.nAdcChannels << " ADC values\n" ;
    }
  for (unsigned int i = 0 ;i<fmin (MAX_ADC_CHANNELS, adcValues.size ()) ;++i)
    {
      treeData.adcBoard[i] = adcValues[i].board ;
      treeData.adcChannel[i] = adcValues[i].channel ;
      treeData.adcData[i] = adcValues[i].adcReadout ;
    }

  treeData.nTdcChannels = tdcValues.size () ;
  if (DEBUG_UNPACKER)
    {
      cout << "[H4Event][fillTreeData]      | FILLING " << treeData.nTdcChannels << " TDC values\n" ;
    }
  for (unsigned int i = 0 ;i<fmin (MAX_TDC_CHANNELS, tdcValues.size ()) ;++i)
    {
      treeData.tdcBoard[i] = tdcValues[i].board ;
      treeData.tdcChannel[i] = tdcValues[i].channel ;
      treeData.tdcData[i] = tdcValues[i].tdcReadout ;
    }

  treeData.nDigiSamples = digiValues.size () ;
  if (DEBUG_UNPACKER)
    {
      cout << "[H4Event][fillTreeData]      | FILLING " << treeData.nDigiSamples << " DIGI values\n" ;
    }

  for (unsigned int i = 0 ; i < fmin (MAX_DIGI_SAMPLES, digiValues.size ()) ;++i)
    {
      treeData.digiFrequency[i] = digiValues[i].frequency ;
      treeData.digiStartIndexCell[i] = digiValues[i].startIndexCell ;
      treeData.digiGroup[i] = digiValues[i].group ;
      treeData.digiChannel[i] = digiValues[i].channel ;
      treeData.digiBoard[i] = digiValues[i].board ;
      treeData.digiSampleIndex[i] = digiValues[i].sampleIndex ;
      treeData.digiSampleValue[i] = digiValues[i].sampleValue ;
      treeData.digiSampleTime[i] = digiValues[i].sampleTime;
    }

  treeData.nScalerWords = scalerWords.size () ;
  if (DEBUG_UNPACKER)
    {
      cout << "[H4Event][fillTreeData]      | FILLING " << treeData.nScalerWords << " Scaler words\n" ;
    }

  for (unsigned int i = 0 ; i < fmin (MAX_SCALER_WORDS, scalerWords.size()) ;++i)
    {
      treeData.scalerWord[i] = scalerWords[i].scalerValue ;
      treeData.scalerBoard[i] = scalerWords[i].board ;
    }

  treeData.nPatterns = patterns.size () ;
  if (DEBUG_UNPACKER)
    {
      cout << "[H4Event][fillTreeData]      | FILLING " << treeData.nPatterns << " Patterns\n" ;
    }

  for (unsigned int i = 0 ; i < fmin (MAX_PATTERNS, patterns.size()) ;++i)
    {
      treeData.pattern[i] 	 = patterns[i].patternValue ;
      treeData.patternBoard[i] 	 = patterns[i].board ;
      treeData.patternChannel[i] = patterns[i].channel ;
    }

  treeData.nTriggerWords = triggerWords.size () ;
  if (DEBUG_UNPACKER)
    {
      cout << "[H4Event][fillTreeData]      | FILLING " << treeData.nTriggerWords << " Trigger words\n" ;
      assert(scalerWords.size()<= MAX_TRIG_WORDS);
    }

  for (unsigned int i = 0 ; i < fmin (MAX_TRIG_WORDS, triggerWords.size()) ;++i)
    {
      treeData.triggerWords[i] 		= triggerWords[i].triggerWord ;
      treeData.triggerWordsBoard[i] 	= triggerWords[i].board ;
    }

  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void H4Analyzer::H4Event::Fill ()
{
  fillTreeData (thisTreeEvent_) ;
  outTree_->Fill () ;
  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void H4Analyzer::H4Event::clear ()
{
  id.runNumber = 0 ; // this number is unsigned
  id.spillNumber = 0 ; // this number is unsigned
  id.evtNumber = 0 ; // this number is unsigned
  triggerWords.clear() ;
//  triggerBits.clear () ;
  adcValues.clear () ; 
  tdcValues.clear () ; 
  digiValues.clear () ; 
  scalerWords.clear () ; 
  patterns.clear () ; 
  evtTimeDist = 0 ; // this number is unsigned
  evtTimeStart = 0 ; // this number is unsigned
  evtTimes.clear() ; // these numbers are unsigned
}



void H4Analyzer::GetCommandLineArgs(int argc, char **argv){
  DatAnalyzer::GetCommandLineArgs(argc, argv);

  pixel_input_file_path = ParseCommandLine( argc, argv, "pixel_input_file" );
  if (pixel_input_file_path == ""){
    if (verbose) { cout << "Pixel input file not provided" << endl; }
  }
  else {
    if (verbose) { cout << "Pixel input file: " << pixel_input_file_path.Data() << endl; }
    pixel_file = new TFile( pixel_input_file_path.Data(),"READ");
    if (!pixel_file) {std::cout << "[ERROR]: Pixel file not found" << std::endl; exit(0);}
    TString tree_name = pixel_file->GetListOfKeys()->At(0)->GetName(); //Only works if it the tree is the first key
    pixel_tree = (TTree*)pixel_file->Get(tree_name);
    if (!pixel_tree) {cout << "[ERROR]: Pixel Tree not found\n"; exit(0);}
    entries_px_tree = pixel_tree->GetEntries();
  }

  calibration_file_path = ParseCommandLine( argc, argv, "calibration_file" );
  if(calibration_file_path == ""){
    calibration_file_path = "calibration/v1740";
  }
  if (verbose) { cout << "Calibration file: " << calibration_file_path.Data() << "_bd1_group_[0-3]_[offset-dV].txt" << endl; }

  TString aux = ParseCommandLine( argc, argv, "Max_corruption" );
  if(aux != "") {
    Max_corruption = aux.Atoi();
  }
  if (verbose) { cout << "[INFO] Max corruption tollerated: " << Max_corruption << endl; }

  for(unsigned int i = 0; i<=Max_corruption; i++) manual_skip.push_back(-1);
  for(unsigned int i = 1; i<=Max_corruption; i++) {
    aux = ParseCommandLine( argc, argv, Form("NSkip%d", i) );
    if(aux == "") {
      manual_skip[i] = -1;
    }
    else {
      manual_skip[i] = aux.Atoi();
    }
  }
}

void H4Analyzer::LoadCalibration(){
  if (verbose) { cout << "---------- Loading calibrations -------------" << endl; }
  for(unsigned int ig = 0; ig < 4; ig++) {
    for(unsigned int is = 0; is < 4; is++) {
      tcal[ig][is] = 0;
      for(unsigned int ic = 0; ic < 9; ic++) off_mean[ig][ic][is] = 0;
    }
  }


  if(calibration_file_path == "ZEROS") return; //Uses defaul values (all 0)

  for( int i = 0; i < 4; i++ ){
      TString f_offset = calibration_file_path + Form("_bd1_group_%d_offset.txt",i);
      TString f_dV = f_offset;
      f_dV.ReplaceAll("_offset.txt", "_dV.txt");

      FILE* fp1 = fopen( f_offset.Data(), "r" );
      for( int k = 0; k < 1024; k++ ) {
          for( int j = 0; j < 9; j++ ){
              fscanf( fp1, "%lf ", &off_mean[i][j][k] );
          }
      }
      fclose(fp1);

      double tcal_dV[1024];
      double dV_sum = 0;
      fp1 = fopen( f_dV.Data(), "r" );
      for( int k = 0; k < 1024; k++) {
        double dummy, aux;
        fscanf( fp1, "%lf %lf %lf %lf %lf ", &dummy, &dummy, &dummy, &dummy, &aux );
        tcal_dV[k] = aux;
        dV_sum += aux;
      }
      fclose(fp1);

      for( int j = 0; j < 1024; j++) {
          tcal[i][j] = 200.0 * tcal_dV[j] / dV_sum;
      }
  }
}

void H4Analyzer::InitLoop()
{
    std::cout << "Initializing input file reader and output tree" << std::endl;
    file = new TFile(output_file_path.Data(), "RECREATE");

    if (!file->IsOpen()){
      cerr << "[ERROR]: Cannot create output file: " << output_file_path.Data() << endl;
      exit(0);
    }
    TTree * outTree = new TTree ("H4tree", "H4 testbeam tree") ;
    event_ = new H4Event(file,outTree);
    
    // Initialize the input file stream
    if ( input_file_path.EndsWith(".dat") )
      bin_file = fopen( input_file_path.Data(), "r" );
    
    if(pixel_input_file_path != "")
      {
	pixel_event = new FTBFPixelEvent;
	pixel_tree->SetBranchAddress("event", pixel_event);
      }
    pixel_tree->GetEntry(0);
    while(pixel_event->trigger <  start_evt && idx_px_tree < entries_px_tree-1) {
      idx_px_tree++;
      pixel_tree->GetEntry( idx_px_tree );
    }
}

int H4Analyzer::FixCorruption(int corruption_grp) {
  N_corr++;
  bool foundEventHeader = false;
  unsigned int long N_byte = 0;
  while (!foundEventHeader) {
    if (feof(bin_file)) return -1;

    //Look for the new header in steps of 8 bits
    N_byte++;
    char tmp = 0;
    fread( &tmp, sizeof(char), 1, bin_file);
    //find the magic word first
    if (tmp == (ref_event_size & 0xff)) {
      if (feof(bin_file)) return -1;
      fread( &tmp, sizeof(char), 1, bin_file);
      if (tmp == ((ref_event_size >> 8) & 0xff)) {
        if (feof(bin_file)) return -1;
        unsigned short tmp2 = 0;
        fread( &tmp2, sizeof(unsigned short), 1, bin_file);
        if ( tmp2 == (0xA000 + ((ref_event_size >> 16) & 0xfff)) ) {
          //now i'm good.
          foundEventHeader=true;
          cout << Form("Found a new event header after %ld bytes", N_byte) << endl;

      	  //**********************************************************
      	  //we need to increment the trigger counter an extra time
      	  //because we're skipping ahead to the next event
          if(manual_skip[N_corr] == -1) {
            if (corruption_grp<0) {
              i_evt++;
              cout << "Since corruption occured at the end of file, INCREMENTING THE i_evt" << endl;
            }
          }
          else {
            i_evt += manual_skip[N_corr];
            cout << "Manual skip set: " << manual_skip[N_corr] << " event skipped"<< endl;
          }


          if(pixel_input_file_path != ""){
            cout << "Resetting the pixel tree" << endl;
            while (idx_px_tree < entries_px_tree && i_evt >= pixel_event->trigger) {
              pixel_tree->GetEntry(idx_px_tree);
              idx_px_tree++;
            }
          }
      	  //**********************************************************

          // Reverse the two bytes read
          fseek(bin_file, -1*sizeof(unsigned int), SEEK_CUR);
          return 1;
        }
        else {
          fseek(bin_file, -1*sizeof(unsigned short), SEEK_CUR);
        }
      }
      else {
        fseek(bin_file, -1*sizeof(char), SEEK_CUR);
      }
    }
  }
  // Should not reach here. In case it does, stop the reading of the file.
  return -1;
}

// Fill tc, raw, time and amplitude
int H4Analyzer::GetChannelsMeasurement() {
    bool is_corrupted = false;
    int corruption_grp = -1;
    ResetAnalysisVariables();
    // Initialize the output variables
    for(int j = 0; j < NUM_CHANNELS; j++) {
      if(j<NUM_TIMES){ tc[j] = 0; }
      for ( int i = 0; i < NUM_SAMPLES; i++ ) {
        raw[j][i] = 0;
      }
    }

    unsigned int event_header;
    bool half_event = false;

    // first header word
    fread( &event_header, sizeof(unsigned int), 1, bin_file);
    unsigned int magicWord1 = (event_header >> 28) & 0xf;
    unsigned int eventSize =  (event_header) & 0xfffffff;
    // cout << "Event size " << eventSize << endl;
    // cout << "Ref event size " << ref_event_size << endl;

    bool ref_check = (ref_event_size == 13836);
    bool size_check = (eventSize == 6920);
    bool check =  (ref_check && size_check);

    if(i_evt == 0) {
      ref_event_size = event_header & 0xfffffff;
      if ( verbose ) { cout << "[INFO] Setting the event size to " << ref_event_size << endl; }
    }
    else if ( check ) {
      N_false++;
      half_event = true;
      // cout << "N_false " << N_false << endl;
    }
    else if ( ref_event_size != eventSize ) {
      cout << "Unexpected non-matching event size" << endl;
      cout << "Ref event size " << ref_event_size << endl;
      cout << "Event size " << eventSize << endl;
      return -1;
    }

    // second header word
    fread( &event_header, sizeof(unsigned int), 1, bin_file);
    unsigned int group_mask = event_header & 0x0f; // 4-bit channel group mask

    // third and fourth header words
    fread( &event_header, sizeof(unsigned int), 1, bin_file);
    triggerNumber = (event_header&0x3fffff);
    // if(N_false%2==1 || half_event) {cout << "triggerNumber" << triggerNumber << endl;}
    triggerNumber -= N_false/2;
    // if(N_false%2==1 || half_event) {cout << "triggerNumber" << triggerNumber << endl;}

    // cout << "triggernumber = " << triggerNumber << "\n";
    if( triggerNumber != i_evt && (N_false%2) == 0 && !half_event) {
      N_corr++;
      cout << "\tDetected missing event: N_trg " << triggerNumber << " != i_evt " << i_evt << endl;
      i_evt = triggerNumber;
      cout << "\tResetting pixel tree" << endl;
      while (idx_px_tree < entries_px_tree && i_evt > pixel_event->trigger) {
        pixel_tree->GetEntry(idx_px_tree);
        idx_px_tree++;
      }
    }

    // if(N_false%2==1 || half_event) {cout << "i_evt" << i_evt << endl;}


    fread( &event_header, sizeof(unsigned int), 1, bin_file);
    if(i_evt == start_evt) event_time_tag = event_header;
    //cout << "Evt Time: " << event_header-event_time_tag << endl;

    // check again for end of file
    if (feof(bin_file)) return -1;

    vector<unsigned int> active_groups;
    for(unsigned int k = 0; k<4; k++){
      if(group_mask & (0x1 << k)) active_groups.push_back(k);
    }

    //************************************
    // Loop over channel groups
    //************************************
    for(auto k : active_groups)
    {
      fread( &event_header, sizeof(unsigned int), 1, bin_file);
      tc[k] = (event_header >> 20) & 0xfff; // trigger counter bin

      // Check if all channels were active (if 8 channels active return 3072)
      int nsample = (event_header & 0xfff) / 3;
      if(nsample != 1024) {
        is_corrupted = true;
        corruption_grp = k;
        cout << "[WARNING]: Corruption detected at beginning of group " << k << endl;
        cout << "[WARNING]: Corruption supposed by unexpected numeber of events. Events from header " << nsample << endl;
        break;
      }

      // Define time coordinate
      time[k][0] = 0.0;
      for( int j = 1; j < 1024; j++ ){
      	time[k][j] = float(j);
      	time[k][j] = float(tcal[k][(j-1+tc[k])%1024]) + time[k][j-1];
      }

      //************************************
      // Read sample info for group
      //************************************
      unsigned short samples[9][1024];
      unsigned int temp[3];
      for ( int j = 0; j < nsample; j++ ) {
      	fread( &temp, sizeof(unsigned int), 3, bin_file);
      	samples[0][j] =  temp[0] & 0xfff;
      	samples[1][j] = (temp[0] >> 12) & 0xfff;
      	samples[2][j] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
      	samples[3][j] = (temp[1] >>  4) & 0xfff;
      	samples[4][j] = (temp[1] >> 16) & 0xfff;
      	samples[5][j] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
      	samples[6][j] = (temp[2] >>  8) & 0xfff;
      	samples[7][j] =  temp[2] >> 20;
      }

      // Trigger channel
      for(int j = 0; j < nsample/8; j++){
      	fread( &temp, sizeof(unsigned int), 3, bin_file);
      	samples[8][j*8+0] =  temp[0] & 0xfff;
      	samples[8][j*8+1] = (temp[0] >> 12) & 0xfff;
      	samples[8][j*8+2] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
      	samples[8][j*8+3] = (temp[1] >>  4) & 0xfff;
      	samples[8][j*8+4] = (temp[1] >> 16) & 0xfff;
      	samples[8][j*8+5] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
      	samples[8][j*8+6] = (temp[2] >>  8) & 0xfff;
      	samples[8][j*8+7] =  temp[2] >> 20;
      }

      //************************************
      // Loop over channels 0-8
      //************************************
      for(int jj = 0; jj < 9; jj++) {
        int j_gl = k*9 + jj;
        if ( !config->hasChannel(j_gl) ) {
          for ( int j = 0; j < 1024; j++ ) {
            raw[j_gl][j] = 0;
            channel[j_gl][j] = 0;
          }
        }
        else{
          for ( int j = 0; j < 1024; j++ ) {
            raw[j_gl][j] = samples[jj][j];
            channel[j_gl][j] = (float)(samples[jj][j]) - off_mean[k][jj][(j+tc[k])%1024];
          }
        }
      }

      fread( &event_header, sizeof(unsigned int), 1, bin_file);
      if(i_evt == start_evt) group_time_tag = event_header;
      // cout << k << ": " << event_header-group_time_tag << endl;
    } //loop over groups

    if (half_event && (N_false%2) == 1) {
      i_evt--;
      return 1;
    }

    if (feof(bin_file)) return 0;
    // Check if the following bytes corresponds to an event header. Important for skipping the event when the corruption happens during the last group;



    fread( &event_header, sizeof(unsigned int), 1, bin_file);
    magicWord1 = (event_header >> 28) & 0xf;
    // unsigned int ref_event_size =  (event_header) & 0xfffffff;
    if (feof(bin_file)) return 0;
    fread( &event_header, sizeof(unsigned int), 1, bin_file);
    int boardID = (event_header >> 27) & 0x1f;
    int pattern = (event_header >> 8) & 0xffff;
    group_mask = event_header & 0x0f; // 4-bit channel group mask

    //reverse by 2 lines
    fseek(bin_file, -2*sizeof(unsigned int), SEEK_CUR);

    if ( magicWord1 != 10 || pattern != 0 ) {
      is_corrupted = true;
      cout << "[WARNING] Following bits not matching the expected header" << endl;
    }

    if(is_corrupted) {
      cout << "Found data Corruption at end of event " << i_evt << endl;
      if(N_corr >= Max_corruption) {
        cout << "[ERROR] Corruption number over threshold (" << Max_corruption << "). Stopping acquisition." << endl;
        return -1;
      }
      cout << "Trying to skip to next event header...\n";
      return FixCorruption(corruption_grp);
    }




    return 0;
}

/*
**************************************************
Speficic analysis of VME, including telescope data
then calls main analyzer DatAnalyzer::Analyze()
**************************************************
*/
void H4Analyzer::Analyze(){

  if(pixel_input_file_path != "")
    {
      xIntercept = -999;
      yIntercept = -999;
      xSlope = -999;
      ySlope = -999;
      chi2 = -999.;
      ntracks = 0;

      while (idx_px_tree < entries_px_tree && i_evt >= (pixel_event->trigger+0)) {
	pixel_tree->GetEntry(idx_px_tree);
	if ((pixel_event->trigger+0) == i_evt) {
	  if(ntracks == 0) {
	    xIntercept = 1e-3*pixel_event->xIntercept; //um to mm
	    yIntercept = 1e-3*pixel_event->yIntercept;
	    xSlope = pixel_event->xSlope;
	    ySlope = pixel_event->ySlope;
	    chi2 = pixel_event->chi2;
	  }
	  ntracks++;
	  idx_px_tree++;
	}
	else if (i_evt > (pixel_event->trigger+0)) {
	  cout << "[ERROR] Pixel tree not ordered" << endl;
	  exit(0);
	}
      }

    }

  event_->clear();

  event_->id.runNumber = pixel_event->runNumber;
  event_->id.spillNumber = 1;
  event_->id.evtNumber = i_evt;

  timeData timeStamp; //timestamp from pixel readout
  timeStamp.board=1;
  timeStamp.time=pixel_event->timestamp;
  event_->evtTimes.push_back(timeStamp);


  //encode tracks into an adcBoard....
  adcData trackData;
  trackData.board=1;
  trackData.channel=1; //ch1 ntracks
  trackData.adcReadout=ntracks;
  event_->adcValues.push_back(trackData);
  trackData.channel=2; //ch2 chi2
  trackData.adcReadout=castFloatToUInt(chi2);
  event_->adcValues.push_back(trackData);
  trackData.channel=3; //ch3 x
  trackData.adcReadout=castFloatToUInt(xIntercept);
  event_->adcValues.push_back(trackData);
  trackData.channel=4; //ch4 y
  trackData.adcReadout=castFloatToUInt(yIntercept);
  event_->adcValues.push_back(trackData);
  trackData.channel=5; //ch5 xSlope
  trackData.adcReadout=castFloatToUInt(xSlope);
  event_->adcValues.push_back(trackData);
  trackData.channel=6; //ch6 ySlope
  trackData.adcReadout=castFloatToUInt(ySlope);
  event_->adcValues.push_back(trackData);

  for (int ich=0;ich<36;++ich)
    for (int isample=0;isample<1024;++isample)
      {	  
	digiData aDigiSample ;
	aDigiSample.board = 1 ;
	aDigiSample.channel = ich%9 ;
	aDigiSample.group = ich/9 ;
	aDigiSample.frequency = 0 ;
	aDigiSample.startIndexCell = tc[aDigiSample.group] ;
	aDigiSample.sampleIndex = isample ;
	
	// if (!saveraw_)
	aDigiSample.sampleValue = channel[ich][isample]; //apply no conversion to mV
	// else
	// aDigiSample.sampleValue = (fnalEvent_->raw[ich][isample]);
	
	// if (!saveraw_)
	aDigiSample.sampleTime = (time[aDigiSample.group][isample])*204.8/200.; //fix for bug in Caltech ntuple producer
	// else
	// aDigiSample.sampleTime = 0.2*isample;
	
	event_->digiValues.push_back (aDigiSample) ;
      }

  event_->Fill();    
}
