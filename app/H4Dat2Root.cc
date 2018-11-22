
// C++ includes
#include <string>
#include <assert.h>

// ROOT includes
#include <TROOT.h>

//LOCAL INCLUDES
#include "H4Analyzer.hh"

using namespace std;

int main(int argc, char **argv) {
  gROOT->SetBatch();

  H4Analyzer* analyzer = new H4Analyzer();
  analyzer->GetCommandLineArgs(argc, argv);
  analyzer->LoadCalibration();
  analyzer->RunEventsLoop();

  return 0;
}
