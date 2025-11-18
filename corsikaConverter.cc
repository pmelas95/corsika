
#include "TVector3.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TDatabasePDG.h"
#include "TObjString.h"

#define CM2M 0.01
#define M2CM 100

using namespace std;

enum class Format { UNDEFINED, NORMAL, COMPACT };
static constexpr unsigned _fbsize_words = 5733 + 2;


//static constexpr unsigned obsv = 1100;
static constexpr unsigned obsv = 455;

static constexpr float x_border[2] = { -1.8,  1.8 };   // half of 3.6 m
static constexpr float y_border[2] = { -3.6,  3.6 };   // half of 7.2 m
static constexpr float z_border[2] = { -3.05, 3.05 };  // half of 6.1 m
static constexpr float box_width[3] = {x_border[1]-x_border[0], y_border[1]-y_border[0], z_border[1]-z_border[0]};

//static constexpr float ceiling     = 470000;
static constexpr float ceiling     = 377000;
//static constexpr float ceiling   = 378500;

union {
	float fl[_fbsize_words];
	unsigned in[_fbsize_words];
	char ch[sizeof(fl[0])*_fbsize_words];
} buf;

const std::map<unsigned int, int> corsikaToPdgId = {
	{1, 22},     // gamma
	{2, -11},    // e+
	{3, 11},     // e-
	{5, -13},    // mu+
	{6, 13},     // mu-
	{7, 111},    // pi0
	{8, 211},    // pi+
	{9, -211},   // pi-
	{10, 130},   // K0_L
//	{11, 321},   // K+
//	{12, -321},  // K-
	{13, 2112},  // n
	{14, 2212},  // p
	{15, -2212}, // pbar
	{16, 310},   // K0_S
	{17, 221},   // eta
	{18, 3122},  // Lambda
	{19, 3222},  // Sigma+
	{20, 3212},  // Sigma0
	{21, 3112},  // Sigma-
	{22, 3322},  // Cascade0
	{23, 3312},  // Cascade-
	{24, 3334},  // Omega-
	{25, -2112}, // nbar
	{26, -3122}, // Lambdabar
	{27, -3112}, // Sigma-bar
	{28, -3212}, // Sigma0bar
	{29, -3222}, // Sigma+bar
	{30, -3322}, // Cascade0bar
	{31, -3312}, // Cascade+bar
	{32, -3334}, // Omega+bar

	{50, 223},    // omega
	{51, 113},    // rho0
	{52, 213},    // rho+
	{53, -213},   // rho-
	{54, 2224},   // Delta++
	{55, 2214},   // Delta+
	{56, 2114},   // Delta0
	{57, 1114},   // Delta-
	{58, -2224},  // Delta--bar
	{59, -2214},  // Delta-bar
	{60, -2114},  // Delta0bar
	{61, -1114},  // Delta+bar
	{62, 10311},  // K*0
	{63, 10321},  // K*+
	{64, -10321}, // K*-
	{65, -10311}, // K*0bar
	{66, 12},     // nu_e
	{67, -12},    // nu_ebar
	{68, 14},     // nu_mu
	{69, -14},    // nu_mubar

	{116, 421},  // D0
	{117, 411},  // D+
	{118, -411}, // D-bar
	{119, -421}, // D0bar
	{120, 431},  // D+_s
	{121, -431}, // D-_sbar
	{122, 441},  // eta_c
	{123, 423},  // D*0
	{124, 413},  // D*+
	{125, -413}, // D*-bar
	{126, -423}, // D*0bar
	{127, 433},  // D*+_s
	{128, -433}, // D*-_s

	{130, 443}, // J/Psi
	{131, -15}, // tau+
	{132, 15},  // tau-
	{133, 16},  // nu_tau
	{134, -16}, // nu_taubar

	{137, 4122},  // Lambda+_c
	{138, 4232},  // Cascade+_c
	{139, 4132},  // Cascade0_c
	{140, 4222},  // Sigma++_c
	{141, 4212},  // Sigma+_c
	{142, 4112},  // Sigma0_c
	{143, 4322},  // Cascade'+_c
	{144, 4312},  // Cascade'0_c
	{145, 4332},  // Omega0_c
	{149, -4122}, // Lambda-_cbar
	{150, -4232}, // Cascade-_cbar
	{151, -4132}, // Cascade0_cbar
	{152, -4222}, // Sigma--_cbar
	{153, -4212}, // Sigma-_cbar
	{154, -4112}, // Sigma0_cbar
	{155, -4322}, // Cascade'-_cbar
	{156, -4312}, // Cascade'0_cbar
	{157, -4332}, // Omega0_cbar
	{161, 4224},  // Sigma*++_c
	{162, 1214},  // Sigma*+_c
	{163, 4114},  // Sigma*0_c

	{171, -4224},     // Sigma*--_cbar
	{172, -1214},     // Sigma*-_cbar
	{173, -4114},     // Sigma*0_cbar
	{176, 511},       // B0
	{177, 521},       // B+
	{178, -521},      // B-bar
	{179, -511},      // B0bar
	{180, 531},       // B0_s
	{181, -531},      // B0_sbar
	{182, 541},       // B+_c
	{183, -541},      // B-_cbar
	{184, 5122},      // Lambda0_b
	{185, 5112},      // Sigma-_b
	{186, 5222},      // Sigma+_b
	{187, 5232},      // Cascade0_b
	{188, 5132},      // Cascade-_b
	{189, 5332},      // Omega-_b
	{190, -5112},     // Lambda0_bbar
	{191, -5222},     // Sigma+_bbar
	{192, -5112},     // Sigma-_bbar
	{193, -5232},     // Cascade0_bbar
	{194, -5132},     // Cascade+_bbar
	{195, -5332},      // Omega+_bbar

	{201, 1000010020}, // Deuteron
	{301, 1000010030}, // Tritium
	{402, 1000020040}, // alpha
	{5626, 1000260560}, // Iron
	{1206, 1000080120}, // Carbon
	{1407, 1000070140}, // Nitrogen
	{1608, 1000080160}, // Oxygen
	{2713, 1000130270}, // Aluminum
	{3216, 1000160320}, // Sulfur
	{2814, 1000140280}, // Silicon
	{9900, 22}          // Cherenkov gamma
};

// Structure to represent each muon and its mother


int main(int argc, char *argv[]) {
	if (argc == 1) {
		std::cout << "Usage: corsikaConverter INPUT" << std::endl;
		return 0;
	}

	ifstream *input = new ifstream(argv[1]);;
    int current_event_number = -1;
    int event_count = 0;
    int run_number = -1;
    Format _infmt = Format::UNDEFINED;
/*
                                std::vector<float> x;
                                std::vector<float> z;
                                std::vector<float> px;
                                std::vector<float> py;
                                std::vector<float> pz;
                                std::vector<int> pdg;
                                std::vector<int> startingx;
                                std::vector<int> startingy;
                                std::vector<int> time;
                                std::vector<int> startingz;
                                std::vector<int> Gene;
*/
	TDatabasePDG *pdg_database = TDatabasePDG::Instance();
	pdg_database->ReadPDGTable();

	char outputFilename[256];
	strcpy(outputFilename, argv[1]);
	strcat(outputFilename,".root");

	TFile *gRooTrackerFile = new TFile(outputFilename,"RECREATE");
	TTree *gRooTracker2 = new TTree("gRooTracker2", "gRooTracker2");
        TTree *gRooTracker = new TTree("gRooTracker", "gRooTracker");
	static constexpr unsigned kMaxParticles = 39;
	TRandom *xyzRandom = new TRandom();
	//TRandom *particleIndex = new TRandom();

	Int_t EvtNum;
	TBits *EvtFlags = new TBits();
	TObjString *EvtCode = new TObjString("");
	Double_t EvtXSec = 1.;
	Double_t EvtDXSec = 1.;
	Double_t EvtWght = 1.;
	Double_t EvtProb = 1.;
	Int_t StdHepPdg[kMaxParticles];
	Int_t StdHepN;
        Double_t EvtVtx [4];
	Int_t StdHepStatus[kMaxParticles];
	Int_t StdHepRescat[kMaxParticles];
	Double_t StdHepP4[kMaxParticles][4];
	Double_t StdHepX4[kMaxParticles][4];
	Double_t StdHepPolz[kMaxParticles][3];
	Int_t StdHepFd[kMaxParticles];
	Int_t StdHepLd[kMaxParticles];
	Int_t StdHepFm[kMaxParticles];
	Int_t StdHepLm[kMaxParticles];
        Int_t Has_siblings[kMaxParticles];
        Double_t  Has_costheta[kMaxParticles];
        Double_t  Has_DT_MIN[kMaxParticles];
        Double_t  Has_MIN[kMaxParticles];
        Double_t  Has_DT[kMaxParticles];
        Double_t  Has_MIN_DT[kMaxParticles];
        Double_t  Has_Y[kMaxParticles];
        Double_t  Has_Gene[kMaxParticles];
	Int_t NuParentPdg;
	Int_t NuParentDecMode;
	Double_t NuParentDecP4[4];
	Double_t NuParentDecX4[4];
	Double_t NuParentProX4[4];
	Double_t NuParentProP4[4];
	Int_t NuParentProNVtx;

        std::vector<float> DT_MIN,DT,Generation,ANGLE,Xmin,Ymin,Zmin,ENERGY,MIN,MIN_DT;
        gRooTracker->Branch("Generation", &Generation);
        gRooTracker->Branch("ANGLE", &ANGLE);
        gRooTracker->Branch("Xmin", &Xmin);
        gRooTracker->Branch("Ymin", &Ymin);
        gRooTracker->Branch("Zmin", &Zmin);
        gRooTracker->Branch("ENERGY", &ENERGY);


	gRooTracker2->Branch("EvtNum", &EvtNum, "EvtNum/I");
        gRooTracker2->Branch("EvtFlags", "TBits", &EvtFlags, 32000, 1);
        gRooTracker2->Branch("EvtCode", "TObjString", &EvtCode, 32000, 1);
	gRooTracker2->Branch("EvtXSec", &EvtXSec, "EvtXsec/D");
	gRooTracker2->Branch("EvtDXSec", &EvtDXSec, "EvtDXsec/D");
	gRooTracker2->Branch("EvtWght", &EvtWght, "EvtWght/D");
	gRooTracker2->Branch("EvtProb", &EvtProb, "EvtProb/D");
	gRooTracker2->Branch("StdHepN", &StdHepN, "StdHepN/I");
        gRooTracker2->Branch("EvtVtx", EvtVtx, "EvtVtx[4]/D");
	gRooTracker2->Branch("StdHepPdg", StdHepPdg, "StdHepPdg[StdHepN]/I");
	gRooTracker2->Branch("StdHepStatus", StdHepStatus, "StdHepStatus[StdHepN]/I");
	gRooTracker2->Branch("StdHepRescat", StdHepRescat, "StdHepRescat[StdHepN]/I");
	gRooTracker2->Branch("StdHepP4", StdHepP4, "StdHepP4[StdHepN][4]/D");
	gRooTracker2->Branch("StdHepX4", StdHepX4, "StdHepX4[StdHepN][4]/D");
	gRooTracker2->Branch("StdHepPolz", StdHepPolz, "StdHepPolz[StdHepN][3]/D");
	gRooTracker2->Branch("StdHepFd", StdHepFd, "StdHepFd[StdHepN]/I");
	gRooTracker2->Branch("StdHepLd", StdHepLd, "StdHepLd[StdHepN]/I");
	gRooTracker2->Branch("StdHepFm", StdHepFm, "StdHepFm[StdHepN]/I");
	gRooTracker2->Branch("StdHepLm", StdHepLm, "StdHepLm[StdHepN]/I");

	gRooTracker2->Branch("NuParentPdg", &NuParentPdg, "NuParentPdg/I");
	gRooTracker2->Branch("NuParentDecMode", &NuParentDecMode, "NuParentDecMode/I");
	gRooTracker2->Branch("NuParentDecP4", NuParentDecP4, "NuParentDecP4[4]/D");
	gRooTracker2->Branch("NuParentDecX4", NuParentDecX4, "NuParentDecX4[4]/D");
	gRooTracker2->Branch("NuParentProP4", NuParentProP4, "NuParentProP4[4]/D");
	gRooTracker2->Branch("NuParentProX4", NuParentProX4, "NuParentProX4[4]/D");
	gRooTracker2->Branch("NuParentProNVtx", &NuParentProNVtx, "NuParentProNVtx/I");
        //gRooTracker2->Branch("Has_siblings", &Has_siblings, "Has_siblings[StdHepN]/I");
        //gRooTracker2->Branch("Has_costheta", &Has_costheta, "Has_costheta[StdHepN]/D");
        gRooTracker2->Branch("Has_MIN_DT", &Has_MIN_DT, "Has_MIN_DT[StdHepN]/D");
        gRooTracker2->Branch("Has_MIN", &Has_MIN, "Has_MIN[StdHepN]/D");
        gRooTracker2->Branch("Has_DT", &Has_DT, "Has_DT[StdHepN]/D");
        gRooTracker2->Branch("Has_DT_MIN", &Has_DT_MIN, "Has_DT_MIN[StdHepN]/D");
        gRooTracker2->Branch("Has_Y", &Has_Y, "Has_Y[StdHepN]/D");
        gRooTracker2->Branch("Has_Gene", &Has_Gene, "Has_Gene[StdHepN]/D");

bool found_marker = false;

    while( input->read(buf.ch, 4)) {
      unsigned reclen = buf.in[0];
      // CORSIKA records are in units of 4 bytes
      if(reclen % 4) {
        throw std::runtime_error("Error: record size not a multiple of 4");
      }

      // We will be looking at at least 8 bytes to determine the
      // input file format, and all real CORSIKA records are longer
      // than that.
      if(reclen < 2*4) {
        throw std::runtime_error("Error: reclen too small");
      }

      // Read the full record
      if(!input->read(buf.ch, reclen)) {
        break;
      }

      // Determine the format and and store the decision for future blocks.
      // We are starting file read, so should see the RUNH marker
      // In COMPACT format each block is preceded by 4 bytes
      // giving the size of the block in words.

      if(!strncmp(buf.ch+0, "RUNH", 4)) {
        std::cout<<"Reading NORMAL format"<<std::endl;
        _infmt = Format::NORMAL;
      }
      else if(!strncmp(buf.ch+4, "RUNH", 4)) {
        std::cout<<"Reading COMPACT format"<<std::endl;
        _infmt = Format::COMPACT;
      }
      else {
        throw std::runtime_error("Error: did not find the RUNH record to determine COMPACT flag");
      }

      unsigned iword = 0;
      if(_infmt == Format::COMPACT) {
        // Move to the beginning of the actual block
        ++iword;
      }

      break;

    }

	input->clear();
    input->seekg(0, ios::beg);
    unsigned particleCounter = 0;
	// bool eventWithParticle = false;
	// FORTRAN sequential records are prefixed with their length
	// in a 4-byte word
	while( input->read(buf.ch, 4)) {

		unsigned reclen = buf.in[0];
		// CORSIKA records are in units of 4 bytes
		if(reclen % 4) {
			throw std::runtime_error("Error: record size not a multiple of 4");
		}

		// We will be looking at at least 8 bytes to determine the
		// input file format, and all real CORSIKA records are longer
		// than that.
		if(reclen < 2*4) {
			throw std::runtime_error("Error: reclen too small");
		}

		if(reclen > 4*_fbsize_words) {
			throw std::runtime_error("Error: reclen too big");
		}

		// Read the full record
		if(!input->read(buf.ch, reclen)) {
			break;
		}

		//================================================================

/*
Generation.clear();
DT_MIN.clear();
iDT.clear();
ANGLE.clear();
Xmin.clear();
Ymin.clear();
Zmin.clear();
ENERGY.clear();
MIN.clear();
MIN_DT.clear();
*/

		// Go over blocks in the record
		for(unsigned iword = 0; iword < reclen/4; ) {

			unsigned block_words = (_infmt == Format::COMPACT) ? buf.in[iword] : 273;

			if(!block_words) {
				throw std::runtime_error("Got block_words = 0\n");
			}

			if(_infmt == Format::COMPACT) {
				// Move to the beginning of the actual block
				++iword;
			}

			std::string event_marker = (_infmt == Format::NORMAL || !event_count) ? "EVTH" : "EVHW";

			// Determine the type of the data block
			if(!strncmp(buf.ch+4*iword, "RUNH", 4)) {
				run_number = lrint(buf.fl[1+iword]);
				float lowE = float(buf.fl[16+iword]);
				float highE = float(buf.fl[17+iword]);
				std::cout << "Start Run " << run_number << std::endl;
				std::cout << "Primary proton energy between " << lowE << " and " << highE << std::endl;
			}
			else if (!strncmp(buf.ch+4*iword, "RUNE", 4)) {
				int end_run_number = lrint(buf.fl[1+iword]);
				int end_event_count = lrint(buf.fl[2+iword]);
				std::cout << "End Run " << run_number << std::endl;

				if(end_run_number != run_number) {
					throw std::runtime_error("Error: run number mismatch in end of run record\n");
				}

				if(event_count != end_event_count) {
					std::cerr<<"RUNE: _event_count = "<<event_count<<" end record = "<<end_event_count<<std::endl;
					throw std::runtime_error("Error: event count mismatch in end of run record\n");
				}

                                //gRooTracker->Write();
				gRooTracker2->Write();
				gRooTrackerFile->Close();

				// Exit the read loop at the end of run
				return 0;
			}
			else if(!strncmp(buf.ch+4*iword, event_marker.data(), 4)) {

				current_event_number = lrint(buf.fl[1+iword]);
  
                          //   x.clear();  z.clear(); px.clear(); py.clear(); pz.clear();pdg.clear(); Gene.clear(); time.clear();

				++event_count;
			}
			else if(!strncmp(buf.ch+4*iword, "EVTE", 4)) {
				int end_event_number = lrint(buf.fl[1+iword]);

				if(end_event_number != current_event_number) {
					throw std::runtime_error("Error: event number mismatch in end of event record\n");
				}
			}
			else {
				// const float x_shift = zxRandom->Uniform(-box_width[0]/2,box_width[0]/2);
				// const float z_shift = zxRandom->Uniform(-box_width[0]/2,box_width[0]/2);
                            
                                  std::vector<float> x;
                                std::vector<float> z;
                                std::vector<float> px;
                                std::vector<float> py;
                                std::vector<float> pz;
                                std::vector<int> pdg;
                                std::vector<int> startingx;
                                std::vector<int> startingy;
                                std::vector<int> time;
                                std::vector<int> startingz;
                                std::vector<int> Gene;
                               
                    float x_prod = -999;
                    float z_prod = -999;
                    float y_prod = -999;

				for (unsigned i_part = 0; i_part < block_words; i_part+=7) {
					unsigned id = buf.fl[iword + i_part] / 1000;
unsigned offset = iword + i_part;
int full_id = static_cast<int>(buf.fl[offset]);
int generation = (abs(full_id) % 1000) / 10;


//if(generation!=52){
//if(id==5||id==6){
/////////////////////////////////////////////////////////////////////////////////////
/*
                float t = buf.fl[offset + 6];
                float x = buf.fl[offset + 4];
                float y = buf.fl[offset + 5];

                float px = buf.fl[offset + 1];
                float py = buf.fl[offset + 2];
                float pz = buf.fl[offset + 3];

if (id == 75 || id == 76) {
                    std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  <<std::endl;
std::cout << "id: " << id
          << ", full_id: " << buf.fl[iword + i_part]
          << ", t: " << t
          << ", x: " << x << ", y: " << y
          << ", px: " << px << ", py: " << py << ", pz: " << pz
          << std::endl;
std::cout << "----------------------------------" << std::endl;

                }else{

std::cout << "id: " << id
          << ", full_id: " << buf.fl[iword + i_part]
          << ", t: " << t
          << ", x: " << x << ", y: " << y
          << ", px: " << px << ", py: " << py << ", pz: " << pz
          << std::endl;
}


}
*/
///////////////////////////////////////////////////////////////////
                if (id == 75 || id == 76) {
                    
                     x_prod = buf.fl[iword + i_part + 4] * CM2M;
                     z_prod = buf.fl[iword + i_part + 5] * CM2M;
                     y_prod = buf.fl[iword + i_part + 6] * CM2M;
                     continue;
                }


                                        if (id == 0) continue;
					x.push_back(buf.fl[iword + i_part + 4] * CM2M);
					z.push_back(buf.fl[iword + i_part + 5] * CM2M);
                                        time.push_back(buf.fl[iword + i_part + 6] * CM2M);
					px.push_back(buf.fl[iword + i_part + 1]);
					pz.push_back(buf.fl[iword + i_part + 2]);
					py.push_back(-buf.fl[iword + i_part + 3]);
                                        auto it = corsikaToPdgId.find(id);
					pdg.push_back(it->second);
if (id == 5 || id == 6)  
    Gene.push_back(generation);
else
    Gene.push_back(-999);                                        
                                        startingx.push_back(x_prod);
                                        startingz.push_back(z_prod);
                                        startingy.push_back(y_prod);
                                        x_prod = -999;
                                        z_prod = -999;
                                        y_prod = -999;                                      
}



Generation.clear();
MIN.clear();
DT.clear();
ANGLE.clear();
Xmin.clear();
Ymin.clear();
Zmin.clear();
ENERGY.clear();


//////////////// # === Bundle Tagging Based on Proximity ===
int counter =-999;

// === Loop only over muons with generation != 52 ===
for (size_t ik = 0; ik < startingx.size(); ++ik) {

 counter =0;
    float closest_dt = std::numeric_limits<float>::max();
    float best_min_dist = std::numeric_limits<float>::max();
    float best_cos_theta = -2.0;

    float closest_dt2 = std::numeric_limits<float>::max();
    float best_min_dist2 = std::numeric_limits<float>::max();
    float best_cos_theta2 = -2.0;

   Ymin.push_back(startingy[ik]);



if(abs(pdg[ik]) == 13) Generation.push_back(Gene[ik]);
    if (abs(pdg[ik]) == 13){
//std::cout<< "Generation= "<< Gene[ik] << " x= "<< startingx[ik] << " y= "<< startingy[ik] << " z= "<< startingz[ik]<< " t= "<< time[ik] <<std::endl; 

//if (Gene[ik] == 52 || Gene[ik] <= 2 ) continue;
// if (Gene[ik] != 52 || Gene[ik] > 2 ) continue;
//    if (Gene[ik] == 52) continue;
counter =counter+1;


    TVector3 pos1(startingx[ik], startingy[ik], startingz[ik]);
    //TVector3 dir1(px[ik], py[ik], pz[ik]);
    TVector3 dir1(x[ik], obsv, z[ik]);
    dir1 = dir1.Unit();


//DT.push_back(time[ik]);
//Xmin.push_back(startingx[ik]);
//Ymin.push_back(startingy[ik]);
//Zmin.push_back(startingz[ik]);

/////////////////////////////////////////////////////////////
    for (size_t iz = 0; iz < startingx.size(); ++iz) {
        if (iz == ik) continue;
        if (abs(pdg[iz]) != 13) continue;
        //if (Gene[iz] == 52) continue;

        float ax = startingx[ik] - startingx[iz];
        float ay = startingy[ik] - startingy[iz];
        float az = startingz[ik] - startingz[iz];

        float minimum_distance = std::sqrt(ax*ax + ay*ay + az*az);
        float dt = abs(time[ik]-time[iz]);

///////////////////////////////////////////////////////////////////////
        TVector3 pos2(startingx[iz], startingy[iz], startingz[iz]);
        //TVector3 dir2(px[iz], py[iz], pz[iz]);
        TVector3 dir2(x[iz], obsv, z[iz]);
        dir2 = dir2.Unit();

        TVector3 cross = dir1.Cross(dir2);
        double denom = cross.Mag();
        double min_dist;

        if (denom < 1e-8) {
            min_dist = (pos2 - pos1).Cross(dir1).Mag();  // nearly parallel
        } else {
            TVector3 w0 = pos2 - pos1;
            min_dist = std::abs(w0.Dot(cross)) / denom;
        }


TVector3 dir3(x[ik] - startingx[ik], obsv - startingy[ik], z[ik] - startingz[ik]);
TVector3 dir4(x[iz] - startingx[iz], obsv - startingy[iz], z[iz] - startingz[iz]);
float cos_theta = dir3.Angle(dir4) * (180.0 / TMath::Pi());

  //    float cos_theta = dir1.Dot(dir2);
//float cos_theta = dir1.Angle(dir2) * (180.0 / TMath::Pi());
 
// if (minimum_distance < best_min_dist) {
// Keep only the closest in time
        if (dt < closest_dt) {
            closest_dt = dt;
            best_min_dist = minimum_distance;  // or use min_dist if you prefer skew line distance
            best_cos_theta = cos_theta;
        }


        if (minimum_distance < best_min_dist2) {
            closest_dt2 = dt;
            best_min_dist2 = minimum_distance;  // or use min_dist if you prefer skew line distance
            best_cos_theta2 = cos_theta;
        }


/////////////////////////////////////////////////////////////////////////////

}}
const float momentum2 = sqrt(px[ik]*px[ik] + py[ik]*py[ik] + pz[ik]*pz[ik]);
const float energy2 = sqrt(0.1057*0.1057 + momentum2*momentum2);

if (closest_dt2 == std::numeric_limits<float>::max() && best_min_dist2 != std::numeric_limits<float>::max()) {
    std::cout << "Warning: MIN_DT is max but MIN is not! This shouldn't happen." << std::endl;
}

    DT.push_back(closest_dt);
    DT_MIN.push_back(best_min_dist);
    ANGLE.push_back(best_cos_theta);
    ENERGY.push_back(energy2);
MIN.push_back(best_min_dist2);
MIN_DT.push_back(closest_dt2);

}

//gRooTracker->Fill();
//std::cout<<"DTsize= "<< Ymin.size() << "  x.size()= "<< x.size()<<std::endl;
//std::cout<<"startingx.size()= "<< Gene.size() << "  x.size()= "<< x.size()<<std::endl;

				for (unsigned i_part = 0; i_part < x.size(); i_part++) {
					const float momentum = sqrt(px[i_part]*px[i_part] + py[i_part]*py[i_part] + pz[i_part]*pz[i_part]);
					const std::vector<float> dir_vector{ px[i_part]/momentum, py[i_part]/momentum, pz[i_part]/momentum };
					const float random_y = xyzRandom->Uniform(y_border[0],y_border[1]);
					const float t = (random_y-ceiling)/dir_vector[1];
					const float x_particle = x[i_part] + dir_vector[0]*t;
					const float z_particle = z[i_part] + dir_vector[2]*t;
					const float deltax[2] = {x_particle-x_border[1],x_particle-x_border[0]};
					const float deltaz[2] = {z_particle-z_border[1],z_particle-z_border[0]};
					const float x_shift = xyzRandom->Uniform(*std::min_element(deltax,deltax+2), *std::max_element(deltax,deltax+2));
					const float z_shift = xyzRandom->Uniform(*std::min_element(deltaz,deltaz+2), *std::max_element(deltaz,deltaz+2));
                                        //if(i_part== x.size()-1) std::cout<<"t= "<< t <<std::endl;
                                       //                     for (int nx = -5; nx <= 5; nx++) {
                                       //                   for (int nz = -5; nz <= 5; nz++) {

							memset(StdHepPdg, 0, sizeof(StdHepPdg));
							memset(StdHepP4, 0, sizeof(StdHepP4));
							memset(StdHepX4, 0, sizeof(StdHepX4));
							memset(StdHepStatus, 0, sizeof(StdHepStatus));
							memset(EvtVtx, 0, sizeof(EvtVtx));
                                                        memset(Has_siblings, 0, sizeof(Has_siblings));
                                                        memset(Has_costheta, 0, sizeof(Has_costheta));


                                                        memset(Has_MIN, 0, sizeof(Has_MIN));
                                                        memset(Has_MIN_DT, 0, sizeof(Has_MIN_DT));
                                                        memset(Has_DT_MIN, 0, sizeof(Has_DT_MIN));
                                                        memset(Has_DT, 0, sizeof(Has_DT));
                                                        memset(Has_Gene, 0, sizeof(Has_Gene));
                                                        memset(Has_Y, 0, sizeof(Has_Y));

							StdHepN = 0;
							StdHepPdg[StdHepN] = pdg[i_part];
							StdHepStatus[StdHepN] = 1;
							StdHepP4[StdHepN][0] = px[i_part];
							StdHepP4[StdHepN][1] = py[i_part]; // Here we switch z and y because in the GDML
							StdHepP4[StdHepN][2] = pz[i_part]; // y is the vertical coordinate
							const float mass = pdg_database->GetParticle(pdg[i_part])->Mass();
							const float energy = sqrt(mass*mass + momentum*momentum);
							StdHepP4[StdHepN][3] = energy;
							StdHepX4[StdHepN][0] = 0;
							StdHepX4[StdHepN][1] = 0;
							StdHepX4[StdHepN][2] = 0;
							StdHepX4[StdHepN][3] = 0;
                                                        Has_Y[StdHepN] = Ymin[i_part];
                                                        Has_Gene[StdHepN] = Gene[i_part];
                                                        Has_MIN[StdHepN] = MIN[i_part];
                                                        Has_MIN_DT[StdHepN] = MIN_DT[i_part];
                                                        Has_DT[StdHepN] = DT[i_part];
                                                        Has_DT_MIN[StdHepN] = DT_MIN[i_part];

                                                        //EvtVtx[0] = x[i_part] - x_shift + nx*(x_border[1]-x_border[0]);
							//EvtVtx[1] = ceiling;
							//EvtVtx[2] = z[i_part] - z_shift + nz*(z_border[1]-z_border[0]);

                                                        EvtVtx[0] = x[i_part]* M2CM;
                                                        EvtVtx[1] = t;
                                                        EvtVtx[2] = z[i_part]* M2CM;
							EvtNum = particleCounter;
							StdHepN++;
							gRooTracker2->Fill();
							particleCounter++;
                                                        //}}

			

}

			}


			// Move to the next block
			iword += block_words;

		} // loop over blocks in a record


                // gRooTracker->Fill();

		// Here we expect the FORTRAN end of record padding,
		// read and verify its value.
		if(!input->read(buf.ch, 4)) {
			break;
		}

		if(buf.in[0] != reclen) {
			throw std::runtime_error("Error: unexpected FORTRAN record end padding");
		}

	} // loop over records

	return 0;
}
