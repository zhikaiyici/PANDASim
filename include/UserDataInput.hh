#ifndef UserDataInput_h
#define UserDataInput_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include <vector>

using namespace std;

class UserDataInput
{
public:
	UserDataInput();
	~UserDataInput();

	static void ReadInputData(); //读取外部输入文件

	static inline G4int GetNumberOfThreads() { return threadNumber; }
	static inline G4int GetNumberOfEvents() { return eventNumber; }
	static inline G4int GetSizeOfArray() { return arraySize; }
	static inline G4int GetPositionOfNeutrino() { return neutrinoPosition; }

	static inline G4bool GetStatusOfUI() { return uiStatus; }
	static inline G4bool GetOpticalPhysicsStatus() { return opticalPhysics; }

	static inline G4double GetDistanceBetweenModules() { return distanceBetweenModules * cm; }
	static inline G4double GetDectorDimensionX() { return dtctrX * cm; }
	static inline G4double GetDectorDimensionY() { return dtctrY * cm; }
	static inline G4double GetDectorDimensionZ() { return dtctrZ * cm; }
	static inline G4double GetNeutrinoPercentage() { return neutrinoPercentage; }
	static inline G4double GetGdFilmThickness() { return gdFilmThickness * um; }

	static inline G4String GetSoureType() { return sourceType; }
	static inline G4String GetSourePosition() { return sourcePosition; }

	static inline vector<G4double> GetPositronEnergy() { return positronEnergy; }
	static inline vector<G4double> GetPositronCDFSpectrum() { return positronCDFSpectrum; }
	static inline vector<G4double> GetNeutronEnergy() { return neutronEnergy; }
	static inline vector<G4double> GetNeutronCDFSpectrum() { return neutronCDFSpectrum; }

private:

	static G4int threadNumber;
	static G4int eventNumber; //要模拟的粒子数
	static G4int arraySize;
	static G4int neutrinoPosition;

	static G4bool uiStatus;
	static G4bool opticalPhysics; //是否打开光学过程

	static G4double distanceBetweenModules; //模块间的空隙
	static G4double dtctrX, dtctrY, dtctrZ; //探测器的长宽高
	static G4double neutrinoPercentage; //总谱中neutrino事件占比
	static G4double gdFilmThickness; //Gd膜厚度

	static G4String sourceType;
	static G4String sourcePosition;

	static vector<G4double> positronEnergy; // 正电子能量
	static vector<G4double> positronCDFSpectrum; //正电子归一化累计能谱
	static vector<G4double> neutronEnergy; //中子能量
	static vector<G4double> neutronCDFSpectrum; //中子归一化累计能谱

	static void ReadSpectra(G4String spectrumName, vector<G4double>& energy, vector<G4double>& cdfSpectrum); //读取能谱文件
};

#endif // !UserDataInput_h
