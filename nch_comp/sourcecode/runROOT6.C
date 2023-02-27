
const int LHC16qRuns[] = {265309, 265332, 265334, 265336, 265338, 265339, 265377, 265378, 265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525};


class AliAnalysisGrid;

#include <vector>
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "AliCentralitySelectionTask.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisPseudoRapidityDensityTemp.h"

R__ADD_INCLUDE_PATH($ALICE_PHYSICS/OADB/macros)

#include "AddTaskPhysicsSelection.C"
#include "AddTaskCentrality.C"
R__ADD_INCLUDE_PATH($ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/) 
#include "AddTaskMultSelection.C"
R__ADD_INCLUDE_PATH($ALICE_ROOT/ANALYSIS/macros/)
#include "AddTaskPIDResponse.C"

void runROOT6(
	const char *taskname = "nch_trig_mc"
	, 
	const char *option = "LHC16qMC" // when scanning AOD, add "AOD"
	,
	const char *gridmode = "terminate" // or "terminate" to merge, or "test" to test
	,
	const char *mode = "grid") //or local
{

    bool gridMerge = false;

	gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY -g");
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libSTEERBase. so");
	gSystem->Load("libESD.so");
	gSystem->Load("libAOD.so");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libOADB.so");
	gSystem->Load("libANALYSISalice.so");
	//gSystem->Load("libqpythia.so");
	//gSystem->Load("libAliPythia6.so");
	gSystem->Load("libPWGTools.so");
	gSystem->Load("libpythia6_4_21.so");

	gInterpreter->ProcessLine(".include $ROOTSYS/include");
	gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
	gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");

	// analysis manager
	AliAnalysisManager *mgr = new AliAnalysisManager(Form("%s%smanager", taskname, option));
	TString foption = option;

	// create the alien handler and attach it to the manager
	AliAnalysisAlien *plugin = new AliAnalysisAlien();
	plugin->SetNtestFiles(5); //Chiara addition
	plugin->SetRunMode(gridmode);
	plugin->SetAPIVersion("V1.1x");
	//plugin->SetAliPhysicsVersion("vAN-20181124-1");
	plugin->SetAliPhysicsVersion("vAN-20200910_ROOT6-1");
	plugin->SetDropToShell(0);
	if (!foption.Contains("MC"))
		plugin->SetRunPrefix("000");
	plugin->SetNrunsPerMaster(1);
	plugin->SetOutputToRunNo();
	plugin->SetMergeViaJDL(gridMerge);

	bool isaa = kFALSE;

	plugin->SetSplitMaxInputFileNumber(100);

	if (foption.Contains("LHC16q"))
	{
		plugin->SetGridDataDir("/alice/data/2016/LHC16q/");
		for (int i = 0; i < 17; i++)
			plugin->AddRunNumber(LHC16qRuns[i]);
			plugin->SetDataPattern("pass2_FAST/*/*ESDs.root");
		if (foption.Contains("MC"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17f2b_fast/");
			// plugin->SetGridDataDir("/alice/sim/2015/LHC15g6d/1/");
			plugin->SetDataPattern("*/AliESDs.root");
		}
		isaa = true;
	}

	// plugin->SetGridWorkingDir(Form("CORRECT_NCH_COMPUTATION%s%s", taskname, option));
	plugin->SetGridWorkingDir(Form("CORRECT_NCH_COMPUTATION%s%s", taskname, option));
	plugin->SetGridOutputDir("out");
	plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include   -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS   -I$ALICE_PHYSICS/OADB/macros");
	plugin->SetAnalysisSource("AliAnalysisPseudoRapidityDensityTemp.cxx");
	plugin->SetAdditionalLibs("libpythia6_4_21.so AliAnalysisPseudoRapidityDensityTemp.cxx AliAnalysisPseudoRapidityDensityTemp.h");
	//plugin->SetAdditionalLibs("libqpythia.so libAliPythia6.so libPWGTools.so");

	plugin->SetDefaultOutputs(kFALSE);
	//plugin->SetOutputFiles("AnalysisResults.root RecTree.root");
	plugin->SetOutputFiles("AnalysisResults.root");
	plugin->SetMasterResubmitThreshold(90);
	plugin->SetFileForTestMode("file.text");
	//plugin->SetUseSubmitPolicy();


	// Optionally set time to live (default 30000 sec)
	plugin->SetTTL(20000);
	// Optionally set input format (default xml-single)
	plugin->SetInputFormat("xml-single");
	// Optionally modify the name of the generated JDL (default analysis.jdl)
	plugin->SetJDLName(Form("%s%s.jdl", taskname, option));
	// Optionally modify the executable name (default analysis.sh)
	plugin->SetExecutable(Form("%s%s.sh", taskname, option));
	// Optionally modify job price (default 1)
	plugin->SetPrice(1);
	// Optionally modify split mode (default 'se')
	plugin->SetSplitMode("se");

	mgr->SetGridHandler(plugin);
	AliInputEventHandler *handler = 0x0;
	if (foption.Contains("AOD"))
		handler = new AliAODInputHandler();
	else if (foption.Contains("ITSRec"))
		handler = new AliESDInputHandlerRP();
	else
		handler = new AliESDInputHandler();

	handler->SetNeedField(1);
	mgr->SetInputEventHandler(handler);

	if (foption.Contains("MC"))
	{
		AliMCEventHandler *mchandler = new AliMCEventHandler();
		// Not reading track references
		mchandler->SetReadTR(kFALSE);
		mgr->SetMCtruthEventHandler(mchandler);
	}

	AliPhysicsSelectionTask *physSelTask = 0x0;
	//gROOT->ProcessLine(Form(".include %s","$ALICE_PHYSICS/OADB/macros/"));
	//gInterpreter->Declare("#include \"AddTaskPhysicsSelection.C\"");
	if (foption.Contains("MC"))
		physSelTask = AddTaskPhysicsSelection(1);
	else
		physSelTask = AddTaskPhysicsSelection(0);
	if (!physSelTask)
	{
		Printf("no physSelTask");
		return;
	}

	AliMultSelectionTask *multtask = AddTaskMultSelection(false);
	if (!multtask)
	{
		Printf("no MultSlection");
		return;
	}
	if (foption.Contains("MC") && foption.Contains("LHC16q"))
		multtask->SetAlternateOADBforEstimators("LHC16q");

	gInterpreter->LoadMacro("AliAnalysisPseudoRapidityDensityTemp.cxx+g");

	AliAnalysisPseudoRapidityDensityTemp *task =
		reinterpret_cast<AliAnalysisPseudoRapidityDensityTemp *>(gInterpreter->ExecuteMacro(
			Form("AddTaskTemp.C(\"%s\",\"%s%s\")", taskname, taskname, option)));

	task->SetIsAA(isaa);


	// enable debug printouts
	mgr->SetDebugLevel(2);
	if (!mgr->InitAnalysis())
		return;
	mgr->PrintStatus();

	// start analysis
	Printf("Starting Analysis....");
	mgr->StartAnalysis(mode, 1234567890, 0);
}
