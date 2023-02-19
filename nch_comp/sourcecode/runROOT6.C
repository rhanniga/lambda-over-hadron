// Macro for analisys task of preselected-central-diffractive events
//------------------------------------------------------------------
// When there is no time to wait for legotrain results,
// please use this macro to run jobs on grid.
// Recommended for pp or pA runs, not for AA runs.
//------------------------------------------------------------------
// Author: Beomkyu Kim
// email:  kimb@cern.ch
//

//pp 0.9 TeV

// pp 7 TeV 33 runs
const int LHC10bRuns[] = {117222, 117220, 117116, 117112, 117099, 117063, 117060, 117059, 117053, 117052, 117050, 117048, 116645, 116643, 116574, 116571, 116562, 116403, 116402, 115521, 115401, 115399, 115393, 115345, 115335, 115193, 115186, 114931, 114930, 114924, 114918, 114798, 114786};
// pp  7 TeV 31 runs
const int LHC10cRuns[] = {120829,120825,120824,120823,120822,120821,120758,120750,120741,120671,120617,120616,120505,120504,120503,120244,120079,120076,120073,120072,120069,120067,119862,119859,119856,119853,119849,119846,119845,119844,119842,119841};

const int LHC11aRuns[] = {146805};
// pp  7 TeV 55 runs
//V0M is not flat
//const int LHC10dRuns[]= { 124187, 124355, 124358, 124360, 124362, 124367, 124371, 124374, 124378, 124380, 124381, 124383, 124385, 124388, 124600, 124604, 124605, 124606, 125083, 125085, 125097, 125100, 125101, 125133, 125134, 125139, 125140, 125156, 125186, 125295, 125296, 125628, 125630, 125632, 125633, 125842, 125843, 125844, 125847, 125848, 125849, 125850, 125851, 125855, 126004, 126007, 126008, 126073, 126078, 126081, 126082, 126086, 126088, 126090, 126097, 126158, 126160, 126167, 126168, 126177, 126283, 126284, 126285, 126351, 126352, 126359, 126403, 126404, 126405, 126406, 126407, 126408, 126409, 126422, 126424, 126425, 126432, 126437};

//126437 125133 125156 removed because V0M is not flat
//const int LHC10dRuns[] = {125083, 125085, 125097, 125100, 125139, 125140, 125186, 125295, 125296, 125628, 125630, 125632, 125633, 125842, 125843, 125844, 125847, 125848, 125849, 125850, 125851, 125855, 126004, 126007, 126008, 126073, 126078, 126081, 126082, 126086, 126088, 126090, 126097, 126158, 126160, 126167, 126168, 126177, 126283, 126284, 126285, 126351, 126352, 126359, 126403, 126404, 126405, 126406, 126407, 126408, 126409, 126422, 126424, 126425, 126432};
const int LHC10dRuns[] = {126407}; //for diffraction study

const int LHC10dRunsFlat[] = {126404, 126284, 125085};

const int LHC11dRuns[] = {159535,159532,159318,159286,158790,158788,158784,158781,158780,158779,158729,158722,158717,158533,158526,158520,158516,158492,158304,158303,158301,158293,158288,158287,158285,158201,158200,158196,158192,158135,158118,158115,158112,158086,157976,157975,157819,157818,157766,157569,157567,157564,157562,157560,157496,157475,157277,157275,157262,157261,157257,157227,157220,157203,157003,156896,156891,156889,156829,156797,156794,156629,156626,156620};

// pp  7 TeV 113 runs
const int LHC10eRuns[] = {130850, 130848, 130847, 130844, 130842, 130840, 130834, 130804, 130803, 130802, 130799, 130798, 130795, 130793, 130704, 130696, 130628, 130623, 130621, 130620, 130609, 130608, 130601, 130524, 130520, 130519, 130517, 130481, 130480, 130479, 130375, 130356, 130354, 130343, 130342, 130178, 130172, 130157, 130151, 130149, 129966, 129962, 129961, 129960, 129959, 129744, 129742, 129738, 129736, 129735, 129734, 129729, 129726, 129725, 129723, 129659, 129653, 129652, 129651, 129650, 129647, 129641, 129639, 129587, 129586, 129540, 129536, 129528, 129527, 129525, 129524, 129523, 129521, 129520, 129514, 129513, 129512, 129042, 128913, 128855, 128853, 128850, 128843, 128836, 128835, 128834, 128833, 128824, 128823, 128820, 128819, 128778, 128777, 128678, 128677, 128621, 128615, 128611, 128609, 128605, 128596, 128594, 128592, 128590, 128582, 128506, 128504, 128503, 128498, 128495, 128494, 128486, 128366};
// pp 13 TeV  64 runs
//const int LHC15fRuns[]={225011,225016,225026,225031,225035,225037,225041,225043,225050,225051,225052,225105,225106,225305,225307,225310,225313,225314,225315,225322,225576,225578,225579,225580,225582,225586,225587,225589,225709,225710,225716,225717,225763,225768,226062,226085,226170,226175,226176,226177,226183,226210,226225,226444,226445,226466,226468,226472,226495,226500,226532,226543,226551,226554,226569,226573,226591,226593,226596,226600,226602,226603,226605,226606};
const int LHC15fRuns[] = {225580, 225582, 226500, 226532, 226605, 226606};

//const int LHC15fRuns[]={225106,225322,225310,225717,225708,225710,225707,225705,225582,226495,226483,226220};

//pA 5.02 TeV
const int LHC13cRuns[] = {195677, 195675, 195673, 195644, 195635, 195633, 195596, 195593, 195592, 195568, 195567, 195566, 195531, 195529};

//AA 2.76 TeV 87 runs
const int LHC10hRuns[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137243, 137236, 137235, 137232, 137231, 137162, 137161};
//AA 2.76 TeV 68 runs
const int LHC11hRuns[] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};

const int LHC12dRuns[] = {184127, 184132, 184135, 184137, 184138};

//const int LHC15nRuns[]={244364, 244359, 244355, 244351, 244340};

const int LHC15nRunsHighMul[] = {244626, 244627, 244628};																																   // for test
//const int LHC15nRuns[] = {244411, 244416, 244418, 244421, 244453, 244456, 244480, 244481, 244482, 244483, 244484, 244531, 244540, 244542, 244617, 244618, 244619, 244626, 244627, 244628}; // yw: add runs
const int LHC15nRuns[]={244411,244416,244418,244421,244453,244456,244480,244481,244482,244483,244484,244531,244540,244542,244617,244618,244619,244626,244627,244628}; // yw: add runs
// 16l 259395 was removed because of non uniform phi dist of hybrid tracks
// 16l 260010 was removed because of big difference of ITS pure SA tracks to other runs
const int LHC16lRuns[] = {259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259164, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962};

const int LHC17jRuns[] = {274671, 274669, 274653, 274657, 274596, 274595, 274594, 274593};

//pp 5.02 TeV 
//const int LHC17pRuns[] = {282343, 282342, 282341, 282340, 282314, 282313, 282312, 282309, 282307, 282306, 282305, 282304, 282303, 282302, 282247, 282230, 282229, 282227, 282224, 282206, 282189, 282147, 282146, 282127, 282126, 282125, 282123, 282122, 282120, 282119, 282118, 282099, 282098, 282078, 282051, 282050, 282031, 282025, 282021, 282017, 282008};
const int LHC17pRuns[] = {282343, 282342, 282341, 282340, 282314, 282313, 282312, 282309, 282307, 282306};
//const int LHC17pRuns[] = {282343};

const int LHC12hRuns[] = {189122, 189146, 189228, 189229, 189231, 189306, 189310, 189350, 189351, 189352, 189353, 189400, 189407, 189409, 189410, 189411, 18957/7, 189602, 189603, 189605, 189610, 189611, 189612, 189616, 189621, 189623, 189650, 189654, 189656, 189658, 189659, 189696, 189697, 189698, 190212, 190214, 190215, 190216, 190337, 190340, 190341, 190342, 190386, 190389, 190390, 190392, 190393, 190417, 190418, 190419, 190421, 190422, 190425, 190898, 190904, 190979, 190983, 192004, 192073, 192095, 192140, 192172, 192199, 192201, 192205, 192246, 192347, 192348, 192349, 192453, 192468, 192471, 192492, 192499, 192732};

const int LHC15iRuns[] = {236557, 236459, 236453, 236444, 236443, 236441, 236395, 236393, 236389, 236386, 236360, 236357, 236354, 236348, 236337, 236281, 236248, 236246, 236244, 236242, 236238, 236234, 236227, 236222, 236203, 236163, 236161, 236159};

const int LHC16kRuns[] = {258537, 258499, 258477, 258456, 258454, 258426, 258393, 258387, 258359, 258336, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107};
const int LHC17oRuns[] = {281961, 281956, 281953, 281940, 281939, 281932, 281931, 281928, 281920, 281918, 281916, 281915, 281895, 281894, 281893, 281892, 281633, 281592, 281583, 281574, 281569, 281568, 281563, 281562, 281557, 281511, 281509, 281477, 281475, 281450, 281449, 281446, 281444, 281443, 281441, 281415, 281321, 281301, 281277, 281275, 281273, 281271, 281244, 281243, 281242, 281241, 281240, 281213, 281212, 281191, 281190, 281189, 281181, 281180, 281179, 281081, 281080, 281062, 281061, 281060, 281036, 281035, 281033, 281032, 280999, 280998, 280997, 280996, 280994, 280990, 280947, 280943, 280940, 280936, 280897, 280880, 280856, 280854, 280849, 280848, 280847, 280844, 280842, 280793, 280792, 280787, 280786, 280768, 280767, 280766, 280765, 280764, 280763, 280762, 280761, 280757, 280756, 280755, 280754, 280753, 280729, 280706, 280705, 280681, 280679, 280671, 280647, 280645, 280639, 280637, 280636, 280634, 280613, 280583, 280581, 280574, 280551, 280550, 280547, 280546, 280519, 280518, 280499, 280490, 280448, 280447, 280446, 280445, 280443, 280419, 280415, 280412, 280406, 280405, 280403, 280375, 280374, 280351, 280350, 280349, 280348, 280312, 280310, 280290, 280286, 280285, 280284, 280282};
const int LHC17hRuns[] = {273103, 273100, 273099, 273077, 273010, 273009, 272985, 272983, 272976, 272949, 272947, 272939, 272935, 272934, 272933, 272932, 272905, 272903, 272880, 272873,
						  272871, 272870, 272836, 272834, 272833, 272829, 272828, 272784, 272783, 272782,
						  272764, 272763, 272760, 272749, 272747, 272712, 272691, 272690, 272620, 272610,
						  272608, 272607, 272585, 272577, 272575, 272574, 272521, 272468, 272466, 272463,
						  272462, 272461, 272413, 272411, 272400, 272399, 272395, 272394, 272389, 272388,
						  272360, 272359, 272340, 272335, 272194, 272156, 272155, 272154, 272153, 272152,
						  272151, 272123, 272101, 272100, 272076, 272042, 272040, 272039, 272038, 272036,
						  272020, 272018, 271886, 271880, 271874, 271873, 271871, 271870};
const int LHC16oRuns[] = {264035, 264033, 263985, 263984, 263981, 263978, 263977, 263923, 263920, 263917, 263916, 263905, 263866, 263863, 263810, 263803, 263793, 263792, 263790, 263787, 263786, 263785, 263784, 263744, 263743, 263741, 263739, 263738, 263737, 263691, 263690, 263682, 263663, 263662, 263657, 263654, 263652, 263647, 263529, 263497, 263496, 263490, 263487, 263332, 263331, 262858, 262855, 262853, 262849, 262847, 262844, 262842, 262841, 262778, 262777, 262776, 262768, 262760, 262727, 262725, 262723, 262719, 262717, 262713, 262708, 262706, 262705, 262428, 262426, 262425, 262424};

const int LHC16pRuns[] = {264347, 264346, 264345, 264341, 264336, 264312, 264306, 264305, 264281, 264279, 264277, 264273, 264267, 264266, 264265, 264264, 264262, 264261, 264260, 264259, 264238, 264235, 264233, 264232, 264198, 264197, 264194, 264190, 264188, 264168, 264164, 264139, 264138, 264137, 264129, 264110, 264109, 264086, 264085, 264082, 264078, 264076};

//const int LHC16qRuns[] = {265309, 265332, 265334, 265335, 265336, 265338, 265339, 265342, 265343, 265344, 265377, 265378, 265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525};

const int LHC16qRuns[] = {265309, 265332, 265334, 265335, 265336, 265338, 265339, 265377, 265378, 265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525};


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
	const char *taskname = "Mul", const char *option = "LHC16q" // when scanning AOD, add "AOD"
	,
	const char *gridmode = "full" // or "terminate" to merge, or "test" to test
	,
	const char *mode = "grid") //or local
{

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
	plugin->SetMergeViaJDL(1);

	bool isaa = kFALSE;

	plugin->SetSplitMaxInputFileNumber(100);

	if (foption.Contains("LHC10b"))
	{
		plugin->SetGridDataDir("/alice/data/2010/LHC10b/");
		for (int i = 0; i < 33; i++)
			plugin->AddRunNumber(LHC10bRuns[i]);
		plugin->SetDataPattern("pass4/*/*ESDs.root");
	}
	if (foption.Contains("LHC11a"))
	{
		plugin->SetGridDataDir("/alice/data/2011/LHC11a/");
		for (int i = 0; i < 1; i++)
			plugin->AddRunNumber(LHC11aRuns[i]);
		plugin->SetDataPattern("/pass4_with_SDD/*/*ESDs.root");
	}

	if (foption.Contains("LHC10c"))
	{
		plugin->SetGridDataDir("/alice/data/2010/LHC10c/");
		for (int i = 0; i < sizeof(LHC10cRuns) / sizeof(LHC10cRuns[0]); i++)
			plugin->AddRunNumber(LHC10cRuns[i]);
		plugin->SetDataPattern("pass4/*/*ESDs.root");
		if (foption.Contains("MCPYTHIA6"))
		{
			plugin->SetGridDataDir("/alice/sim/2014/LHC14j4c");
			plugin->SetDataPattern("*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(400);
		}
	}
	if (foption.Contains("LHC10d"))
	{
		plugin->SetGridDataDir("/alice/data/2010/LHC10d/");
		//for (int i=0; i<55; i++) plugin->AddRunNumber(LHC10dRuns[i]);
		for (int i = 0; i < sizeof(LHC10dRuns) / sizeof(LHC10dRuns[0]); i++)
			plugin->AddRunNumber(LHC10dRuns[i]);
		//for (int i=0; i<20; i++) plugin->AddRunNumber(LHC10dRuns[i]);
		plugin->SetDataPattern("pass4/*/*ESDs.root");
		if (foption.Contains("MCPYTHIA6"))
		{
			//plugin->SetGridDataDir("/alice/sim/2014/LHC14j4d");
			plugin->SetGridDataDir("/alice/sim/2011/LHC11h4a");
			plugin->SetDataPattern("*/AliESDs.root");
			//plugin->SetSplitMaxInputFileNumber(400);
			plugin->SetSplitMaxInputFileNumber(100);
		}
	}
	if (foption.Contains("LHC11d"))
	{
		plugin->SetGridDataDir("/alice/data/2011/LHC11d/");
		//for (int i=0; i<55; i++) plugin->AddRunNumber(LHC10dRuns[i]);
		for (int i = 0; i < sizeof(LHC11dRuns) / sizeof(LHC11dRuns[0]); i++)
			plugin->AddRunNumber(LHC11dRuns[i]);
		//for (int i=0; i<20; i++) plugin->AddRunNumber(LHC10dRuns[i]);
		plugin->SetDataPattern("pass2_muon/*/*ESDs.root");
		if (foption.Contains("MCPYTHIA6"))
		{
			plugin->SetGridDataDir("/alice/sim/2014/LHC14j4d");
			//plugin->SetGridDataDir("/alice/sim/2015/LHC15g6d/1/");
			if (foption.Contains("Flat"))
				plugin->SetGridDataDir("/alice/sim/2014/LHC14j2a"); // flat
			plugin->SetDataPattern("*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(400);
		}
	}
	if (foption.Contains("LHC12h"))
	{
		plugin->SetGridDataDir("/alice/data/2012/LHC12h/");
		//for (int i=0; i<55; i++) plugin->AddRunNumber(LHC10dRuns[i]);
		for (int i = 0; i < sizeof(LHC12hRuns) / sizeof(LHC12hRuns[0]); i++) plugin->AddRunNumber(LHC12hRuns[i]);
		plugin->SetDataPattern("/pass2/*/*ESDs.root");
		plugin->SetSplitMaxInputFileNumber(100);
		if (foption.Contains("PYTHIA8"))
		{
			plugin->SetGridDataDir("/alice/sim/2015/LHC15h1h");
			plugin->SetDataPattern("*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(300);
		}
	}

	if (foption.Contains("LHC10e"))
	{
		plugin->SetGridDataDir("/alice/data/2010/LHC10e/");
		//in total 113
		for (int i = 60; i < 113; i++)
			plugin->AddRunNumber(LHC10eRuns[i]);
		plugin->SetDataPattern("pass4/*/*ESDs.root");
		if (foption.Contains("MCPYTHIA6"))
		{
			plugin->SetGridDataDir("/alice/sim/2014/LHC14j4e");
			//plugin->SetGridDataDir("/alice/sim/2015/LHC15g6d/1/");
			plugin->SetDataPattern("*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(400);
		}
	}
	if (foption.Contains("LHC13c"))
	{
		plugin->SetGridDataDir("/alice/data/2013/LHC13c/");
		for (int i = 0; i < sizeof(LHC13cRuns) / sizeof(LHC13cRuns[0]); i++)
			plugin->AddRunNumber(LHC13cRuns[i]);
		plugin->SetDataPattern("pass4/*/*ESDs.root");
	}
	if (foption.Contains("LHC16q"))
	{
		plugin->SetGridDataDir("/alice/data/2016/LHC16q/");
		//for (int i = 0; i < sizeof(LHC16qRuns) / sizeof(LHC16qRuns[0]); i++)
		for (int i = 0; i < 10; i++)
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
	if (foption.Contains("LHC11h"))
	{ //68 runs
		plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/");
		//for (int i=0; i<sizeof(LHC11hRuns)/sizeof(LHC11hRuns[0]); i++)
		for (int i = 30; i < 45; i++)
			plugin->AddRunNumber(LHC11hRuns[i]);
		//for (int i=0; i<10; i++) plugin->AddRunNumber(LHC10bRuns[i]);
		//plugin->SetDataPattern("pass2/*/*ESDs.root");
		plugin->SetDataPattern("pass2/AOD115/*/AliAOD.root");
		isaa = kTRUE;
	}
	if (foption.Contains("LHC10h"))
	{ // 87 runs
		plugin->SetGridDataDir("/alice/data/2010/LHC10h/");
		//for (int i=istart; i<iend; i++) {
		//if (i>=sizeof(LHC10hRuns)/sizeof(LHC10hRuns[0])) continue;
		//plugin->AddRunNumber(LHC10hRuns[i]);
		//}
		//plugin->SetDataPattern("pass2/*/*ESDs.root");
		//for (int i=16; i<20; i++) plugin->AddRunNumber(LHC10hRuns[i]);
		for (int i = 0; i < 1; i++)
			plugin->AddRunNumber(LHC10hRuns[i]);
		plugin->SetDataPattern("pass2/AOD160/*/AliAOD.root");
		isaa = kTRUE;
	}
	if (foption.Contains("LHC12d"))
	{
		plugin->SetGridDataDir("/alice/data/2012/LHC12d/");
		//for (int i=0; i<sizeof(LHC11hRuns)/sizeof(LHC11hRuns[0]); i++)
		for (int i = 0; i < 5; i++)
			plugin->AddRunNumber(LHC12dRuns[i]);
		//for (int i=0; i<10; i++) plugin->AddRunNumber(LHC10bRuns[i]);
		plugin->SetDataPattern("pass2/*/*ESDs.root");
	}
	if (foption.Contains("LHC15n"))
	{
		plugin->SetGridDataDir("/alice/data/2015/LHC15n/");
		for (int i = 0; i < sizeof(LHC15nRuns) / sizeof(LHC15nRuns[0]); i++)
			plugin->AddRunNumber(LHC15nRuns[i]);
		plugin->SetDataPattern("/pass4/*/AliESDs.root");
		if (foption.Contains("PYTHIA6"))
		{
			plugin->SetGridDataDir("/alice/sim/2016/LHC16h8b/");
			plugin->SetSplitMaxInputFileNumber(100);
		}
		if (foption.Contains("PYTHIA8"))
		{
			//plugin->SetGridDataDir("/alice/sim/2018/LHC18j3/");
			plugin->SetGridDataDir("/alice/sim/2017/LHC17e2/");
			//plugin->SetGridDataDir("/alice/sim/2016/LHC16h8a/");
			plugin->SetSplitMaxInputFileNumber(100);
			//plugin->SetSplitMaxInputFileNumber(400);
		}
		if (foption.Contains("MC"))
			plugin->SetDataPattern("/*/AliESDs.root");
	}
	if (foption.Contains("LHC15f"))
	{
		plugin->SetGridDataDir("/alice/data/2015/LHC15f/");
		for (int i = 0; i < sizeof(LHC15fRuns) / sizeof(LHC15fRuns[0]); i++)
			plugin->AddRunNumber(LHC15fRuns[i]);
		plugin->SetDataPattern("/pass2/*/AliESDs.root");
		if (foption.Contains("EPOS"))
			plugin->SetGridDataDir("/alice/sim/2016/LHC16d3/");
		if (foption.Contains("PYTHIA6"))
		{
			plugin->SetGridDataDir("/alice/sim/2015/LHC15g3c3/");
			plugin->SetSplitMaxInputFileNumber(400);
		}
		if (foption.Contains("PYTHIA8"))
		{
			plugin->SetGridDataDir("/alice/sim/2015/LHC15g3a3/");
			plugin->SetSplitMaxInputFileNumber(100);
		}
		if (foption.Contains("MC"))
			plugin->SetDataPattern("*/AliESDs.root");
	}
	if (foption.Contains("LHC17j"))
	{
		plugin->SetGridDataDir("/alice/data/2017/LHC17j/");
		for (int i = 0; i < sizeof(LHC17jRuns) / sizeof(LHC17jRuns[0]); i++)
			plugin->AddRunNumber(LHC17jRuns[i]);
		plugin->SetDataPattern("/pass1/*/AliESDs.root");
		if (foption.Contains("EPOS"))
			plugin->SetGridDataDir("/alice/sim/2016/LHC16d3/");
		if (foption.Contains("PYTHIA8"))
			plugin->SetGridDataDir("/alice/sim/2017/LHC17h11/");
		if (foption.Contains("PYTHIA6"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17h7a/");
		}
		if (foption.Contains("PHOJET"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17h7b/");
		}
		if (foption.Contains("MC"))
			plugin->SetDataPattern("*/AliESDs.root");
	}
	if (foption.Contains("LHC17h"))
	{
		plugin->SetGridDataDir("/alice/data/2017/LHC17h/");
		for (int i = 0; i < sizeof(LHC17hRuns) / sizeof(LHC17hRuns[0]); i++)
			plugin->AddRunNumber(LHC17hRuns[i]);
		plugin->SetDataPattern("/pass1/*/AliESDs.root");
		plugin->SetSplitMaxInputFileNumber(200);
	}
	if (foption.Contains("LHC16o"))
	{
		plugin->SetGridDataDir("/alice/data/2016/LHC16o/");
		for (int i = 0; i < sizeof(LHC16oRuns) / sizeof(LHC16oRuns[0]); i++)
			plugin->AddRunNumber(LHC16oRuns[i]);
		plugin->SetDataPattern("/pass1/*/AliESDs.root");
		plugin->SetSplitMaxInputFileNumber(200);
		if (foption.Contains("PYTHIA8"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17d16_extra2/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(100);
		}
	}
	if (foption.Contains("LHC16p"))
	{
		plugin->SetGridDataDir("/alice/data/2016/LHC16p/");
		for (int i = 0; i < sizeof(LHC16pRuns) / sizeof(LHC16pRuns[0]); i++)
			plugin->AddRunNumber(LHC16pRuns[i]);
		plugin->SetDataPattern("/pass1/*/AliESDs.root");
		plugin->SetSplitMaxInputFileNumber(200);
		if (foption.Contains("PYTHIA8"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17d18_extra2/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(100);
		}
	}

	if (foption.Contains("LHC15o"))
	{
		plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
		plugin->AddRunNumber(245064);
		plugin->SetDataPattern("/pass_lowint_firstphys/*/AliESDs.root");
		plugin->SetSplitMaxInputFileNumber(20);
		if (foption.Contains("HIJING"))
			plugin->SetGridDataDir("/alice/sim/2015/LHC15k1a1_plus/");
		if (foption.Contains("MC"))
		{
			plugin->SetDataPattern("*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(50);
		}
		isaa = true;
	}

	if (foption.Contains("LHC17n"))
	{
		plugin->SetGridDataDir("/alice/data/2017/LHC17n/");
		plugin->AddRunNumber(280235);
		plugin->SetDataPattern("/pass1_itsrecpoints/*/AliESDs.root");
		if (foption.Contains("Tracks"))
			plugin->SetDataPattern("/pass1/*/AliESDs.root");
		plugin->SetSplitMaxInputFileNumber(10);
		if (foption.Contains("HIJING"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17j6/");
			if (foption.Contains("Tracks"))
				plugin->SetGridDataDir("/alice/sim/2017/LHC17j7/");
		}
		if (foption.Contains("MC"))
		{
			plugin->SetDataPattern("*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(20);
		}
		isaa = true;
	}

	if (foption.Contains("LHC16l"))
	{
		plugin->SetGridDataDir("/alice/data/2016/LHC16l/");
		plugin->SetDataPattern("/pass2/*/AliESDs.root");
		//		plugin->SetDataPattern("/pass2/AOD234/*/AliAOD.root");
		//plugin->SetDataPattern("/AOD208/*/AliAOD.root");
		plugin->SetSplitMaxInputFileNumber(100);
		if (foption.Contains("PYTHIA8"))
		{
		  plugin->SetGridDataDir("/alice/sim/2018/LHC18d8/");
		  plugin->SetDataPattern("/*/AliESDs.root");
		  //		  plugin->SetGridDataDir("/alice/sim/2019/LHC19h11c/");
		  //		  plugin->SetDataPattern("/*/AliESDs.root");
		  //  plugin->SetDataPattern("/*/AliAOD.root");
			plugin->SetSplitMaxInputFileNumber(50
			);
		}
		if (foption.Contains("EPOS"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17d20b2/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(100);
		}

		//		for (int i = 0; i < sizeof(LHC16lRuns) / sizeof(LHC16lRuns[0]); i++)
		for (int i = 0; i < sizeof(LHC16lRuns) / sizeof(LHC16lRuns[0]); i++)
		//		for (int i = 0; i < 1; i++)
			plugin->AddRunNumber(LHC16lRuns[i]);
	}
	if (foption.Contains("LHC16k"))
	{
		plugin->SetGridDataDir("/alice/data/2016/LHC16k/");
		plugin->SetDataPattern("/pass2/*/AliESDs.root");
		if (foption.Contains("PYTHIA8"))
		{
			plugin->SetGridDataDir("/alice/sim/2018/LHC18f1/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(300);
		}
		if (foption.Contains("EPOS"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17d20b2/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(100);
		}

		for (int i = 0; i < sizeof(LHC16kRuns) / sizeof(LHC16kRuns[0]); i++)
			plugin->AddRunNumber(LHC16kRuns[i]);
	}

	if (foption.Contains("LHC15i"))
	{
		plugin->SetGridDataDir("/alice/data/2015/LHC15i/");
		plugin->SetDataPattern("/pass2/*/AliESDs.root");
		if (foption.Contains("PYTHIA8"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17i4_2/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(100);
		}

		for (int i = 0; i < sizeof(LHC15iRuns) / sizeof(LHC15iRuns[0]); i++)
			plugin->AddRunNumber(LHC15iRuns[i]);
	}

	if (foption.Contains("LHC17p"))
	{
		plugin->SetGridDataDir("/alice/data/2017/LHC17p/");
		plugin->SetDataPattern("/pass1_CENT_woSDD/*/AliESDs.root");
		plugin->SetSplitMaxInputFileNumber(200);
		if (foption.Contains("PYTHIA8"))
		{
		  //			plugin->SetGridDataDir("/alice/sim/2017/LHC17l4_cent_lastG4fix/");
			plugin->SetGridDataDir("/alice/sim/2017/LHC17l3b_cent_woSDD/");
			//plugin->SetGridDataDir("/alice/sim/2018/LHC18j2_fast/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(900);
		}
		if (foption.Contains("EPOS"))
		{
			plugin->SetGridDataDir("/alice/sim/2017/LHC17l4b_fast/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(100);
		}

		for (int i = 0; i < sizeof(LHC17pRuns) / sizeof(LHC17pRuns[0]); i++)
			plugin->AddRunNumber(LHC17pRuns[i]);
	}

	if (foption.Contains("LHC17o"))
	{ // for 17p highmult MC
		plugin->SetGridDataDir("/alice/data/2017/LHC17o/");
		plugin->SetDataPattern("/pass1/*/AliESDs.root");
		plugin->SetSplitMaxInputFileNumber(200);
		if (foption.Contains("PYTHIA8"))
		{
			plugin->SetGridDataDir("/alice/sim/2018/LHC18a9/");
			plugin->SetDataPattern("/*/AliESDs.root");
			plugin->SetSplitMaxInputFileNumber(300);
		}
		for (int i = 0; i < sizeof(LHC17oRuns) / sizeof(LHC17oRuns[0]); i++)
			plugin->AddRunNumber(LHC17oRuns[i]);
	}

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
