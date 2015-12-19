//============================================================================
// Name        : vocoder.cpp
// Author      : Chenyu
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdio.h>
#include <string>
#include <vector>
#include <complex>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>
#include <algorithm>
#include "public.h"
#include "Filter.h"
#include "LinearEstimate.h"
#include "RandomSeq.h"
#include "Fixed.h"
#include "Resample.h"
#include "fitsutil.h"
#include "mixer.h"
#include "delay.h"
#include "Correlation.h"
#include "SourceSink.h"
#include "FFT.h"
using namespace std;

extern "C" {
#include "typedef.h"
#include "interf_enc.h"
#include "interf_dec.h"
}

typedef complex<double> ComplexDouble;


//pass band ripple =0.06, stop band ripple = 0.0015
const double l_1700_1800[] = {
		0.000830736,0.000346737,-0.000804337,-0.000729008,0.000609432,0.00107665,-0.000241184,-0.00130975,-0.000270373,0.00135325, 0.000859143,-0.00115328,
		-0.0014289,0.00069233,0.00186638,-6.85813e-018,-0.00206015,-0.000843627,0.00192247,0.0017137,-0.00141051,-0.00245588,0.000542707,0.00290982,0.000593542,
		-0.0029377,-0.00184565,0.00245341,0.00301214,-0.00144709,-0.00387036,-4.03886e-018,0.0042125,0.00171435,-0.00388465,-0.0034451,0.00282259,0.00489455,
		-0.00107779,-0.00576136,-0.00117229,0.00579094,0.00363321,-0.00482564,-0.00592321,0.00284667,0.00762125,-1.39233e-017,-0.00832817,-0.0033997,0.00773322,
		0.00689022,-0.00567654,-0.00990746,0.00219806,0.0118514,0.00243523,-0.0121645,-0.0077287,0.0104123,0.0129871,-0.00635529,-0.0173646,1.62303e-017,0.019927,
		0.00837761,-0.0197055,-0.0182429,0.0157061,0.028847,-0.00679315,-0.0393018,-0.0087881,0.0486753,0.0351992,-0.056097,-0.0880693,0.0608592,0.312134,0.4375,
		0.312134,0.0608592,-0.0880693,-0.056097,0.0351992,0.0486753,-0.0087881,-0.0393018,-0.00679315,0.028847,0.0157061,-0.0182429,-0.0197055,0.00837761,0.019927,
		1.62303e-017,-0.0173646,-0.00635529,0.0129871,0.0104123,-0.0077287,-0.0121645,0.00243523,0.0118514,0.00219806,-0.00990746,-0.00567654,0.00689022,0.00773322,
		-0.0033997,-0.00832817,-1.39233e-017,0.00762125,0.00284667,-0.00592321,-0.00482564,0.00363321,0.00579094,-0.00117229,-0.00576136,-0.00107779,0.00489455,
		0.00282259,-0.0034451,-0.00388465,0.00171435,0.0042125,-4.03886e-018,-0.00387036,-0.00144709,0.00301214,0.00245341,-0.00184565,-0.0029377,0.000593542,
		0.00290982,0.000542707,-0.00245588,-0.00141051,0.0017137,0.00192247,-0.000843627,-0.00206015,-6.85813e-018,0.00186638,0.00069233,-0.0014289,-0.00115328,
		0.000859143,0.00135325,-0.000270373,-0.00130975,-0.000241184,0.00107665,0.000609432,-0.000729008,-0.000804337,0.000346737,0.000830736
};
const double l_1800_1900[] = {
		0.00084114,0.000211517,-0.000907576,-0.000468052,0.000912079,0.000756839,-0.00083918,-0.00105961,0.000677175,0.00135325,-0.00041976,-0.0016109,6.74689e-005,
		0.00180357,0.000371247,-0.00190214,-0.000879399,0.00187965,0.00143143,-0.0017137,-0.0019938,0.00138892,0.00252629,-0.000899184,-0.00298394,0.00024948,0.00331953,
		0.000542773,-0.00348666,-0.00144709,0.00344303,0.00242007,-0.00315394,-0.00340648,0.00259564,0.00434108,-0.00175846,-0.00515144,0.000649342,0.00576136,0.000706276,
		-0.00609488,-0.00226347,0.00608066,0.00395776,-0.00565643,-0.0057061,0.00477333,0.00740865,-0.0033997,-0.00895147,0.00152434,0.0102096,0.000841376,-0.0110504,
		-0.00366228,0.011336,0.00687959,-0.0109248,-0.0104123,0.00966991,0.0141599,-0.0074123,-0.0180063,0.00396373,0.0218243,0.000930443,-0.0254817,-0.00767369,
		0.028847,0.0170141,-0.0317958,-0.0305775,0.0342167,0.0526793,-0.0360165,-0.0993732,0.0371254,0.316043,0.4625,0.316043,0.0371254,-0.0993732,-0.0360165,0.0526793,
		0.0342167,-0.0305775,-0.0317958,0.0170141,0.028847,-0.00767369,-0.0254817,0.000930443,0.0218243,0.00396373,-0.0180063,-0.0074123,0.0141599,0.00966991,-0.0104123,
		-0.0109248,0.00687959,0.011336,-0.00366228,-0.0110504,0.000841376,0.0102096,0.00152434,-0.00895147,-0.0033997,0.00740865,0.00477333,-0.0057061,-0.00565643,
		0.00395776,0.00608066,-0.00226347,-0.00609488,0.000706276,0.00576136,0.000649342,-0.00515144,-0.00175846,0.00434108,0.00259564,-0.00340648,-0.00315394,0.00242007,
		0.00344303,-0.00144709,-0.00348666,0.000542773,0.00331953,0.00024948,-0.00298394,-0.000899184,0.00252629,0.00138892,-0.0019938,-0.0017137,0.00143143,0.00187965,
		-0.000879399,-0.00190214,0.000371247,0.00180357,6.74689e-005,-0.0016109,-0.00041976,0.00135325,0.000677175,-0.00105961,-0.00083918,0.000756839,0.000912079,
		-0.000468052,-0.000907576,0.000211517,0.00084114
};
const double l_1900_2000[] = {
		0.000846358,7.10892e-005,-0.000960662,-0.00016128,0.00107587,0.000272047,-0.00118985,-0.000404736,0.00130023,0.000560535,-0.00140437,-0.000740449,0.0014994,
		0.000945276,-0.00158224,-0.00117559,0.00164957,0.00143171,-0.00169785,-0.0017137,0.00172337,0.00202133,-0.00172221,-0.00235409,0.00169026,0.00271118,-0.00162324,
		-0.00309148,0.00151667,0.00349358,-0.00136584,-0.00391576,0.00116584,0.00435603,-0.000911467,-0.00481212,0.000597153,0.00528149,-0.000216893,-0.00576136,-0.00023591,
		0.00624875,0.000768647,-0.00674046,-0.00138978,0.00723317,0.00210924,-0.00772341,-0.00293899,0.00820761,0.00389381,-0.00868219,-0.00499249,0.0091435,0.00625954,
		-0.00958797,-0.00772788,0.010012,0.00944299,-0.0104123,-0.0114697,0.0107855,0.0139039,-0.0111285,-0.0168933,0.0114384,0.0206778,-0.0117126,-0.0256734,0.0119488,
		0.0326684,-0.0121449,-0.0433551,0.0122992,0.0621394,-0.0124104,-0.105186,0.0124776,0.318003,0.4875,0.318003,0.0124776,-0.105186,-0.0124104,0.0621394,0.0122992,
		-0.0433551,-0.0121449,0.0326684,0.0119488,-0.0256734,-0.0117126,0.0206778,0.0114384,-0.0168933,-0.0111285,0.0139039,0.0107855,-0.0114697,-0.0104123,0.00944299,
		0.010012,-0.00772788,-0.00958797,0.00625954,0.0091435,-0.00499249,-0.00868219,0.00389381,0.00820761,-0.00293899,-0.00772341,0.00210924,0.00723317,-0.00138978,
		-0.00674046,0.000768647,0.00624875,-0.00023591,-0.00576136,-0.000216893,0.00528149,0.000597153,-0.00481212,-0.000911467,0.00435603,0.00116584,-0.00391576,
		-0.00136584,0.00349358,0.00151667,-0.00309148,-0.00162324,0.00271118,0.00169026,-0.00235409,-0.00172221,0.00202133,0.00172337,-0.0017137,-0.00169785,0.00143171,
		0.00164957,-0.00117559,-0.00158224,0.000945276,0.0014994,-0.000740449,-0.00140437,0.000560535,0.00130023,-0.000404736,-0.00118985,0.000272047,0.00107587,-0.00016128,
		-0.000960662,7.10892e-005,0.000846358
};
void test_transform(int argc, char **argv)
{
	AudioSource au;
	AudioSink as, as2;
	string filesrc, filedst, filerec;
	AudioFormat in_fmt=WAV;
	int *enstate;
	int *destate;
	unsigned char serial_data[32];
	short speech[160];
	short synth[160];
	enum Mode req_mode = MR122;

	for (int i=1; i<argc; i++) {
		if (argv[i][0]=='-') {
			if (argv[i][1]=='r' || argv[i][1]=='R')
				in_fmt=RAW;
		} else {
			if (filesrc.empty())
				filesrc = argv[i];
			else
				filedst = argv[i];
		}
	}

	if (filesrc.empty())
		filesrc = "C:\\ChenYu\\work\\vocoder\\wav\\rc_8.wav";
		//filesrc = "C:\\ChenYu\\work\\vocoder\\wav\\sweep_100_4000.wav";
	if (filedst.empty())
		filedst = "C:\\ChenYu\\work\\vocoder\\wav\\transform_8.wav";
	filerec = filesrc;
	filerec.replace(filerec.end()-4, filerec.end(), "_r.wav");

	au.open(filesrc.c_str(), in_fmt, 1, 8000, 2);
	as.open(filedst.c_str(), WAV, 1, 8000, 2);
	//as2.open(filerec.c_str(), WAV, 1, 8000, 2);

	vector<double> f1, f2, f3;
	f1.assign(l_1700_1800, l_1700_1800 + sizeof(l_1700_1800)/sizeof(l_1700_1800[0]));
	f2.assign(l_1800_1900, l_1800_1900 + sizeof(l_1800_1900)/sizeof(l_1800_1900[0]));
	f3.assign(l_1900_2000, l_1900_2000 + sizeof(l_1900_2000)/sizeof(l_1900_2000[0]));
	FreqSwap<double> tx_swap(f1, f2), rx_swap(f2, f3);
	vector<int> x, y;
	int frame=0;
	enstate = (int*) Encoder_Interface_init(0);
	destate = (int*) Decoder_Interface_init();
	while (1) {
		if (!au.getpcm(x, 160))
			break;
		
		CCfits::fill(x, speech);
		Encoder_Interface_Encode(enstate, req_mode, speech, serial_data, 0);
		Decoder_Interface_Decode(destate, serial_data, synth, 0);
		CCfits::fill(synth, y, 0, 160);

		as.putpcm(y, 160);
		//as2.putpcm(x, 160);
	}
	as.close();
	au.close();
}


int main(int argc, char **argv)
{
	test_transform(argc, argv);
}