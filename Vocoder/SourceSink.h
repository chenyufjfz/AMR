/*
 * SourceSink.h
 *
 *  Created on: Dec 4, 2015
 *      Author: yuchen
 */

#ifndef SOURCESINK_H_
#define SOURCESINK_H_
#include <string>
#include <vector>
#include "public.h"
using namespace std;

enum AudioFormat{
	RAW,
	WAV
};
typedef   struct   {
	char         fccID[4];       /*   should   be   "RIFF "   */
	unsigned   long       dwSize;           /*   byte_number   behind   it   */
	char         fccType[4];   /*   should   be   "WAVE "   */
} WAVEHEADER;
typedef   struct   {
	char         fccID[4];       /*   should   be   "fmt   "   */
	unsigned   long       dwSize;           /*   should   be   0x10   */
	unsigned   short     wFormatTag;   /*   should   be   1   */
	unsigned   short     wChannels;
	unsigned   long       dwSamplesPerSec;
	unsigned   long       dwAvgBytesPerSec;
	unsigned   short     wBlockAlign;
	unsigned   short     uiBitsPerSample;
} WAVEFMT;
typedef   struct   {
	char         fccID[4];         /*   should   be   "data "   */
	unsigned   long       dwSize;             /*   byte_number   of   PCM   data   in   byte*/
} WAVEDATA;
class AudioSource {
private:
	FILE * fp;
	string filename;
	AudioFormat fmt;
	int channel_num;
	int sample_rate;
	int byte_num;
	vector<int> buf;
public:
	AudioSource() {
		fp =NULL;
	}
	~AudioSource() {
		if (fp!=NULL)
			close();
	}
	void open(string source, AudioFormat format=RAW, int channel=-1, int sr=-1, int bytes=2) {
		filename = source;
		fmt = format;
		channel_num = channel;
		byte_num = bytes;
		sample_rate = sr;
		fp = fopen(filename.c_str(), "rb");
		if (fp==NULL)
			throw MyException("source file" + filename +" not exit");
		if (format==WAV) {
			WAVEHEADER header;
			WAVEFMT format;
			WAVEDATA data;
			fread(&header, sizeof(header), 1, fp);
			fread(&format, sizeof(format), 1, fp);
			fread(&data, sizeof(data), 1, fp);
			if (memcmp(header.fccID, "RIFF", 4)!=0 || memcmp(header.fccType, "WAVE", 4)!=0 ||
				memcmp(format.fccID, "fmt ", 4)!=0 || memcmp(data.fccID, "data", 4)!=0 ||
				format.dwSize!=0x10 || format.wFormatTag!=1)
				throw MyException("not wav format");
			if ((int)format.wChannels!=channel && channel!=-1)
				throw MyException("wav channel number mismatch");
			channel_num = format.wChannels;
			if ((int)format.dwSamplesPerSec !=sr && sr!=-1)
				throw MyException("wav sample rate mismatch");
			sample_rate = format.dwSamplesPerSec;
			if (((int)format.dwAvgBytesPerSec != sr * bytes * channel || (int) format.wBlockAlign!=bytes
					|| (int) format.uiBitsPerSample!=bytes*8) && bytes!=-1)
				throw MyException("wav PCM byte mismatch");
			byte_num = format.wBlockAlign;
		}
		buf.clear();
	}
	void close() {
		fclose(fp);
	}

	int getpcm(vector<int> & v, int pcm_num) {
		if (fmt==RAW || fmt==WAV) {
			void * buf;
			int pcm_read;
			buf =malloc(pcm_num *byte_num * channel_num);
			pcm_read=fread(buf, byte_num* channel_num, pcm_num, fp);
			v.assign(pcm_num * channel_num, 0);
			if (byte_num==2) {
				short * buf_16 = (short*) buf;
				for (int i=0; i<pcm_read * channel_num; i++)
					v[i] = buf_16[i];
			} else
			if (byte_num==4) {
				int * buf_32 = (int*) buf;
				for (int i=0; i<pcm_read * channel_num; i++)
					v[i] = buf_32[i];
			}
			free(buf);
			return pcm_read;
		}
		return -1;
	}
	bool get_bufpcm(vector<int> & v, int pcm_num, int buf_num) {
		vector<int> s;
		int buf_size, i;
		ASSERT(buf_num<pcm_num && (int) buf.size() <pcm_num);
		if (feof(fp))
			return false;
		v.resize(pcm_num);
		buf_size = buf.size();
		for (i=0; i<buf_size; i++)
			v[i] = buf[i];
		getpcm(s, pcm_num-buf_size);
		for (i=buf_size; i<pcm_num; i++)
			v[i] = s[i-buf_size];
		buf.resize(buf_num);
		for (i=0; i<buf_num; i++)
			buf[i] = v[pcm_num-buf_num+i];
		return true;
	}
	int get_samplerate() {
		return sample_rate;
	}
};

class AudioSink {
private:
	FILE * fp;
	string filename;
	AudioFormat fmt;
	int channel_num;
	int byte_num;
	int sample_rate;
	unsigned write_bytes;
	vector<int> buf;
public:
	AudioSink() {
		fp =NULL;
	}
	~AudioSink() {
		if (fp!=NULL)
			close();
	}
	void open(string source, AudioFormat format=RAW, int channel=1, int sr=48000, int bytes=2) {
		filename = source;
		fmt = format;
		channel_num = channel;
		sample_rate = sr;
		byte_num = bytes;
		fp = fopen(filename.c_str(), "wb");
		if (fp==NULL)
			throw MyException("source file" + filename +" not exit");
		buf.clear();
		write_bytes=0;
		if (fmt==WAV) {
			fseek(fp,sizeof(WAVEHEADER) + sizeof(WAVEFMT) + sizeof(WAVEDATA),1);
		}
	}
	void close() {
		if (fmt==WAV) {
			WAVEHEADER header;
			WAVEFMT format;
			WAVEDATA data;
			memcpy(header.fccID, "RIFF", 4);
			memcpy(header.fccType, "WAVE", 4);
			header.dwSize = write_bytes+44;

			memcpy(format.fccID, "fmt ", 4);
			format.wChannels = channel_num;
			format.dwSize = 0x10;
			format.wFormatTag = 1;
			format.dwSamplesPerSec = sample_rate;
			format.dwAvgBytesPerSec = sample_rate * byte_num * channel_num;
			format.wBlockAlign = byte_num;
			format.uiBitsPerSample = byte_num *8;

			memcpy(data.fccID, "data", 4);
			data.dwSize = write_bytes;
			rewind(fp);
			fwrite(&header,sizeof(header),1,fp);
			fwrite(&format,sizeof(format),1,fp);
			fwrite(&data, sizeof(data), 1, fp);
		}
		fclose(fp);
		fp = NULL;
	}
	void putpcm(vector<int> & v, int pcm_num=0) {
		if (pcm_num==0)
			pcm_num = v.size() / channel_num;
		if (fmt==RAW || fmt==WAV) {
			void * buf;
			ASSERT(pcm_num * channel_num <= (int)v.size());
			buf = malloc(pcm_num *byte_num * channel_num);
			if (byte_num==2) {
				short * buf_16 = (short*)buf;
				for (int i=0; i<pcm_num * channel_num; i++)
					buf_16[i] = (v[i]>32767) ? 32767 : (v[i]<-32768 ? -32768 : v[i]);
			}
			if (byte_num==4) {
				int * buf_32 = (int*)buf;
				for (int i=0; i<pcm_num * channel_num; i++)
					buf_32[i] = v[i];
			}
			if ((int) fwrite(buf, byte_num* channel_num, pcm_num, fp) !=pcm_num)
				throw MyException("write bytes i/o error");
			write_bytes+=byte_num* channel_num*pcm_num;
			free(buf);
		}
	}
	void put_bufpcm(vector<int> & v, int buf_num) {
		int i;
		vector<int> s;

		ASSERT(buf_num <(int)v.size() && buf.size() <v.size());
		s.assign(v.size(), 0);
		for (i=0; i<(int) buf.size(); i++)
			s[i] = buf[i];
		for (i=0; i<(int)v.size(); i++)
			s[i] += v[i];
		buf.resize(buf_num);
		for (i=0; i<buf_num; i++)
			buf[i] = s[v.size()-buf_num+i];
		putpcm(s, v.size()-buf_num);
	}
};


#endif /* SOURCESINK_H_ */
