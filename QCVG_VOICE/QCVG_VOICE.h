#ifndef QCVG_VOICE_h
#define QCVG_VOICE_h

#include "Arduino.h"
#include "QCVG_MAIN.h"

extern class QCVG_MAIN QCVG;
extern class QCVG_VOICE *Voices[];
extern int Voicing;
extern int Repetitions;
extern double Glissando;
extern double TimingOffset;
extern int Layers[];
extern boolean SingleVoice;
extern boolean Strata;
extern boolean Debug;
extern boolean RingChanges;
extern boolean ChangeRoot;
extern boolean ShiftRegisterMode;
extern boolean Quantize;
extern boolean Phasing;
extern boolean Modulate;
extern boolean Flicker;
extern boolean BZZZ;
extern boolean ProximityMode;
extern boolean FillIn;
extern boolean ChangeTrigOrder;
extern boolean InputMode;
extern const int TrigPins[];

void Swap(double*, int, int);

class QCVG_VOICE {
public:
	QCVG_VOICE(int, int , char, int);
	void GenerateTriggerSequence(int);
	void PrintTriggerSequence();
	void PrintCVSequence();
	void Advance();
	bool GetTrig(long,bool&);
	void FlipT();
	void Invert();
	void Reverse();
	void OrderSequence();
	void Transpose(double);
	void Shift(int, int);
	double ArrayMin(double[]);
	double ArrayMax(double[]);
	void UpdateNote();
	void sendTrigger();
	bool Play(int);
	int length;
	int TriggerSubdivisions;
	long TriggerStart;
	double TimeFactor;
	int counter;
	long clock;
	long advanceclock;
	int CurrentPosition;
	int iteration;
	double Glide;

private:
	bool *_TrigArray;
	double *_NoteArray;
	char _Rule;
	double _roof;
	double _floor;
	double CV_prev;
	double _Interval;
	int _n;
	bool direction;
	int _channel;
	unsigned long _GlideLength;
	bool _Play;
	void Bjorklund(int,int);
	void ShiftBoolean(int);
	double GetCV(int = -1);
	void FixInterval(double&);
	bool RandomStep(double);
	void SetupSingleVoice();
	int PreviousVoice();
	void EvaluateProximity();
	void RandomizeInitialOrder();
};

#endif