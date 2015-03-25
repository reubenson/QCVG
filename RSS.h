// #ifndef RSS_h
// #define RSS_h

#include "Arduino.h"

//#define PRINT 1
extern double scale[];
extern int s_length;
extern volatile long base_period;
extern char Tuning;
extern  boolean SingleVoice;
extern boolean Strata;
extern int Layers[4];
extern int HarmonicCount;
extern const boolean Debug;
extern boolean RingChanges;
extern boolean Hocket;
extern double Root;
extern boolean ChangeRoot;
extern const double pi;
extern boolean Glide;
extern boolean Sync;
extern int NoteDuration;
extern class RSS * Voices[];
extern boolean ShiftRegisterMode;
extern const boolean Quantize;
extern const boolean Phasing;
extern boolean Detune;
extern boolean Modulate;
extern boolean Flicker;
extern boolean BZZZ;
extern boolean ProximityMode;
extern boolean FillIn;
extern const int TrigPins[];
extern const boolean InputMode;
extern const int LeftIn;
extern volatile boolean knobswitch;
extern void usDelay(unsigned long);
extern void Send2DAC(int, int);
extern unsigned long ModDelay(double);
extern void Swap(double*, int, int);
extern double Stack[];
extern double TrigOrder[];
extern double LFO(double, double, double, boolean);
extern void DepleteNoteDuration();
extern double ReadLeft();

class RSS {
public:
	RSS(int N, int , char, int);
	void GenerateTriggerSequence(int);
	void bjorklund(int,int);
	void ShiftBoolean(int);
	void PrintTriggerSequence();
	void Advance();
	double FixInterval(double);
	double DifferenceInterval(double, int);
	double InharmonicInterval();
	double HarmonicInterval();
	double JustInterval(int, int);
	//void Swap();
	double Saw(double, double, double, bool);
	//double LFO(double, double, double, bool);
	double GetCV(int);
	bool GetTrig(int);
	void FlipT();
	void Invert();
	void Reverse();
	void OrderSequence();
	void Transpose(double);
	void Shift(int, int);
	void SetNoteDuration(int);
	double ArrayMin(double[]);
	double ArrayMax(double[]);
	void sendNote(int);
	void sendTrigger();
	double QuantizeNote(double);
	//unsigned long ModDelay(boolean, double);
	int flicker(int);
	//void usDelay(unsigned long);
	// void DepleteNoteDuration();
	// int GetLength();
	bool Play(int);
	double ApplyShiftRegister();
	void EvaluateProximity();
	void RandomizeInitialOrder();
	int SubD;
	bool G;
	bool B;
	bool _D;
	int _N;
	int _NoteDuration;
	int TriggerSubdivisions;


private:
	bool *_TrigArray;
	double *_CV;
	char _Rule;
	double _roof;
	double _floor;
	int _pos;
	int _prevCV;
	bool _Prox;
	int _HarmonicCount;
	double _Interval;
	int _n; // scale index
	bool _Dir; // direction of CV movement
	int _channel;
	int _iteration;
	bool _Play;
	int _counter;
	bool RandomStep(double);
	void UpdateInterval();
};

//#endif