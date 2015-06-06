#ifndef QCVG_MAIN_h
#define QCVG_MAIN_h

#include "Arduino.h"
#include <SPI.h>

// ---------------------- Arduino Nano Pin Definitions -------------------------
const int CS1 = 8;						  // DAC Chip Select 1 (CV Output 1 & 2) 
const int CS2 = 9;               		  // DAC Chip Select 2 (CV Output 3 & 4)
const int TrigPins[] = {4,5,6,7};         // Trigger Output Pins
const int LeftIn = 21, RightIn = 20;      // Input Left (A7) & Right (A6)
const int XPin = 3, YPin = 12, ZPin = 2;  // Rotary Encoder Pins

extern class QCVG_VOICE * Voices[];
extern volatile boolean SwitchedOn;
extern boolean Debug;
extern boolean SingleVoice;
extern boolean ET;
extern char Rule;
extern int TriggerLength;
extern boolean ChangeRootMode;
extern double TrigOrder[];
extern double root_list[4];

const double pi = 3.14159;

void DepleteNoteDuration();
void UpdateTiming();
void PressRelease();
void KnobSpin();
void Send2DAC(const int, const int);
double ReadInput(bool);
double LFO(double,double,double,int,bool);

class QCVG_MAIN {
public:
	QCVG_MAIN(unsigned long);
	bool GetTrigs(bool*);
	void loop();
	void TickOver();
	double Root;
	double TrigOrder[4];
	unsigned long Clock();
	unsigned long clock_period;
	void ModulateClock(double,double);

private:
	unsigned long t;
	unsigned long StartTime;
	double a1;
	double f1;
};

#endif