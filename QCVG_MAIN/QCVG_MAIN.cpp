#include "QCVG_MAIN.h"
#include "Arduino.h"
#include "QCVG_VOICE.h"

/*
  QCVG_MAIN consists of functions that handle the main loop of the program,
  communicate with the two MC4822 DACs via SPI, and handle interrupts from
  user input via the QCVG's rotary encoder and momentary switch.
*/

QCVG_MAIN::QCVG_MAIN(unsigned long input) {

  SPI.begin();

  clock_period = input;

  pinMode(CS1,OUTPUT);     digitalWrite(CS1, HIGH);       
  pinMode(CS2,OUTPUT);     digitalWrite(CS2, HIGH);
  for (int i=0; i<4; i++) {  
    pinMode(TrigPins[i],OUTPUT);   
    digitalWrite(TrigPins[i],LOW);
  }
  if (InputMode)  {
    pinMode(LeftIn,INPUT); 
    digitalWrite(RightIn,INPUT); 
  }
  pinMode(XPin,INPUT);     digitalWrite(XPin, HIGH);      
  pinMode(YPin,INPUT);     digitalWrite(YPin, HIGH);
  pinMode(ZPin,INPUT);     digitalWrite(ZPin, HIGH);    

  attachInterrupt(1,KnobSpin,CHANGE);
  attachInterrupt(0,PressRelease,CHANGE); 

  randomSeed( analogRead(A3) );

  Root = 1.0; 

  // initialize TrigOrder (order in which Voices are played)
  for (int i=0; i<4; i++) {
    TrigOrder[i] = i;
  }

  // initialize amplitude and frequency of clock modulation
  a1 = 0.; f1 = 1.; 

  StartTime = micros();
}

void QCVG_MAIN::loop() {
  if (SwitchedOn) {

    bool TriggerList[4] = {false,false,false,false};

    // determine which (if any) triggers need to be sent
    bool increment = GetTrigs(TriggerList);

    // for each trigger, determine the accompanying CV value
    for ( int l=0; l<4; l++) {
      int v = TrigOrder[l];
      if (TriggerList[v] || ProximityMode) {
        Voices[v]->sendTrigger();
      }
    }

    // determine if any triggers need to be set to LOW
    DepleteNoteDuration();

    // update CV
    for (int i=0; i<4; i++) {
      Voices[i]->UpdateNote();
    }

    // end of tick routines
    if (increment) {
      TickOver();
    }
  }
}

bool QCVG_MAIN::GetTrigs(bool *TriggerList) {
  if (ProximityMode) {return false;}

// determine which voices need to be triggered
  bool Trigger = false;
  unsigned long time = micros();

  for (int i=0; i<4; i++) {
    int j = TrigOrder[i];
    bool AdvanceMe;
    if ( Voices[j]->GetTrig(time,AdvanceMe) ) {
      Trigger = true;
      TriggerList[j] = true;

      // Advance() takes the most time, so it's processed before trigger and CV are updated
      if (AdvanceMe) {
        Voices[j]->Advance();
      }
    }
  }

  // if there are any voices to trigger, return true
  return Trigger;
}


void QCVG_MAIN::TickOver() {
  t++; 

  if (BZZZ) {
    clock_period *= 0.975;
    if (clock_period < 50) {
      clock_period = random(50,200); 
    }
  }

  if ( t%(Voices[0]->length*Repetitions) == 0) {
    // Repetitions = random(20,25);
  }
  
  if (ChangeRootMode) {
    if ( t%(Voices[0]->length*Repetitions) == 0) {
      static int root_i;
      static int root_list_length = sizeof(root_list)/sizeof(double);
      root_i = (root_i+1) % root_list_length; 
      Root = root_list[root_i];
    }
  }
  
  if (ET && t%(Voices[0]->length*Repetitions) == 0) {
    for ( int l=0; l<4; l++) { 
      // Voices[l]->Reverse(); 
    }  
  }

  if (ChangeTrigOrder && t%4==0) { 
    Swap(TrigOrder,t/4-1,4); 
  }
    
}

void DepleteNoteDuration() {
// Turn trigger off, based on TriggerLength
// i.e. TriggerLength==10 yields a 10ms ON/OFF
// non-instantaneous triggers are useful for more complex envelope generation
  unsigned long time = micros();
  for (int i=0; i<4; i++) { 
    if (time-Voices[i]->TriggerStart >= 1000.*TriggerLength) {
      digitalWrite(TrigPins[i],LOW); 
    }
  }
}

void UpdateTiming() {
  unsigned long now = micros();
  for (int i=0; i<4; i++) {
    Voices[i]->TriggerStart = now;
    Voices[i]->clock = now;
  }

  // Call loop to stabilize?
  loop();
}

void PressRelease() {
// Toggle start/stop by pressing rotary encoder button

  static unsigned long last_press, t_press;

  bool state = digitalRead(ZPin); 
  if (millis() - last_press > 10) {  
    last_press = millis();
    if (!state) { 
      t_press = millis(); //??
    }
    else { 
      if (millis()-t_press < 1000) {
        SwitchedOn = !SwitchedOn;
      }
    //else {PO++; Rule = PatchOrder[PO % 2]; t = -1; }
    }
  UpdateTiming();
  } 
}

void KnobSpin() {
// Speed up or slow down clock speed based on rotation of rotary encoder

  static unsigned long last_rotation;

  // ignore if interrupts come faster than 600ms
  if (millis() - last_rotation > 600) {  
    if (Debug) {Serial.println("KnobSpin");}
    last_rotation = millis();   boolean clockwise;
    if (digitalRead(XPin) == HIGH) {
      if (digitalRead(YPin) == LOW) { clockwise = false; }
      else                          { clockwise = true; } 
    }
    else {                             
      if (digitalRead(YPin) == LOW) { clockwise = false; }
      else                          { clockwise = true; } 
    }
    if (clockwise) { QCVG.clock_period *= 0.75; }
    else           { QCVG.clock_period *= 1.25; }
    UpdateTiming();
  } 
}


void Send2DAC(const int channel, const int val) {
// Convert pitch interval to CV via DAC

  // The 12-bit MCP4822 chip allows a maximum of 4096 mV
  int value = constrain(val,0,4096);

  // send value to DAC channel via SPI as 2 byte word
  if (channel<2) {
    digitalWrite(CS1, LOW);
    if (channel == 0)  { SPI.transfer(B00010000 | (value >> 8) ); }
    else               { SPI.transfer(B10010000 | (value >> 8) ); }
    SPI.transfer(value | B00000000);
    digitalWrite(CS1, HIGH); }
  else { // channel == 2 || channel == 3
    digitalWrite(CS2, LOW);
    if (channel == 2)  { SPI.transfer(B00010000 | (value >> 8) ); }
    else               { SPI.transfer(B10010000 | (value >> 8) ); }
    SPI.transfer(value | B00000000);
    digitalWrite(CS2, HIGH); }
}

double ReadInput(bool Left) {
// Read signal at left or right input of QCVG
// return value normalized to 1
  
  // average over five readings of input value (1-1024)
  int Pin;
  if (Left) { Pin = LeftIn; }
  else      { Pin = RightIn; }

  double value = 0;
  for (int i=0; i<5; i++) {
    value += analogRead(Pin)+1;
  }
  value /= 5.;

  // return normalized value
  return value/1024.;
}

double LFO( double amplitude , double frequency , double phase = 0.0 , int waveform = 0, bool invert = 0 ) {
  // Low Frequency Oscillator - Modulation Source

  unsigned long period_us = 1.0e6/frequency;           // period of LFO (us)
  unsigned long time_us = micros()+phase*period_us;    // current time (us)
  double time = double(time_us%period_us)/1.0e6;       // current time (s) relative to period

  // Calculate position within waveform
  double wave;
  switch (waveform) {
    // Sine
    case 0: { wave = sin(frequency*time*2.0*pi+phase*2.0*pi); break;}

    // Saw
    case 1: { wave = 2.0*double(time)*double(frequency) - 1; break; }

    // Square
    case 2: { wave = bool(time>(period_us/2.0e6)); break;}
  }

  if (invert) {
    wave *= -1;
  }

  return amplitude*wave;
}

unsigned long QCVG_MAIN::Clock() {
// returns current value of clock (in microseconds)

  // add layer of modulation
  clock_period *= 1+LFO(a1,f1);

  return clock_period*1000;
}

void QCVG_MAIN::ModulateClock(double amp, double freq) {
  a1 = amp;
  f1 = freq;
}
