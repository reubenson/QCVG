#include "Arduino.h"
#include "QCVG_VOICE.h"
#include "Scale.h"


QCVG_VOICE::QCVG_VOICE(int N , int hits, char Rule, int channel){
/*
  Constructor for a 'voice', which is defined by a sequence of notes.
  QCVG::Advance will build the sequence of notes. With each advance, the next note in
  the sequence will be determined according to the specified Rule.

  In an analogue modular synthesizer, the pitch of an oscillator is voltage-controlled, while 
  the base frequency may be freely changed by hand. Notes are therefore determined relative to the 
  oscillator's base frequency. In the case of Eurorack synthesizers, the voltage-control of 
  an oscillator's frequency is defined by the volt-per-octave standard. For example, an input 
  voltage of 0V produces no change from the base frequency, and an input voltage of 2V produces 
  a two-octave increase in the oscillator's frequency. The DAC chips in the QCVG have a range of
  0-4 volts (4 octaves of pitch control), with mV precision. Here, notes will be defined as
  intervallic ratios that range from 1:1 (base frequency) to 16:1 (4 octaves above the base frequency).

  Along with pitch, each voice has a sequence of triggers determined by the Euclidean Sequencer method, 
  in which the number of triggers in a sequence are distributed to maximize spacing between triggers. 
  These triggers are used in conjunction with envelope generators in the synthesizer, which provide the 
  amplitude contour of the note. The oscillator's output is multiplied with the envelope in a VCA, or
  voltage-controlled amplifier, to produce an articlated note.

  The sequence will repeat for the number of times defined by Repetitions. This implementation
  allows general functionality that allows both a 'free-running' mode of note generation, as well as 
  more traditional sequencing modes if Repetitions>1.
*/

  // Set parameters according to arduino sketch
  length = N;                 // length of sequence
  _TrigArray = new bool[N];   // trigger array
  _NoteArray = new double[N]; // note array
  _Rule = Rule;               // rule for generating notes
  _channel = channel;         // unique ID (0,1,2,3)

  // Default Parameters
  TriggerSubdivisions = 1;    // number of trigger subdivisions
  _Play = true;               // enable voice on
  TimeFactor = 1.0;           // multiply clock speed by this factor
  _n = 0;                     // index within scale
  counter = 0;                // number of times note has been advanced in sequence
  iteration = -1;             // number of times sequence has been updated
  clock = 0;                  // internal clock
  Glide = Glissando;          // set note glide

  if (RingChanges) { _n = Scale::length; }  // start from highest pitch

  // _floor and _roof define the bounds of each voice's pitch
  if (Strata) { 
    // each voice sits in its own octave register
    _floor = pow(2.,Layers[_channel]);
    _roof  = pow(2.,Layers[_channel]+1);
  }
  else {
    // each voice spans entire four octaves
    _floor = pow(2.,0); 
    _roof = pow(2.,4); 
  }

  if (SingleVoice) { SetupSingleVoice(); } // determine relative intervals

  // generate 'Euclidean' trigger sequence
  GenerateTriggerSequence(hits);
 
}

void QCVG_VOICE::GenerateTriggerSequence(int hits){
// Euclidean sequence generator:
// Evenly distribute the number of notes over a sequence of given length
// (only determines timing/triggers, not pitch)

  // hits: number of notes in a sequence of the defined length
  hits = constrain(hits,0,length);

  // initialize the sequence with 1s and then 0s
  for (int i = 0 ; i<length ; i++) {                      
    if (i<hits) { _TrigArray[i] = true; }
    else        { _TrigArray[i] = false; } 
  }

  // obtain Euclidean Sequence
  Bjorklund(length,hits);
}

void QCVG_VOICE::Bjorklund(int length , int h) {
// Implementation of the Bjorklund algo for even distribution of timing events
// For reference, see http://cgm.cs.mcgill.ca/~godfried/publications/banff.pdf 

  int j = 2;
  int iterations = -1+length / min(h,length-h);
  int leftovers =  length % min(h,length-h);

  for (int k = 0; k<iterations; k++) {
    int l = min(h,length-h);
    for (int i = 0; i<l; i++ ) {
      ShiftBoolean( i*j); 
      if (l<h) {
        h--;
      }
      length--;
    }  
    j++;
  }

  int l = min(h,length-h);
  for (int i = 0; i<l; i++ ) {
    ShiftBoolean( i*j*h/l-i*j/l);
    length--;  
  }
}

void QCVG_VOICE::ShiftBoolean( int a ) {
// Subroutine used in conjunction with Bjorklund
// Within trigger array, move the boolean at the end of the array to position after input 'a'
   
   // move values one position down to make room for the transported value
  bool temp = _TrigArray[length-1];
   for (int i = length-1 ; i > a+1 ; i--) {
     _TrigArray[i] = _TrigArray[i-1]; 
   }
   _TrigArray[a+1] = temp;
}


void QCVG_VOICE::PrintTriggerSequence() {
// Print trigger sequence to serial monitor

  Serial.print("Trigger Sequence "); Serial.print(_channel); Serial.println(": ");
  for (int j = 0; j < length ; j++) {
    Serial.print(_TrigArray[j]); Serial.print("  "); 
  }
  Serial.println();  
}

void QCVG_VOICE::PrintCVSequence() {
// Print CV Sequence to serial monitor

  Serial.print("\nVoice"); Serial.print(_channel+1); 
  Serial.print(", Iteration ");
  Serial.print(iteration); Serial.print(": \n");

  for (int i = 0; i < length ; i++) {
    if (_TrigArray[i]) { Serial.print("_____  "); }
    else               { Serial.print("       "); }
  }
  Serial.println();
  for (int i = 0; i < length ; i++) {
    Serial.print(_NoteArray[i]); Serial.print("  "); 
    if (_NoteArray[i]<10) {Serial.print(' ');}
  }
  Serial.println();
}

void QCVG_VOICE::Advance() {
// Determine the next note in the sequence as determined by program mode and rule
  
  // number of times a Voice has been advanced
  int _counter = counter-1; 

  // return early if condition not met
  if (_counter%(Repetitions*length) >= length) { 
    return; 
  }

  double PrevCV = Voices[PreviousVoice()]->GetCV();
  static double ShiftRegister[4] = {1.0,1.0,1.0,1.0};
  double new_note = _NoteArray[CurrentPosition];

  static bool PopulateSequenceOnce = RingChanges;
  static bool SingleVoiceModes = SingleVoice || ShiftRegisterMode || ET;

  if ( (((!SingleVoiceModes) || _channel == 0) && !PopulateSequenceOnce) || (PopulateSequenceOnce && _counter<length)){
    
    // determine new note according to _Rule
    switch (_Rule) {

      case 'A': { 
        static boolean direction = _channel%2;
        double increment = pow(Scale::s[random(1,Scale::length)],pow(-1,direction));
        if (new_note * increment > _roof || new_note * increment < _floor) { 
          direction = !direction; increment = 1./increment;
        }  
        new_note *= increment; 
        break;
      }
                
      case 'B': { 
        if (new_note* PrevCV > _roof || new_note* PrevCV < _floor) {
          PrevCV = 1./PrevCV; FixInterval(PrevCV);
          _roof=16;
        }
        new_note *= PrevCV; 
        break; 
      }
    
      //  LFO mode
      case 'C': {
        static double freq = 0.005;
        double freq_offset = _channel*0.01;
        double sine = 2+LFO(2,freq+freq_offset,0,1,0);
        new_note = pow(2,sine);
        break; 
      }
       
      // generate intervals that produce specific difference tones         
      case 'D': { 
        double modulation = LFO(1,0.01,0,0,0);
        double span = 1.;
        if (InputMode) { span = ReadInput(1); }
        new_note = Scale::DifferenceInterval(_channel,modulation,span); 
        break; 
      }
    
      // generate random note
      case 'E': {
        double low = 1000*log(_floor)/log(2); 
        double high = 1000*log(_roof)/log(2);
        double rand = double(random(low,high++))/1000.;
        new_note = pow(2,rand);
        break;
      }   

      // random walk through scale
      case 'G': { 
        int n_high = (_roof/_floor)*Scale::length/Scale::s[length-1];
        new_note = _floor*Scale::Get(_n);
        // cv = constrain(cv,_floor,_roof);
        _n += random(-1,2);
        _n = constrain(_n,0,n_high);  
        break; 
      }

      // generate harmonic intervals
      case 'H': { 
        double base = 1.;
        if (InputMode) { base = ReadInput(1); }
        new_note = Scale::HarmonicInterval(1,16,base); 
        break; 
      }

      // generate just intonation intervals
      case 'J': {
        static int limit = 7;  // n-limit of just intervals
        static int power = 1;  // just interval power
        // cv = ShiftRegister[_channel]*Scale1.JustInterval(limit,power);
        new_note = PrevCV*Scale::JustInterval(limit,power);
        break;
      }
                
      // derive note from instantaneous value of sawtooth LFO                                                        
      case 'R': {
        double freq = 0.5/(QCVG.Clock()/1.e6);    // define relative to Clock
        if (InputMode) { freq *= ReadInput(1); }  // allow user input to scale

        double Amplitude = 0.5*log(_roof/_floor)/log(2);
        double Offset = _floor*Amplitude;
        new_note = pow(2.,Offset+LFO(Amplitude, freq, _channel/4.0, 1,_channel%2));
        break;
      }
           
      // ascend and descend scale                   
      case 'S': {
        new_note = _floor*Scale::Get(_n);

        // apply reflecting boundary condition
        if (new_note<=_floor) {direction = 1; new_note = _floor;}
        if (new_note-_roof>-1e-3) {direction = 0; new_note = _roof;}

        int stepsize = 1;
        if (RandomStep(0.3)) {stepsize+=_channel;}
        if (direction) {_n+=stepsize;}
        else           {_n-=stepsize;}
        break; 
      }
    }

    // quantize to nearest note in Scale
    if (Quantize) {  Scale::QuantizeNote(new_note);  }

  }

  else if (RingChanges && _counter%length==0) {
    Swap(_NoteArray, iteration, length);
    new_note = _NoteArray[CurrentPosition];
  }

  else if (ET) {
    switch (_channel) {
      case 1: new_note = Voices[0]->GetCV()*2.; break;
      case 2: new_note = Voices[0]->GetCV(length-1-CurrentPosition); break;
      case 3: new_note = Voices[0]->GetCV(length-1-CurrentPosition)*2; break;
    }
  }

  else if (SingleVoice){ 
    // Voices 2-4 pitch fixed with respect to by Voice1

    // allow user input of FM Carrier:Modulator Ratio on the fly via Serial
    if (Serial.available() && _channel==3) {_Interval = Serial.parseFloat();}

    new_note = Voices[0]->GetCV()*_Interval;

  }
  
  else if (ShiftRegisterMode) {
    new_note = ShiftRegister[_channel];
  }

  // apply Root Shift (Root = 1.0 if !ChangeRootMode)
  new_note *= QCVG.Root;
  
  // make sure new_note is between 1 and 16
  FixInterval(new_note);

  // update sequence with new note
  _NoteArray[CurrentPosition] = new_note;

  // enforce temporal cohesion in SingleVoice mode
  if (_channel!=0 && SingleVoice) { TimeFactor = Voices[0]->TimeFactor; }
  // TimeFactor = 1./pow(new_note,0.5);

  // update shift register 
  if (_channel == 0){
    for (int i=3;i>0;i--) {
      ShiftRegister[i] = ShiftRegister[i-1];
    }
    ShiftRegister[0] = new_note;
  }

}

void QCVG_VOICE::FixInterval( double &Interval) {
// constrain Interval to periodic boundaries defined by _roof and _floor
  while (Interval-_roof>1e-3) {Interval /= _roof/_floor;}
  while (_floor-Interval>1e-3) {Interval *= _roof/_floor;}
}

double QCVG_VOICE::GetCV( int i ) {
// return CV corresponding to CurrentPosition

  if (i==-1) {
    // default (no arg)
    return _NoteArray[CurrentPosition];
  }
  else {
    return _NoteArray[i%length];
  }
}

bool QCVG_VOICE::GetTrig( long time, bool &AdvanceMe) {
// Determine whether to trigger the next note

  // offset timing for different voices
  long _Spacing = TimingOffset*QCVG.Clock()*QCVG.TrigOrder[_channel];

  // in Phasing mode, modulate the amount of timing offset with quadrature lfo
  if (Phasing) {
    _Spacing = QCVG.Clock()/2+LFO(QCVG.Clock()/2,0.01,_channel/4.0,0,0);
  }

  time -= _Spacing;

  unsigned long bp = QCVG.Clock()*TimeFactor;

  // Determine duration of note and sequence
  // default: all sequences will have the same time duration,
  // regardless of length (the number of notes in the sequence),
  // i.e. longer sequences will have shorter notes
  long SingleNote= bp*(Voices[0]->length)/length;
  //SingleNote = bp; // Give each note the same duration

  // Length of glissando determined as fraction of SingleNote
  _GlideLength = Glide*SingleNote;

  bool AlwaysAdvance = true;
  // true: each note trigger will correspond to moving to the next note in the sequence
  // false: advance to next note in sequence only after equal note durations (ignores TimeFactor)
  //   allows TimeFactor to work as a trigger multiplier in which the same note is repeated

  AdvanceMe = false;

  bool send_trig = false;
  if (time-clock >= SingleNote) {
    if (_Play) {
      send_trig = _TrigArray[counter%length];

      // add exception to allow RingChanges mode to work for hits<length
      if (RingChanges && counter<length) {send_trig = true;}

      // this definition of clock prevents accumulation of timing error
      clock += SingleNote;

      if (send_trig) {
        TriggerStart = clock+_Spacing;
        if (AlwaysAdvance) {
          CurrentPosition = counter%length;
          AdvanceMe = true;
        }
        else {
          if (time-advanceclock >= SingleNote/TimeFactor) {
            advanceclock+=SingleNote/TimeFactor;
            CurrentPosition = counter%length;
            AdvanceMe = true;
          }
        }
      }
      counter++;
      if (counter%(length*Repetitions)==0 && counter>0) {
        iteration++;
        
        // print sequence to serial monitor at the end of sequence
        if (Debug) {
          if (counter%(Repetitions*length) == 0) {
            PrintCVSequence();
          }
        }
      }
    }
  }

  return send_trig;
}



double QCVG_VOICE::ArrayMax(double * array) {
// Return maximum value in input array
  double m = array[0];
  for (int i = 1; i < 4; i++) {  m = max(m,array[i]);  }
  return m;
}

double QCVG_VOICE::ArrayMin(double * array) {
// Return minimum value in input array
  double m = array[0];
  for (int i = 1; i < 4; i++) {  m = min(m,array[i]);  }
  return m;
}

void QCVG_VOICE::sendTrigger() {    
// Send trigger ON to TrigPin
// (trigger OFF is handled by QCVG_MAIN::DepleteNoteDuration)

  if ( (counter-1)%TriggerSubdivisions!=0) {return;}

  digitalWrite(TrigPins[_channel], HIGH);
}

bool QCVG_VOICE::Play(int t){
// Return true if voice is enabled
  return _Play ? true : false;
}

void QCVG_VOICE::EvaluateProximity(){
// Used for ProximityMode, in which the detected proximity of the current
// voice's pitch to that of at least one other voice is used to send trigger
  double ProximityRange = 0.05;
  bool WithinProximity = false;
  for (int i = 0; i < 4; i++) {
    if (i != _channel) {
      if ( abs( _NoteArray[CurrentPosition] - Voices[i]->GetCV() ) < ProximityRange ) {
        digitalWrite(TrigPins[i],HIGH); 
      }
      else {
        digitalWrite(TrigPins[i],LOW); 
      }
    }
  }
}

void QCVG_VOICE::RandomizeInitialOrder() {
  int a, b;
  for (int i=0; i<length/2; i++) {
    a = random(0,length);
    b = random(0,length);
    while (a==b) { b = random(0,length); }
    double temp = _NoteArray[a];
    _NoteArray[a] = _NoteArray[b];
    _NoteArray[b] = temp;
  }
}

bool QCVG_VOICE::RandomStep(double probability) {
  int r = random(0,1001);
  if (r<1000*probability) {
    return true;
  }
  else { 
    return false;
  }
}

void QCVG_VOICE::SetupSingleVoice() {
// Define 4-note chord for use in SingleVoice mode

  double chord[4];

  switch (Voicing) {
    case 1: {chord[0] = 1; chord[1] = 2; chord[2] = 3; chord[3] = 4; break;}
    case 2: {chord[0] = 1; chord[1] = 2; chord[2] = 4; chord[3] = 6; break;}
    case 3: {chord[0] = 1; chord[1] = 3; chord[2] = 5; chord[3] = 7; break;}
    case 4: {chord[0] = 2; chord[1] = 3; chord[2] = 5; chord[3] = 6; break;}

    case 5: {chord[0] = 63; chord[1] = 62; chord[2] = 56; chord[3] = 42; break;}

    case 6: {chord[0] = 1; chord[1] = 2300./270; chord[2] = 3000./270; chord[3] = chord[2]; break;}
    case 7: {chord[0] = 1; chord[1] = 870./300; chord[2] = 2250./300; chord[3] = chord[0]; break;}
    case 8: {chord[0] = 1; chord[1] = 2000./400; chord[2] = 2250./400; chord[3] = chord[0]; break;}
    case 9: {chord[0] = 1; chord[1] = 1700./660; chord[2] = 2400./600; chord[3] = chord[0]; break;}
    case 10: {chord[0] = 1; chord[1] = 1700./660; chord[2] = 2400./600; chord[3] = chord[0]; break;}

    case 11: {chord[0] = 3./2; chord[1] = 5./3; chord[2] = 8./5; chord[3] = 13./8; break;}
    case 12: {chord[0] = 6./5; chord[1] = 8./7; chord[2] = 12./11; chord[3] = 14./13; break;}
    case 13: {chord[0] = 8; chord[1] = 4; chord[2] = 2; chord[3] = 1; break;}

  }

  // update the intervals by which to multiply Voices2-4
  for (int i = 0; i<4; i++){
    Voices[i]->_Interval = chord[i] / chord[0];
  }

  // restrict pitch of Voice1 such the highest interval of the four-voice chord doesn't exceed 16
  if (SingleVoice) {
    Voices[0]->_roof = 16.*(chord[0]/ArrayMax(chord));
    Voices[0]->_floor = (chord[0]/ArrayMin(chord));
  }
}

int QCVG_VOICE::PreviousVoice() {
// Return the previous voice in TrigOrder

  // Determine position of current voice in TrigOrder
  int i = 0;
  while (_channel!=QCVG.TrigOrder[i]) {
    i++;
  }

  return QCVG.TrigOrder[(i+3)%4];
}

void QCVG_VOICE::UpdateNote() {
// Update pitch information (CV) according to current note and modulation+glissando

  double note = _NoteArray[CurrentPosition];

  // ProximityMode is an alternate mode under development
  if (ProximityMode) {EvaluateProximity();}

  double CV = 1000.*log(note)/log(2);   // frequency ratio to CV (millivolts)

  // apply LFO modulation
  if (Modulate)  {
    double Amp = 500; // (millivolts)
    if (InputMode) {
      Amp *= ReadInput(1);
    }

    CV += log(LFO(Amp,50,_channel/4.0,0,0))/log(2);
  }
  
  // if (CV_prev == CV) {return;} // No need to update CV if previous value == new value

  // apply portamento/glissando
  if (_GlideLength>0) {
    // Time since note start
    unsigned long time_note = micros()-TriggerStart;

    // interpolate between previous CV and target CV
    if (time_note<=_GlideLength) {
      double fraction = double(time_note)/double(_GlideLength);
      CV = CV_prev*(1-fraction)+fraction*CV;
    }
  }

  // send note to DAC 
  Send2DAC(_channel,CV);

  // reference for next updateNote()
  CV_prev = CV;

  // Serial.print("Voice "); Serial.print(_channel); Serial.print(": "); Serial.print(note); Serial.println();
}

void Swap( double *Sequence , int iteration , int NumberOfBells) {
// Swap members of Sequence according to the 'plain bob' method, which is the simplest instance of method ringing.
// In method ringing, the goal is to play each possible permutation of a sequence of notes rung on bells.
// For reference, see http://plus.maths.org/content/ringing-changes
  int pos;
  if (iteration%2 == 0) {pos = 0;}
  else {pos = 1;}
  
  // after NumberOfBells*2, an extra increment on pos is imposed
  if ( (iteration+1)%NumberOfBells*2 == 0 ) {
    pos++;
  }
  
  // starting with pos, swap adjacent members of sequence
  while ( pos < (NumberOfBells-1) ) {
    double storage = Sequence[pos];
    Sequence[pos] = Sequence[pos+1];
    Sequence[pos+1] = storage;
    pos = pos+2;
  }
}

void QCVG_VOICE::Transpose(double TransposeAmount) {
  for (int i=0; i<length; i++) { 
    _NoteArray[i] *= TransposeAmount;
    FixInterval(_NoteArray[i]);
  } 
}

void QCVG_VOICE::OrderSequence() {
// Put CV array in order of ascending pitch
  bool recurse = false; 
  double temp;
  for (int i=1; i<length; i++) { 
    if (_NoteArray[i]<_NoteArray[i-1]) {
      temp = _NoteArray[i]; 
      _NoteArray[i] = _NoteArray[i-1]; 
      _NoteArray[i-1] = temp; 
      recurse = true; 
    }
  }
  if (recurse) { OrderSequence(); } 
}

void QCVG_VOICE::Reverse() {
// Reverse the order of CV array
  double temp;
  for (int i=0; i<(length-1)/2; i++) { 
    temp = _NoteArray[i]; 
    _NoteArray[i] = _NoteArray[length-(i+1)]; 
    _NoteArray[length-(i+1)] = temp; 
  } 
}

void QCVG_VOICE::Invert() {
// Invert CV array with respect to _roof
  for (int i=0; i<length; i++) {
    _NoteArray[i] = _roof/_NoteArray[i]; 
  } 
}

void QCVG_VOICE::FlipT() {
// Flip booleans in trigger array
  for (int i=0; i<length; i++) { 
    _TrigArray[i] = !_TrigArray[i];
  } 
}

void QCVG_VOICE::Shift( int Tshift , int Nshift ) {
// Rotate CV array and trigger array by input amount
  bool TempTrigArray[length]; 
  double TempNoteArray[length];
  for (int i=0; i<length; i++){
    TempTrigArray[i] = _TrigArray[i];
    TempNoteArray[i] = _NoteArray[i];
  }

  for (int i = 0; i < length;  i++) {
    _TrigArray[i] = TempTrigArray[(i+Tshift+length) % length];
    _NoteArray[i]        = TempNoteArray[(i+Nshift+length) % length]; 
  }
}
