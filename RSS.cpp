#include "Arduino.h"
#include "RSS.h"
// #include "MemoryFree.h"
//#include "SPI.h"

//using namespace grubs;

// using namespace grubs;

// namespace grubs{
//   bool c = 1;
// }
//using namespace std;

// namespace grubs{
//   int blue;
// };

RSS::RSS(int N , int hits, char Rule, int channel){
  _TrigArray = new bool[N];
  _CV = new double[N];
  _Rule = Rule;
  _pos = 0;
  SubD=1;
  _Prox = 0;

  if (Strata) {
    _floor = pow(2.,channel) ; _roof = pow(2.,channel+1); }
  else {
    _floor = 1.0; _roof = 16.0; }

  _D = !Sync; // Delay after note
  G = Glide; // Glide On
  _NoteDuration = NoteDuration;
  _N = N;
  B = false;
  _channel = channel;
  _iteration = 0;
  _Play = true;
  _counter = 0;
  _HarmonicCount = (int) _floor;
  _Dir = (channel+1)%1;
  TriggerSubdivisions = 1;

  // Write _CV to initial state (i.e. [1,2,3,4] for N = 4)
  // This initial state will be overwritten if !RingChanges
  for (int i = 0; i < N; i++) {
    _CV[i] = i+1.;
  }

  // UpdateInterval();

  RandomizeInitialOrder();

  GenerateTriggerSequence(hits);
  if (Hocket) {
    Shift(-1*(channel-1),0);
  }


}

void RSS::GenerateTriggerSequence(int hits){
 // hits: number of notes in a bar of the defined length
  hits = constrain(hits,0,_N);

  // initialize the sequence
  for (int i = 0 ; i<_N ; i++) {                      
    if (i<hits) { _TrigArray[i] = true; }
    else        { _TrigArray[i] = false; } }

  // obtain Euclidean Sequence
  bjorklund(_N, hits);   
  //if (Debug) { PrintTriggerSequence(seq); }  
}

void RSS::bjorklund(int length , int h) {
  int j = 2;
  int iterations = -1+length / min(h,length-h);
  int leftovers =  length % min(h,length-h);

  for (int k = 0; k<iterations; k++) {
    int l = min(h,length-h);
    for (int i = 0; i<l; i++ ) {
      ShiftBoolean( i*j); 
      if (l<h) {h--;}
      length--;  }  
  j++;
   }

    int l = min(h,length-h);
    for (int i = 0; i<l; i++ ) {
      ShiftBoolean( i*j*h/l-i*j/l);
      length--;  }
}

void RSS::ShiftBoolean( int a ) {
   // within the boolean array, move the boolean at the end of the array to position after 'a'
   
   bool temp = _TrigArray[_N-1];
   
   // move values one position down to make room for the transported value
   for (int i = _N-1 ; i > a+1 ; i--) {
     _TrigArray[i] = _TrigArray[i-1]; }
     
   _TrigArray[a+1] = temp;
     //if (Debug) { PrintTriggerSequence(seq);
}

void RSS::PrintTriggerSequence() {
  for (int i = 0; i < _N ; i++) {
    Serial.print(_TrigArray[i]); Serial.print("  "); }
  Serial.println();  
}

void RSS::Advance() {

  double Interval2 = Voices[(_channel+3)%4]->GetCV(0), cv = _CV[_pos];

  if ( (((!SingleVoice && !ShiftRegisterMode) || _channel == 0) && !RingChanges) || (RingChanges && _counter<_N)){
    switch (_Rule) {
      case 'A': { double increment = pow(scale[random(1,s_length)],pow(-1,_Dir));
              if (cv * increment > _roof || cv * increment < _floor) { 
                _Dir = !_Dir; increment = 1./increment;
              }  
                cv *= increment; 
                break;}
                
      case 'B': { if (cv* Interval2 > _roof || cv* Interval2 < _floor) {
                    _roof = 2; Interval2 = 1/FixInterval(Interval2); _roof=16;} 
                  cv *= Interval2; 
                  break; }
    
      //  LFO mode
      case 'C': {double sine = 2+LFO(2,0.005+_channel*0.01,0,1);
                cv = pow(2,sine);
                break; }
                
      case 'D': { cv = DifferenceInterval(Interval2,_channel); break; }
    
      case 'E': {
        switch (Tuning) { 
          case 'J': { // Just Intonation
            int JI_limit = 5;  // n-limit of just intervals
            int JI_power = 1;  // just interval power
            cv = JustInterval(JI_limit,JI_power);  
            break;}      
          case 'E': { // Scale with defined number of equidistant notes per octave
            double n_scale = 12.;  // no. of equal intervals per octave (E tuning)
            cv = _floor*(1+random(0,n_scale+1)/n_scale); 
            break;}     
          case 'P': { // Predefined Scale
            cv = _floor*scale[random(0,s_length)]; 
            break;}       
          case 'R': { // Random
            double low = 1000.*log(_floor)/log(2); 
            double high = 1000.*log(_roof)/log(2);
            cv = pow(2,random(low,high++)/1000.);
            break;}
          }
          
          // randomly spread across octave registers
          cv = FixInterval( cv*pow(2,random(0,log(_roof/_floor)/log(2))));

          break;
        }
    
      case 'H': { cv = HarmonicInterval(); break; }
    
      case 'I': { cv = Interval2*pow(InharmonicInterval(),pow(-1,random(0,2))); break; }

      case 'K': {
        // Use with SingeVoice mode, follows CV input
        UpdateInterval();
        cv = _floor+ReadLeft()*(_roof-_floor);
        break;}
    
      case 'G': { cv = scale[_n]*pow(2,_channel); SubD = 16/pow(2,_channel+1);
                _n += random(-1,2);  _n = constrain(_n,0,s_length-1);  break; }
                                                                    
      case 'R': {
        double freq = 2.05*0.25/((base_period)/1000.);
        if (InputMode) { freq *= ReadLeft(); }
        // if (freq == 0) {freq = 1; Serial.println("avoid");}
        cv = _floor*pow(2,Saw(log(_roof/_floor)/log(2), freq , 1*_channel/4. , 1+_channel%1)); 
        break;}
                          
      case 'S': { 
        Serial.println(s_length);
        int n_low,n_high;
        if (Strata) {
          n_low = _channel*s_length;
          n_high = (_channel+1)*s_length;
        }
        else {
          n_low = 0; n_high = 4*s_length;
        }
        _n = constrain(_n,n_low,n_high); 
        cv = scale[_n%s_length];//*pow(2,_channel); //SubD[channel] = 16/pow(2,channel+1); 

        cv *= pow(2,(_n)/s_length); // set octave register

        if (_n<=n_low) {_Dir = 1;}//_n=low+1;}
        if (_n>=n_high) {_Dir = 0;}// _n=high-1;}
        int r = random(1,1000);
        int stepsize = 1;
        if (RandomStep(0.3)) {stepsize++;}
        if (_Dir) {_n+=stepsize;}
        else      {_n-=stepsize;}
        
        //if (_channel!=0) {cv = 4;}
        break; }
                          
      case 'F': { cv = Interval2*pow(1.5,pow(-1,random(0,1))); break; }

      // Just Intonation + Drone
      case 'X': {
        if (_channel == 0) { cv = 4*JustInterval(7,1)*JustInterval(7,1); }
        else { cv = 4;}
        break;}

    }

    cv = FixInterval(cv);
  }

  else if (RingChanges && _counter>=_N) {
    if (_counter%_N==0) {
      Swap(_CV, _iteration, _N);
      _iteration++;
    }
    cv = _CV[_pos];
  }

  // else if (!ShiftRegisterMode && !RingChanges){ //SingleVoice mode
  else if (SingleVoice){ // channel!=0
    // only works perfectly when ChangeTrigOrder = 0;

    // Allow user input of FM Carrier:Modulator Ratio on the fly via Serial
    if (Serial.available() && _channel==3) {_Interval = Serial.parseFloat();}

    //cv = Voices[0].GetCV(0)*_Interval;

    cv = Voices[0]->GetCV(0)*_Interval;

    if (cv < 1) {cv = 1;} // Exception for rounding error
  }
  
  // else if (ShiftRegisterMode && !RingChanges) {
  else if (ShiftRegisterMode) {

    cv = ApplyShiftRegister();
  }

  _CV[_pos] = cv;

  _pos = (_pos+1)%_N;
  _counter++;
  // if (ProximityMode) {EvaluateProximity();}
//}
}

double RSS::FixInterval( double Interval) { 
  if (Interval == 0) {Interval = _floor;}
  //if (Interval > r) {Interval = Interval - 2*(Interval-r);}
  //if (Interval < f) {Interval = Interval + 2*(f-Interval);}

  while (Interval > _roof) {Interval /= _roof/_floor;}
  while (Interval < _floor) {Interval *= _roof/_floor;}
  //Serial.print("Interval"); Serial.print("\t"); Serial.println(Interval);
  return Interval;
}

double RSS::JustInterval( int base , int power) {
  double interval = pow(random(1,base+1),random(1,power+1))/pow(random(1,base+1),random(1,power+1));
  //if (interval < 1) {interval = 1./interval;}
  //return interval; }
  double r = random(2,8);
  r /= (r-1.);
  int posneg = pow(-1,random(-1,1));
  r = pow(r,posneg);
  return r;}

double RSS::HarmonicInterval() {
  // Generates Harmonic Intervals in ascending and then descending order
  double interval = _HarmonicCount;
  _HarmonicCount += random(1,2) ; 
  if (_HarmonicCount>0.5*_roof/_floor) { _HarmonicCount -= (int) _roof/(2*_floor) ;}
  interval = random(112,145)/64.;
  return interval;
}

double RSS::InharmonicInterval() {
  // here for now, will probably deprecate before long
  const int offset = 9; // move and make static?
  const int primes[] = {5,7,11,13,17,19,23,29,31,37,41,43,47,51,53,57}, p_length = 5;
  static int p, lp[4];   int ind = random(0,p_length);
  while (ind == lp[0] || ind == lp[1] || ind == lp[2] || ind == lp[3]) { ind = random(0,p_length); }
  lp[p%4] = ind; p++;
  if (Debug) {Serial.print("index: "); Serial.println(ind);}
  return ( pow(2,(primes[offset+ind]+1.)/primes[offset+ind])); }
  //return (primes[offset+ind]/47.); }
  
double RSS::DifferenceInterval(double Interval, int channel) {
// DifferenceInterval generates a set of four pitches that yield a single fundamental difference tone

  // Register determines the pitch of the cluster
  double Register = 10; // value should like between 2 and 15 ( ? )
  
  // Increment determines the frequency of the difference tone, as a fraction of the cluster pitch
  // Smaller values result in lower difference tones, larger values result in higher difference tones
  double Increment = Register/10.; 
  double t = micros();
  double Sweep = 120.*1000000.;
  Increment *= t/(Sweep);
  if (t>Sweep){_Play = 0;}
  
  // The pitch of the difference tone can be modulated by an LFO
  boolean ModulateDifferenceTone = 0;
  Increment += ModulateDifferenceTone * LFO(-0.1,0.05,0,1);
  

  if (InputMode) {Increment *= ReadLeft();}
  
  /* Register can be modulated while keeping Increment the same to yield a pitch cluster that changes
  while the difference tone remains the same */
  boolean ModulateRegister = 0;
  Register += ModulateRegister * LFO(Register/4,0.01,0,1);
   
  /* The difference tone is generated by constructing a tone cluster in which 
  each tone is offset by equal amounts */
  Interval = Register+Increment*channel*1;
  //Register = Interval;
  return Interval;
}

double RSS::Saw( double amplitude , double frequency , double phase , bool rising) {
  unsigned long period = 1000./frequency; // [ms]
  unsigned long result = millis();
  result += phase*double(period);
  while (result > period) {result -= period;}
  if (!rising) {result = period-result;}
  return amplitude * result/period;
}

// double RSS::LFO( double amplitude , double frequency , double phase , bool sine ) {
//   double result = amplitude*sin(frequency*millis()/1000.*2.*pi+phase*pi);         // sine-wave LFO
//   if (!sine) { if ( result < 0 ) {result = 0;} else { result = amplitude; } }    // half-wave square
//   return result; 
// }

double RSS::GetCV(int i) {
  return _CV[(_pos+i)%_N];
}

bool RSS::GetTrig( int t ) {
  return _TrigArray[t%_N];
}

void RSS::Transpose(double TransposeAmount) {
  for (int i=0; i<_N; i++) { 
    _CV[i] = FixInterval(_CV[i]*TransposeAmount);} }

void RSS::OrderSequence() {
  // NEEDS REWRITE
  bool recurse = false; double temp;
  for (int i=1; i<_N; i++) { 
    if (_CV[i]<_CV[i-1]) {
      temp = _CV[i]; _CV[i] = _CV[i-1]; _CV[i-1] = temp; recurse = true; } }
  if (recurse) {OrderSequence();} }

void RSS::Reverse() {
  double temp;
  for (int i=0; i<(_N-1)/2; i++) { 
    temp = _CV[i]; 
    _CV[i] = _CV[_N-(i+1)]; 
    _CV[_N-(i+1)] = temp; } }

void RSS::Invert() {
  for (int i=0; i<_N; i++) {
    _CV[i] = _roof/_CV[i]; } }

void RSS::FlipT() { 
  for (int i=0; i<_N; i++) { 
    _TrigArray[i] = !_TrigArray[i];} 
  }

void RSS::Shift( int Tshift , int Nshift ) {
  bool TempTrigArray[_N]; double TempNoteArray[_N];
  for (int i=0; i<_N; i++){
    TempTrigArray[i] = _TrigArray[i];
    TempNoteArray[i] = _CV[i];
  }

  for (int i = 0; i < _N;  i++) {
    _TrigArray[i] = TempTrigArray[(i+Tshift+_N) % _N];
    _CV[i]        = TempNoteArray[(i+Nshift+_N) % _N]; }
 }

void RSS::SetNoteDuration(int amt){
  _NoteDuration = amt;
  // EXCEPTION
  if (SubD!=1 && _Rule == 'E'){
    _NoteDuration = 0;
  }
}

double RSS::ArrayMax(double * array) {
  double m = array[0];
  for (int i = 1; i < 4; i++) {  m = max(m,array[i]);  }
  return m;
}

double RSS::ArrayMin(double * array) {
  double m = array[0];
  for (int i = 1; i < 4; i++) {  m = min(m,array[i]);  }
  return m;
}

void RSS::sendNote(int t) {
  double val = _CV[t%_N];
  // overwrite t dependence
  val = _CV[(_pos-1+_N)%_N];
  //_CV[_pos%(_N-1)]??
  
  // Apply Root Shift
  if (ChangeRoot) {
    val *= Root; // need to modify roof to respect root?
    val = FixInterval(val);
  }

  static unsigned long ClockStart;
  if (Quantize) {
    val = QuantizeNote(val);
  }
  double mV = 1000*log(val)/log(2);                // frequency ratio to CV [mV]
  double mV0 = _prevCV;   // previous CV [mV]
  //mV0 = mV-150; // temp
  
  // DelayLength defines the note duration
  unsigned long DelayLength = 1000*base_period; // [us]
  if (_D) { DelayLength = ModDelay( mV );}

  if (Phasing) {DelayLength += LFO(DelayLength,0.01,_channel-1,HIGH); mV += 0*LFO(100,0.005,_channel,HIGH);}

  // GlideLength defines the glide length, if glide is on
  unsigned long GlideLength = DelayLength/2.; // [us]
  //GlideLength = Portamento;
  if (_Rule == 'B') {
    // GlideLength *= 1.25;
    // Serial.println(GlideLength);
    // GlideLength = max(GlideLength,300000);
    // if (_channel==0) {
      GlideLength = 1000*random(100,610);
      //} 
  }
  double GlideRate = (mV-mV0)/GlideLength;  // [mV/us]
  unsigned long BounceTime = 13000+0.15*pow(GlideLength/1000000.,1.5)*1000000.; // changed back from 13000 to 3000?
  //if (Rule == 'A') {BounceTime *=random(1,1)/5.;}
  int i = -1;  unsigned long t_start = micros();
  if (!B && !G) { GlideLength = 0; }
  if (!G) { mV0 = mV; GlideRate = 0; } 
  
  unsigned long CalibrationAmt = micros()-ClockStart;
  unsigned long tStart, GlideStart=micros();
  while ( i == -1 || (i*100 < long(GlideLength) && BounceTime > 6000 ) ) {
    tStart = micros();
    if ( i == -1 || (B && (micros()-t_start)>=BounceTime && BounceTime > 3000) ) {
      sendTrigger(); BounceTime = pow(BounceTime,0.99-LFO(0.0,0.01,1,1)); t_start = micros();
      //sendTrigger(channel); BounceTime = pow(BounceTime,1.01-LFO(0.0,0.01,1,1)); t_start = micros();
    }
    if ((i+1)*100 >= GlideLength) { mV0 = mV;}
    else { mV0 += double(GlideRate*100.);}
    int b = mV0;
    if (Detune)    {b += random(-20,21);}
    //if (_channel==0 && _Rule=='B') {b=random(3800,4001);}
    if (Flicker)   {b = flicker(b);}
    Send2DAC(_channel,b);
 
    if (i>-1) {delayMicroseconds(100-(micros()-tStart)); }
    i++; 
  }
  
  _prevCV = mV;
  if (BZZZ) {return;}
  if (Debug) {Serial.print("Voice "); Serial.print(_channel); Serial.print(": "); Serial.println(val); }

  GlideLength = micros()-GlideStart;
  CalibrationAmt += GlideLength;
  
  // Adjust the post-trigger delay to account for computation time
  if (CalibrationAmt < DelayLength) {
     DelayLength = DelayLength-CalibrationAmt;
  }
         
  //if ( G[channel] && DelayLength>GlideLength) {DelayLength -=GlideLength;}
  //else if (G[channel]) {Serial.println('z'); DelayLength = 0;}
  
  if (_D) {
      usDelay(DelayLength);   
  }

  ClockStart = micros();
}


void RSS::sendTrigger() {    
// send longer pulse, i.e. delayMicroseconds(300), for pinging filters (~300us) and LPGs (~200us) ??
  if (_counter%TriggerSubdivisions!=0) {return;}

  if (ProximityMode){
    EvaluateProximity();

    if (!_Prox) {
      // in case Trig level is high, write low
      digitalWrite(TrigPins[_channel],LOW);
      return;
    }
  }
  //if (SubD[channel] != 1 && Rule == 'E') {Duration = 0;} // EXCEPTION
  _NoteDuration = NoteDuration;
  if (FillIn) { 
    int i = 2; 
    while (_TrigArray[_pos+i]==0 && knobswitch) {//????????
      _NoteDuration++; i++;
    }
  }
  digitalWrite(TrigPins[_channel], HIGH);

  if (ProximityMode) {
    if (_Prox) {_NoteDuration=1;}
  }
  if (_NoteDuration<1) { digitalWrite(TrigPins[_channel],LOW); }
  if (!Sync) { DepleteNoteDuration(); } 
}

double RSS::QuantizeNote( double Note ) {
  int i = 0;
  while (Note>=2) {Note /= 2; i++; }
  
  int j=0;
  while (Note > scale[j] && j<s_length-1) {j++;}
  if (scale[j]-Note > Note-scale[j-1]) {Note = scale[j-1];}
  else {Note = scale[j];}
  return Note*pow(2,i);
}

// unsigned long RSS::ModDelay(boolean Delay , double millivolts) {
//   //(1,1) arguments yield no modification to delay
//   unsigned long DelayLength = Delay*1000*base_period;

//   DelayLength *= pow( max(base_period*(1+LFO(a2,f2,0,HIGH)+LFO(a1,f1,0,HIGH)) , 1) / base_period, D2) * pow(millivolts/2000,D1);
//   if (InputMode) { DelayLength *= analogRead(LeftIn)/1023.; }
//   return DelayLength;  
// }

int RSS::flicker( int note ) {
  static boolean flicker_switch;
  if (flicker_switch) {note+=1000;}
  flicker_switch = !flicker_switch;
  return note;
}

bool RSS::Play(int t){
  if (t%SubD == 0 && _Play) {
    return true;
  }
  else{
    return false;
  }
}

double RSS::ApplyShiftRegister(){
// maybe poorly defined for ChangeTrigOrder mode  
  // determine current voice
  bool found = false;
  int i = 0;
  while (!found) {
    if (_channel == TrigOrder[i]){
      found = true; break;
    }
    if (!found) {i++;}
  }

  // define previous voice
  int j = TrigOrder[i-1];

  // get shift register from previous voice
  return Voices[j]->GetCV(-2);
}

void RSS::EvaluateProximity(){
// evaluate proximity of voice's pitch to that of the other voices
  double ProximityRange = 0.05;
  bool WithinProximity = false;
  for (int i = 0; i < 4; i++) {
    if (i != _channel) {
      if ( abs( _CV[_pos] - Voices[i]->GetCV(0) ) < ProximityRange ) {
        WithinProximity = true;
      }
    }
  }
  if (WithinProximity) {_Prox = true; }
  else {_Prox = false;}
}

void RSS::RandomizeInitialOrder() {
  int a, b;
  for (int i=0; i<_N/2; i++) {
    a = random(0,_N);
    b = random(0,_N);
    while (a==b) { b = random(0,_N); }
    double temp = _CV[a];
    _CV[a] = _CV[b];
    _CV[b] = temp;
  }
}

bool RSS::RandomStep(double probability) {
  int r = random(0,1001);
  if (r<1000*probability) {
    return true;
  }
  else { 
    return false;
  }
}

void RSS::UpdateInterval() {
  bool ModifyStack = true;

  if (ModifyStack) {

  // if (SingleVoice) {
    //double Formants[4] = {1, 2300./270, 3000./270, 3000./270};
    //double Formants[4] = {700,760,700,760};
    // double Formants[4] = {1, 870/300, 2250/300, 1};
    // double Formants[4] = {1, 2000/400, 2550/400, 1};
    //double Formants[4] = {1, 2000/400, 2550/400, 1};
    //double Formants[4] = {1, 1700/660, 2400/660, 2400/660};
    //double Formants[4] = {1.618,1.618,1,1};

    double offset = 112;
    double HarmonicSeparation = 1;
    if (InputMode) {
      HarmonicSeparation = int(10*ReadLeft());
      // Serial.print("HarmonicSeparation: ");Serial.println(HarmonicSeparation);
    }
    for (int i = 0; i < 4; i++ ) {
      // Stack[i] = offset+HarmonicSeparation*i;
      int center = 128;
      int degree = 16;
      int harmonic = random(center-degree,center+degree+1);
      //Stack[_channel] = random(112,145);
      while (harmonic == Stack[0] || harmonic == Stack[1] || harmonic == Stack[2] || harmonic == Stack[3]) {
        harmonic = random(center-degree,center+degree+1);
      }
      Stack[i] = harmonic;
    }
  }

  double StackMin = ArrayMin(Stack);
  // _Interval = Stack[_channel] / Stack[0];
  for (int i = 0; i<4; i++){
    Voices[i]->_Interval = Stack[i] / Stack[0];
  }

  if (_channel==0) {
    _roof = 16.*(Stack[0]/ArrayMax(Stack));
    _floor = (Stack[0]/ArrayMin(Stack));
    // _floor = 1./range;
  }
  // Serial.print("test"); Serial.println(_Interval);
  // }
}

