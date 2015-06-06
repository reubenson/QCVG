#include "Arduino.h"
#include "Scale.h"
#include "QCVG_MAIN.h"

int Scale::length;
double *Scale::s;

void Scale::GenerateScale() {
// Generate scale according to Tuning, NumberOfDivisions, and ScaleChoice

  switch (Tuning) {
    case 'P': {
      s = new double[NumberOfDivisions];
  
      // generate Equal Temperament Scale
      double S[NumberOfDivisions];
      for (int i=0; i<NumberOfDivisions+1; i++) {
        S[i] = pow(2,double(i)/double(NumberOfDivisions));
        s[i] = S[i];
      }
      length = NumberOfDivisions;

      if (NumberOfDivisions==12) {
        switch (ScaleChoice) {

          // Modcan quantizer scale banks below
          case -1:  { Build(S,0,2,4,5,7,9,11); break; }          // Major
          case -2:  { Build(S,0,2,3,5,7,8,10); break; }          // Minor
          case -3:  { Build(S,0,2,3,5,7,9,10); break; }          // Dorian
          case -4:  { Build(S,0,1,3,5,7,8,10); break; }          // Phyrgian
          case -5:  { Build(S,0,2,4,6,7,9,11); break; }          // Lydian
          case -6:  { Build(S,0,2,3,6,7,8,11); break; }          // Aeolian
          case -7:  { Build(S,0,2,4,5,7,9,10); break; }          // Mixolydian
          case -9:  { Build(S,0,3,5,6,7,10);   break; }          // Blues
          case -10: { Build(S,0,2,3,5,6,8,9,11); break; }        // Diminished
          case -11: { Build(S,0,3,4,6,8,11); break; }            // Augmented
          case -12: { Build(S,0,2,5,7,10); break; }              // Pentatonic Neutral
          case -13: { Build(S,0,5,10); break; }                  // Fourths

          case -14: { Build(S,0,4,7); break;}                    // Major
          case -15: { Build(S,0,4,7,9); break;}                  // Major6
          case -16: { Build(S,0,4,7,11); break;}                 // Major7
          case -17: { Build(S,0,4,6,11); break;}                 // Major7b5
          case -18: { Build(S,0,3,7); break;}                    // Minor
          case -19: { Build(S,0,3,7,9); break;}                  // Minor6
          case -20: { Build(S,0,3,7,10); break;}                 // Minor7
          case -21: { Build(S,0,5,7); break;}                    // Sus4
          case -22: { Build(S,0,2,7); break;}                    // Sus2
          case -23: { Build(S,0,2,5,7); break;}                  // Sus4 Sus2
          case -24: { Build(S,0,4,8); break;}                    // Augmented
          case -25: { Build(S,0,3,6); break; }                   // Diminshed
          case -26: { Build(S,0,3,6,9); break; }                 // Diminished 7th
          case -27: { Build(S,0,5,7,10); break; }                // 7 Sus4
          case -28: { Build(S,0,3,6,10); break;}                 // Min 7 b5
          case -29: { Build(S,0,3,8,10); break;}                 // Min 7 b5
  
          case -30: { Build(S,0,2,3,5,6,7,8,11); break;}         // Algerian
          case -31: { Build(S,0,1,3,4,6,8,10); break; }          // Altered
          case -32: { Build(S,0,2,3,5,6,8,9,11); break;}         // Aux Diminished
          case -33: { Build(S,0,1,3,7,8); break;}                // Balinese
          case -34: { Build(S,0,1,4,5,7,8,11); break;}           // Byzantine
          case -35: { Build(S,0,2,4,7,9); break;}                // Diatonic
          case -36: { Build(S,0,1,4,5,7,8,10); break;}           // Byzantine
          case -37: { Build(S,0,1,4,5,7,8,11); break;}           // Double Harmonic
          case -38: { Build(S,0,2,4,5,7,8,10); break;}           // Hindu 
          case -39: { Build(S,0,3,4,5,8,9); break;}              // Sixtone Symmetric
          case -40: { Build(S,0,2,3,4,6,7,8,9,11); break; }      // Nine Tone
          case -41: { Build(S,0,2,4,6,7,9,10); break;}           // Overtone Dominant
          case -42: { Build(S,0,1,3,7,8); break;}                // Pelog
          case -43: { Build(S,0,2,4,6,9,10); break;}             // Prometheus
          case -44: { Build(S,0,1,4,6,8,10,11); break;}          // Enigmantic
          case -45: { Build(S,0,1,3,4,6,7,9,10); break;}         // Octatonic    

          case 0: {  s[0]=1.; s[1]=16./15 ; s[2]= 6./5  ; s[3]=3./2     ; s[4]=8./5; length = 5; break; } // pelog (modcan) 
          case 1: {  s[0]=1.;s[1]=9./8    ; s[2]=5./4   ; s[3]=3./2     ; s[4]=5./3; length = 5; break; } // Sunda
          case 2: { s[0]=1.;s[1]=35./32   ; s[2]=5./4   ; s[3]=21./16   ; s[4]=49./32   ; s[5]=105./64 ; s[6]=7./4; length = 7; break; } // HMC American Gamelan
          case 3: { s[0]=1.;s[1]=8./7     ; s[2]=12./7  ; s[3]=21./17   ; s[4]=147./128 ; length = 5; break; } // Balungan
          case 4: { s[0]=1.;s[1]=567./512 ; s[2]=9./8   ; s[3]=147./128 ; s[4]=21./16   ; s[5]=189./128; // (LMY 'Well' Tuning)  
              s[6]=3./2     ; s[7]=49./32 ; s[8]=7./4   ; s[9]=63./32   ; length=10;  break;}
        }
      }
      break;
    }

    case 'H': {
      double LowestHarmonic, HighestHarmonic;
      bool LinearMode = false;
      int center = 8, span = 4;
      if (LinearMode) {
        LowestHarmonic = center-span;
        HighestHarmonic = center+span;
      }
      else {
        int Register = 0;
        LowestHarmonic = pow(2,Register+0);
        HighestHarmonic = pow(2,Register+4);
      }
        
      length = HighestHarmonic-LowestHarmonic+1;
      s = new double[length];
      for (int i=0; i<length; i++) {
        s[i] = double(LowestHarmonic+i)/double(LowestHarmonic);
      }
      break;
    }
  }
}

void Scale::Build(double *S,int a,int b,int c,int d,int e,int f,int g,int h, int i){
// Subroutine to populate scale array, s, with the inputs
  int array[] = {a,b,c,d,e,f,g,h,i};
  
  int l=0;
  do { s[l] = S[array[l]];  l++; } 
  while ( array[l]>0 );
 
  // length of scale
  length = l;
}

void Scale::QuantizeNote( double & Note ) {
// Quantize Note to the scale defined

  // divide Note by 2 until between 1 and 2
  int i = 0;
  if (s[length-1] < 2) {
    while (Note>2) {Note /= 2; i++; }
  }

  // determine closest note on the scale to Note
  int j=1;
  while (Note > s[j] && j<length) {j++;}
  double NoteBelow = s[j-1], NoteAbove = s[j];
  // double NoteAbove = s[j];
  if (j == length) { NoteAbove=2.0; }

  if (NoteAbove/Note > Note/NoteBelow ) {
    // Note is closer to NoteBelow than NoteAbove
    Note = NoteBelow;
  }
  else  { 
    // Note is closer or equal to NoteAblove than NoteBelow
    Note = NoteAbove;
  }

  // restore Note to its original octave register
  Note *= pow(2,i);
}

void Scale::PrintScale() {
// Print scale to serial monitor
  Serial.print("Scale: ");
  for (int i=0; i<length; i++) {
    Serial.print(s[i]); Serial.print("  "); 
  }
  Serial.println();
}

double Scale::JustInterval( int base , int power) {
// Generate just interval (whole number ratio) 
  return pow(random(2,base+1),random(1,power+1))/pow(random(2,base+1),random(1,power+1));
}

double Scale::Get(int i){
// Return pitch interval corresponding to input
// Go up in register if input is larger than length of scale  
  int Factor = 0;
  while (i >= length) { i -= length; Factor++; }
  return s[i]*pow(2,Factor);
}

double Scale::HarmonicInterval(double LowestHarmonic , double HighestHarmonic , double base) {
// Return intervallic value based on the harmonic series

  bool RisingMode = true;
  // RisingMode: ascending and descending harmonic series
  // !RisingMode: random selection from harmonic series

  // keep track of current set of harmonics
  static double memory[4];
  static int index=0;

  if (RisingMode) {
    static bool Direction;
    // Generates Harmonic Intervals in ascending and then descending order
    memory[index%4] = memory[(index-1)%4]+pow(-1,Direction);
    if (memory[index%4]<LowestHarmonic) {
      memory[index%4] = LowestHarmonic;
      Direction = 0;
    }
    else if (memory[index%4]>HighestHarmonic) {
      memory[index%4] = HighestHarmonic;
      Direction = 1;
    }
    memory[index%4] = constrain(memory[index%4],LowestHarmonic,HighestHarmonic);
  }
  else {
    // randomly select from harmonic series
    int harmonic = random(LowestHarmonic,HighestHarmonic+1);
    // ensure that a new value is chosen
    while (harmonic == memory[0] || harmonic == memory[1] || harmonic == memory[2] || harmonic == memory[3]) {
      harmonic = random(LowestHarmonic,HighestHarmonic+1);
    }
    memory[index%4] = harmonic;
  }

  // calculate intervallic ratio, relative to the lowest harmonic
  double interval = memory[index%4]/LowestHarmonic;

  // scale to base
  double maxvalue = 16.*(LowestHarmonic/HighestHarmonic);
  interval *= (1-base)+(maxvalue)*(base);

  index++;
  return interval;
}

double Scale::DifferenceInterval(int channel, double mod, double span) {
// Generate intervals that yield a single fundamental difference tone
  
  // Increment determines the frequency of the difference tone
  // smaller values result in lower difference tones, larger values result in higher difference tones
  // Increment==0 yields unison tuning
  double Increment = 2.0;
  // multiply by scaling factor (0 to 1)
  Increment *=span;
  
  // Register determines the pitch of the tone cluster. Modulation of register yields a changing
  // pitch cluster with a constant difference tone (i.e. low frequency sine wave to mod input)
  double ModulationRange = 15.0-3.0*Increment;
  double RegisterCenter = 1.+0.5*(ModulationRange);
  double Register = RegisterCenter+0.5*ModulationRange*mod;
   
  // The difference tone is generated by constructing a tone cluster in which 
  // each frequency is offset by an equal amount, Increment.
  return Register+Increment*channel;
}
