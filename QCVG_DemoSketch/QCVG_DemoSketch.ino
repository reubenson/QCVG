/*** Sketch Demonstrating Usage of the QCVG (Quad Control Voltage Generator) Module, Developed by Reuben Son  ***/
/***************  Audio Recording at https://soundcloud.com/reubenson/qcvg-demo-sketch   ******************/

#include <SPI.h>
#include <MemoryFree.h>
#include <QCVG_MAIN.h>
#include <QCVG_VOICE.h>
#include <Scale.h>

// ---------------------------------------- DEBUG MODE ---------------------------------------------------------------------
boolean Debug = 0;

// ---------------------------------------- SETUP QCVG & SET TEMPO ---------------------------------------------------------
QCVG_MAIN QCVG(100);                        // Input base note duration in ms

// ---------------------------------------- DEFINE FOUR VOICES -------------------------------------------------------------
char Rule = 'E';                            // Select rule for note generation, i.e., A,B,C,D,E,G,H,I,J,R,S
QCVG_VOICE Voice1(16,1,Rule,0);
QCVG_VOICE Voice2(16,8,Rule,1);
QCVG_VOICE Voice3(16,5,Rule,2);
QCVG_VOICE Voice4(16,16,Rule,3);
QCVG_VOICE *Voices[4] = {&Voice1,&Voice2,&Voice3,&Voice4};

// --------------------------------------- DEFINE MAIN PARAMETERS-----------------------------------------------------------
volatile boolean SwitchedOn = 1;            // 1: run immediately    0: run on KnobPress (i.e after tuning)
int ScaleChoice = -20;                      // Select Scale
char Tuning = 'P';                          // Select tuning: [P]re-defined, [H]armonic Series, [J]ustIntervals
int NumberOfDivisions = 12;                 // Number of equal-temperament notes in scale
boolean Quantize = true;                    // 1: Quantize notes to the scale    0: Notes not quantized
int TriggerLength = 0;                      // Trigger/Gate duration [ms]
double Glissando = 0.0;                     // Length of glissano as fraction of note length
int Repetitions = 1;                        // No. of times sequence repeats before generating new sequence
double TimingOffset = 0;                    // Amount of timing offset to voices
boolean SingleVoice = false;                // 1: Single Voice Mode     0: Default mode (four independent voices)
 int Voicing = 13;                          // Choose 4-note chord preset for SingleVoice mode
boolean ShiftRegisterMode = false;          // 1: Shift Register Mode   0: Default mode
boolean Strata = true;                      // 0: Each voice spans four octaves 1: Each voice occupies its own octave register
 int Layers[4] = {0,3,1,2};                 // Define the octave register that each of the four voices occupies
boolean RingChanges = false;                // 1: 'Plain Bob' method ringing applied to voices 0: Method ringing not applied
boolean ChangeTrigOrder = false;            // 1: Change TrigOrder after each cycle using method ringing
boolean Modulate = false;                   // 1: Apply LFO modulation to voice pitch    0: Default mode
boolean BZZZ = false;                       // 1: Erratic timer mode / under development    0: Default mode
boolean ET = false;                         // 1: Euclidean Transpose Toggle  0: Default mode
boolean ChangeRootMode = false;              // 1: Change root note, according to root_list  0: Default mode
double root_list[] = {4./3,81./64,9./8,1.}; // Order of root shift
boolean InputMode = true;                   // 1: Parametric input via input pins    0: Ignore input pins
boolean ProximityMode = false;              // 1: Voices triggered based on proximity in pitch to other voices
boolean Phasing = false;                    // 1: Phasing applied to note timing    0: Default mode

// ----------------------------------------- SETUP ---------------------------------------------------------------------

void setup() {
  Serial.begin(14400);   Serial.setTimeout(100); Serial.println("\nStart!"); 
  
  Scale::GenerateScale();
  Scale::PrintScale();
    
  // Tune oscillators to highest pitch
  if (!SwitchedOn) { Send2DAC(0,4000); Send2DAC(1,4000); Send2DAC(2,4000); Send2DAC(3,4000); }
 
}

// ----------------------------------------- MAIN LOOP -----------------------------------------------------------------

void loop() {
  QCVG.loop();
}


