#ifndef Scale_h
#define Scale_h

extern char Tuning;
extern int ScaleChoice;
extern int NumberOfDivisions;

class Scale {
public:
	static void GenerateScale();
	static void QuantizeNote(double &);
	static void PrintScale();
	static double JustInterval(int, int);
	static double HarmonicInterval(double,double,double);
	static double DifferenceInterval(int, double, double);
	static double Get(int);
	static double *s;
	static int length;

private:
	static void Build(double*,int,int b=0,int c=0,int d=0,int e=0,int f=0,int g=0,int h=0,int i=0);
};

#endif