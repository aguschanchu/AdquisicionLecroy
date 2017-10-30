#if ARDUINO >= 100
  #include "Arduino.h"
#else
  #include "WProgram.h"
#endif

#include "PID_AutoTune_v0.h"


PID_ATune::PID_ATune(double* Input, double* Output, double* Setpoint)
{
	input = Input;
	output = Output;
  setpoint = *Setpoint;
	controlType = 1; //Es para que trabaje con PID, en lugar de PI
	noiseBand = 0.25;
	running = false;
	oStep = 255;
	SetLookbackSec(10);
	lastTime = millis();
}



void PID_ATune::Cancel()
{
	running = false;
}

int PID_ATune::Runtime()
{
	justevaled=false;
	if(peakCount>9 && running)
	{
		running = false;
		FinishUp();
		return 1;
	}
	unsigned long now = millis();
  //sampleTime=100;
	Serial.println(sampleTime);
	if((now-lastTime)<sampleTime) return false;
	lastTime = now;
	double refVal = *input;
	justevaled=true;
	if(!running)
	{ //initialize working variables the first time around
		peakType = 0;
		peakCount=0;
		justchanged=false;
		absMax=refVal;
		absMin=refVal;
		running = true;
    estable = false;
		outputStart = *output;
		*output = outputStart+oStep;
	}
	else
	{
		if(refVal>absMax)absMax=refVal;
		if(refVal<absMin)absMin=refVal;
	}

	//oscillate the output base on the input's relation to the setpoint
	if(refVal>setpoint+noiseBand && *output-oStep > -256) {
	  *output = outputStart-oStep;
	}
	else if (refVal<setpoint-noiseBand && *output+oStep<256) {
	  *output = outputStart+oStep;
	}

  // Como a veces se utiliza un setpoint de autotuneo distinto al valor inicial, eso da un valor de amplitud incorrecta. Ver changelog, corrida 19
  // Por eso, medimos la amplitud, a partir de tener dos picos (donde esperamos que el sistema ya se haya estabilizado)
  if (peakCount>2 && estable == false) {
    absMax=refVal;
    absMin=refVal;
    estable = true;
  }

  //bool isMax=true, isMin=true;
  isMax=true;isMin=true;
  //id peaks
  for(int i=nLookBack-1;i>=0;i--)
  {
    double val = lastInputs[i];
    if(isMax) isMax = refVal>val;
    if(isMin) isMin = refVal<val;
    lastInputs[i+1] = lastInputs[i];
  }
  lastInputs[0] = refVal;
  if(nLookBack<9)
  {  //we don't want to trust the maxes or mins until the inputs array has been filled
	return 0;
	}

  if(isMax)
  {
    if(peakType==0)peakType=1;
    if(peakType==-1)
    {
      peakType = 1;
      justchanged=true;
      peak2 = peak1;
    }
    peak1 = now;
    peaks[peakCount] = refVal;

  }
  else if(isMin)
  {
    if(peakType==0)peakType=-1;
    if(peakType==1)
    {
      peakType=-1;
      peakCount++;
      justchanged=true;
    }

    if(peakCount<10)peaks[peakCount] = refVal;
  }

  if(justchanged && peakCount>5)
  { //we've transitioned.  check if we can autotune based on the last peaks
    double avgSeparation = (abs(peaks[peakCount-1]-peaks[peakCount-2])+abs(peaks[peakCount-2]-peaks[peakCount-3])+abs(peaks[peakCount-3]-peaks[peakCount-4])+abs(peaks[peakCount-4]-peaks[peakCount-5])+abs(peaks[peakCount-5]-peaks[peakCount-6]))/5;
    if( avgSeparation < 0.05*(absMax-absMin))
    {
		FinishUp();
      running = false;
	  return 1;

    }
  }
   justchanged=false;
	return 0;
}
void PID_ATune::FinishUp()
{
	  *output = outputStart;
      //we can generate tuning parameters!
      Ku = 4*(2*oStep)/((absMax-absMin)*3.14159);
      Pu = (double)(peak1-peak2) / 1000;
      Serial.println("Ahi van los Ku y Pu");
      Serial.println(Ku);
      Serial.println(Pu);
}

double PID_ATune::GetKp()
{
	return controlType==1 ? 0.6 * Ku : 0.4 * Ku;
}

double PID_ATune::GetKi()
{
	return controlType==1? 1.2*Ku / Pu : 0.48 * Ku / Pu;  // Ki = Kc/Ti
}

double PID_ATune::GetKd()
{
	return controlType==1? 0.075 * Ku * Pu : 0;  //Kd = Kc * Td
}

void PID_ATune::SetOutputStep(double Step)
{
	oStep = Step;
}

double PID_ATune::GetOutputStep()
{
	return oStep;
}

void PID_ATune::SetControlType(int Type) //0=PI, 1=PID
{
	controlType = Type;
}
int PID_ATune::GetControlType()
{
	return controlType;
}

void PID_ATune::SetNoiseBand(double Band)
{
	noiseBand = Band;
}

double PID_ATune::GetNoiseBand()
{
	return noiseBand;
}

void PID_ATune::SetLookbackSec(int value)
{
    if (value<1) value = 1;

	if(value<25)
	{
		nLookBack = value * 4;
		sampleTime = 250;
	}
	else
	{
		nLookBack = 100;
		sampleTime = value*10;
	}
}

int PID_ATune::GetLookbackSec()
{
	return nLookBack * sampleTime / 1000;
}
