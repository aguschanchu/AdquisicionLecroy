#include <Wire.h>
#include "Adafruit_MCP9808.h"
#include "PID_v1.h"
#include "PID_AutoTune_v0.h"

////////Ajustes
double kp=302.61,ki=1.08,kd=21121.99;
double setpoint=5;
int temperaturaMaxima=71;
boolean tuning = false;
// Utiliza el resistor como perturbacion
boolean perturbador = false;
//Settings del perturbador
double dutyCyclePerturbacion = 0.4;
int periodoPerturbacion = 15; //en minutos
int esperaInicialPerturbado = 60; //Cuanto espera antes de iniciar el perturbador en minutos
//Settings de los pins del puente H
int pincalentado = 5;
int pinenfriado = 6;
//Numero del sensor ambiente. Contando desde 0.
int sensorAmbiente = 5;
////////Fin de ajustes
// Create the MCP9808 temperature sensor object
Adafruit_MCP9808 tempsensor = Adafruit_MCP9808();
byte ATuneModeRemember=2;
double input=25, output=0;
double outputStart=0;
double aTuneStep=255, aTuneNoise=0.2, aTuneStartValue=0;
unsigned int aTuneLookBack=20;
unsigned long ultimoPeriodo = -periodoPerturbacion*60*100;
unsigned long contadorEsperaInicial = 0;
bool dutyModificacion = false;
bool retenerPerturbador = false;
unsigned long  modelTime, serialTime, now, tiempo;
float ambiente;
double newSetpoint;

//Utilizado para la comunicacion
const byte buffSize = 40;
char inputBuffer[buffSize];
const char startMarker = '<';
const char endMarker = '>';
byte bytesRecvd = 0;
char messageFromPC[buffSize] = {0};
boolean readInProgress = false;
boolean newDataFromPC = false;
unsigned long prevReplyToPCmillis = 0;
unsigned long replyToPCinterval = 1000;
PID myPID(&input, &output, &setpoint,kp,ki,kd, DIRECT);
PID_ATune aTune(&input, &output, &setpoint);


void setup()
{
  int R_en = 2;
  int L_en = 1;
  pinMode(pincalentado, OUTPUT);
  pinMode(pinenfriado, OUTPUT);
  pinMode(7, OUTPUT);
  pinMode(R_en, OUTPUT);
  pinMode(L_en, OUTPUT);
  digitalWrite(R_en, 1);
  digitalWrite(L_en, 1);
  digitalWrite(7, 1);

  //Setup the pid
  myPID.SetMode(AUTOMATIC);

  if(tuning)
  {
    tuning=false;
    changeAutoTune();
    tuning=true;
  }

  serialTime = 0;
  Serial.begin(9600);

}

void loop()
{
  now = millis();
  // Reloj en segundos
  tiempo = now/1000;

  //pull the input in from the real world
  if (!tempsensor.begin(0x18)) {
    Serial.println("Sensor 0x18 no encontrado");
    while (1);
  }
  else {
    tempsensor.wake();
    double c;
     if (tuning) {
      // Hay que tomar menos puntos si promediamos, por el NLoookBack    if (tunning) {
        for (int j = 0; j < 10; j++) {
        c = c + tempsensor.readTempC();
        delay(100);
      }
      tempsensor.shutdown();
      input = c/10;
      }
    else {
      for (int j = 0; j < 500; j++) {
        c = c + tempsensor.readTempC();
        delay(20);
      }
      tempsensor.shutdown();
      input = c/500;
    }
    Serial.print("[7,");
    Serial.print(input);
    Serial.print(",");
    Serial.print(tiempo);
    Serial.print("]\n");
  }


  if(tuning)
  {
    byte val = (aTune.Runtime());
    //Si esta tuneado, no queremos perturbarlo, de modo que almacenamos
    //la intencion de hacerlo, en una variable, para despues reactivarlo
    if (perturbador){
      retenerPerturbador = true;
      perturbador = false;
    }
    if (val!=0)
    {
      tuning = false;
    }
    if(!tuning)
    { //we're done, set the tuning parameters
      kp = aTune.GetKp();
      ki = aTune.GetKi();
      kd = aTune.GetKd();
      Serial.println("Ahi van las salidas del autotune, como estado");

      Serial.print("['e',");
      Serial.print(tiempo);
      Serial.print(",");
      Serial.print(kp);
      Serial.print(",");
      Serial.print(ki);
      Serial.print(",");
      Serial.print(kd);
      Serial.print(",");
      Serial.print(setpoint);
      Serial.print(",");
      Serial.print("'autotune finalizado, PIDando'");
      Serial.print("]\n");

      myPID.SetTunings(kp,ki,kd);
      //Agrego un tiempo de control mayor a 0.1s
      myPID.SetSampleTime(1000);
      AutoTuneHelper(false);
      //Reactivamos el perturbador en caso de que haya que hacerlo
      //y reiniciamos el contadorEsperaInicial
      contadorEsperaInicial=now;
      if (retenerPerturbador){
        perturbador=true;
      }


    }
  }
  else myPID.Compute();


  //Ejecutamos una funcion que verifique que no se este friendo la caja
  if (input > temperaturaMaxima) {
    output=0;
  }

  //Aplicamos los cambios del PID
  if (output >= 0) {
  analogWrite(pinenfriado,0);
  analogWrite(pincalentado,output);
  }
  if (output < 0) {
  analogWrite(pincalentado,0);
  analogWrite(pinenfriado,-1*output);
  }

  //Para almacenar el valor de potencia
  Serial.print("[6,");
  Serial.print(output);
  Serial.print(",");
  Serial.print(tiempo);
  Serial.print("]\n");

  // Perturbamos a la caja con el resistor
  if (perturbador && (now-contadorEsperaInicial)>esperaInicialPerturbado*60*1000){
    if ((now-ultimoPeriodo)>=periodoPerturbacion*60*1000){
      CambiarVoltaje(6);
      dutyModificacion=false;
      ultimoPeriodo=now;
    }
    if ((now-ultimoPeriodo)>=periodoPerturbacion*60*1000*dutyCyclePerturbacion && dutyModificacion == false){
      CambiarVoltaje(3);
      dutyModificacion=true;
    }

  }
  //send-receive with processing if it's time
  if(millis()>serialTime)
  {
    SerialSend();
    //Enviamos OK a PC, y esperamos respuesta
    Serial.flush();
    Serial.println("READY");
    delay(100);
    getDataFromPC();
    update();
    delay(100);
    replyToPC();
    delay(100);
    serialTime+=60*1000;
  }

  // Enviamos la temperatura del resto de los sensores
  int numeroDeSensores = 6;
  for (int sensor = 24; sensor < 24 + numeroDeSensores; sensor++) {
    if (!tempsensor.begin(sensor)) {
      Serial.print("Sensor ");
      Serial.print(sensor-24);
      Serial.print("no encontrado\n");
      while (1);
    }
    else {
      tempsensor.wake();
      float c = tempsensor.readTempC();
      Serial.print("[");
      Serial.print(sensor-24);
      Serial.print(",");
      Serial.print(c);
      Serial.print(",");
      Serial.print(tiempo);
      Serial.print("]\n");
      tempsensor.shutdown();
    }
    delay(10);
  }

}

void changeAutoTune()
{
  if(!tuning)
  {
    //Set the output to the desired starting frequency.
    output=aTuneStartValue;
    aTune.SetNoiseBand(aTuneNoise);
    aTune.SetOutputStep(aTuneStep);
    aTune.SetLookbackSec((int)aTuneLookBack);
    AutoTuneHelper(true);
    tuning = true;
  }
  else
  { //cancel autotune
    aTune.Cancel();
    tuning = false;
    AutoTuneHelper(false);
  }
}

void AutoTuneHelper(boolean start)
{
  if(start)
    ATuneModeRemember = myPID.GetMode();
  else
    myPID.SetMode(ATuneModeRemember);
}


void SerialSend()
{
  Serial.print("setpoint: ");
  Serial.print(setpoint);
  Serial.print(" ");
  Serial.print("input: ");
  Serial.print(input);
  Serial.print(" ");
  Serial.print("output: ");
  Serial.print(output);
  Serial.print(" ");
  if(tuning){
    Serial.println("tuning mode");
  }
  else {
    Serial.print("kp: ");
    Serial.print(myPID.GetKp());
    Serial.print(" ");
    Serial.print("ki: ");
    Serial.print(myPID.GetKi());
    Serial.print(" ");
    Serial.print("kd: ");
    Serial.print(myPID.GetKd());
    Serial.println();
  }
}

void CambiarVoltaje(double voltaje)
{
  //Envia el voltaje en voltios al arduino, para que lo cambie en la fuente de DC
  Serial.print("['v',");
  Serial.print(voltaje);
  Serial.print(",");
  Serial.print(millis()/1000);
  Serial.print("]\n");

}


void getDataFromPC() {

    // La paja de enviar cosas al arduino, es que te toma de a bytes
    // Entonces, es necesario usar delimitadores para enviar mensajes mas largos
    // Basado en http://forum.arduino.cc/index.php?topic=225329.msg1810764#msg1810764

  if(Serial.available() > 0) {

    char x = Serial.read();

      // the order of these IF clauses is significant

    if (x == endMarker) {
      readInProgress = false;
      newDataFromPC = true;
      inputBuffer[bytesRecvd] = 0;
      parseData();
    }

    if(readInProgress) {
      inputBuffer[bytesRecvd] = x;
      bytesRecvd ++;
      if (bytesRecvd == buffSize) {
        bytesRecvd = buffSize - 1;
      }
    }

    if (x == startMarker) {
      bytesRecvd = 0;
      readInProgress = true;
    }
  }
}

//=============

void parseData() {

    // split the data into its parts

  char * strtokIndx; // this is used by strtok() as an index

  strtokIndx = strtok(inputBuffer,",");      // get the first part - the string
  strcpy(messageFromPC, strtokIndx); // copy it to messageFromPC

  strtokIndx = strtok(NULL, ","); // this continues where the previous call left off
  newSetpoint = atoi(strtokIndx);     // convert this part to an integer

}

//=============

void replyToPC() {

  if (newDataFromPC) {
    newDataFromPC = false;
    Serial.print("<Msg ");
    Serial.print("OK");
    Serial.println(">");
  }
}

//============

void update() {

   // this illustrates using different inputs to call different functions
  if (strcmp(messageFromPC, "setpoint") == 0) {
      // Nos dijeron de cambiar el setpoint, es valido?
     if (newSetpoint > -1 && newSetpoint < 10) {
       setpoint = newSetpoint;
       // Hay que cambiar las Kes, cual es la temp ambiente?
      if (!tempsensor.begin(sensorAmbiente+24)) {
         Serial.print("Sensor ambiente no encontrado\n");
         while (1);
       }
       else {
         tempsensor.wake();
         float ambiente = tempsensor.readTempC();
         tempsensor.shutdown();
       }
      kp=myPID.CalculaKp(setpoint,ambiente);
      ki=myPID.CalculaKi(setpoint,ambiente);
      kd=myPID.CalculaKd(setpoint,ambiente);
      myPID.SetTunings(kp,ki,kd);
   }
  }

}
