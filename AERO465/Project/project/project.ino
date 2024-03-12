/* Example program for using TCS3200 color sensor
 * author: Panjkrc
 * date: 12/14/2019
 * url: https://github.com/Panjkrc/TCS3200_library
 */


#include <tcs3200.h>

const char checkColor = 'b';
const int threshold = 100;
const int speedScale = 200; //scale the speed of operation
const int lockoutMillis = 100; // a new turn will not start if it has been less than this time since the last turn

const int pwmA = 5;
const int in1A = 7;
const int in2A = 8;
 
const int pwmB = 6;
const int in1B = 9;
const int in2B = 10;

const int led = 13;

// Motor Speed Values - Start at zero

int MotorSpeed1 = 0;
int MotorSpeed2 = 0;
int lastTurn = 0;

// State
bool turnDir = 0; // if 0 it was last going left, if 1 it was going right

// int red, green, blue, white;
int color;
bool tape;

tcs3200 tcs(0, 1, 2, 3, 4); // (S0, S1, S2, S3, output pin)  //

void setup() {
  pinMode(pwmA, OUTPUT);
  pinMode(pwmB, OUTPUT);
  pinMode(in1A, OUTPUT);
  pinMode(in2A, OUTPUT);
  pinMode(in1B, OUTPUT);
  pinMode(in2B, OUTPUT);
  pinMode(led, OUTPUT);
  Serial.begin(9600);
  lastTurn = millis();
}

void loop() {
  color = tcs.colorRead('b');    //reads color value for blue
  if (color > threshold){
    tape = 1;
  }else{
    tape = 0;
  }
  // Serial.print("color = ");
  // Serial.print(color);
  
  digitalWrite(in1A, LOW);
  digitalWrite(in2A, HIGH);
  // Set Motor B forward
  digitalWrite(in1B, HIGH);
  digitalWrite(in2B, LOW);

  // if(turnDir){
  //   analogWrite(pwmA, speedScale);
  //   analogWrite(pwmB, 0);
  // }else{
  //   analogWrite(pwmB, speedScale);
  //   analogWrite(pwmA, 0);
  // }

  // if(!tape && (millis() - lastTurn) > lockoutMillis){
  //   lastTurn = millis();
  //   turnDir = !turnDir;
  // }

  if(tape){
    analogWrite(pwmA, speedScale);
    analogWrite(pwmB, 0);
  }else{
    analogWrite(pwmB, speedScale);
    analogWrite(pwmA, 0);
  }

  digitalWrite(led,tape);
  Serial.print(color);
  Serial.println();

  delay(20);

}
