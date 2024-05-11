/* Credit for using TCS3200 color sensor Library:
 * author: Panjkrc
 * date: 12/14/2019
 * url: https://github.com/Panjkrc/TCS3200_library
 */


#include <tcs3200.h>

const char checkColor = 'b'; // choose what color to look for, we could use green tape by setting this to 'g'
const int threshold = 100; // choose a threshold, this takes a lot of tuning based on lighting conditions
const int speedScale = 200; //scale the speed of operation

const int pwmA = 5; // IO pin constants for Motor A
const int in1A = 7;
const int in2A = 8;
 
const int pwmB = 6; // IO pin constants for Motor B
const int in1B = 9;
const int in2B = 10;

const int led = 13; // IO pin constant for the LED

// init some vars to track state
int color;
bool tape;

// init the color sensor library of choice
tcs3200 tcs(0, 1, 2, 3, 4); // (S0, S1, S2, S3, output pin)  

void setup() {
  pinMode(pwmA, OUTPUT); // set up all the pins
  pinMode(pwmB, OUTPUT);
  pinMode(in1A, OUTPUT);
  pinMode(in2A, OUTPUT);
  pinMode(in1B, OUTPUT);
  pinMode(in2B, OUTPUT);
  pinMode(led, OUTPUT);
  // Set Motor A forward
  digitalWrite(in1A, LOW);
  digitalWrite(in2A, HIGH);
  // Set Motor B forward
  digitalWrite(in1B, HIGH);
  digitalWrite(in2B, LOW);
  // the highs and lows are flipped because they are installed backwards
  Serial.begin(9600); // init the serial library
}

void loop() {
  color = tcs.colorRead('b');    //reads color value for blue
  if (color > threshold){ // compare value to threshold
    tape = 1;
  }else{
    tape = 0;
  }

  if(tape){ // if the car detects that it is over tape
    analogWrite(pwmA, speedScale); // turn on Motor A
    analogWrite(pwmB, 0); // turn off Motor B
  }else{
    analogWrite(pwmB, speedScale); // turn on Motor B
    analogWrite(pwmA, 0); // turn off Motor A
  }

  digitalWrite(led,tape); // if the car detects that it is over tape, turn on the LED
  Serial.print(color); // for debugging
  Serial.println();

  delay(20); // a little break for the poor thing

}
