/* Example program for using TCS3200 color sensor
 * author: Panjkrc
 * date: 12/14/2019
 * url: https://github.com/Panjkrc/TCS3200_library
 */


#include <tcs3200.h>


const int pwmA = 5;
const int in1A = 7;
const int in2A =8;
 
const int pwmB = 6;
const int in1B = 9;
const int in2B = 10;

const int led = 13;

const char checkColor = 'b';
const int threshold = 10;


// Motor Speed Values - Start at zero

int MotorSpeed1 = 0;
int MotorSpeed2 = 0;

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
}

void loop() {
  //red = tcs.colorRead('r', 0);    //scaling can also be put to 0%, 20%, and 100% (default scaling is 20%)   ---    read more at: https://www.mouser.com/catalog/specsheets/TCS3200-E11.pdf
  //red = tcs.colorRead('r', 20);
  //red = tcs.colorRead('r', 100);

  color = tcs.colorRead('b');    //reads color value for blue
  if (color > threshold){
    tape = 1;
  }else{
    tape = 0;
  }
  // Serial.print("color = ");
  // Serial.print(color);
  
  digitalWrite(in1A, HIGH);
  digitalWrite(in2A, LOW);
  // Set Motor B forward
  digitalWrite(in1B, HIGH);
  digitalWrite(in2B, LOW);
  analogWrite(pwmA, tape*255);
  analogWrite(pwmB, tape*255);
  digitalWrite(led,tape);
  Serial.print(tape);
  Serial.println();

  delay(20);

}
