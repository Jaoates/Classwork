#include <Servo.h>

// set up vars
int in;
int i = 0;
Servo servo;

// set up constants
const int servopin = 9;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600); //set up serial
  pinMode(servopin,OUTPUT); //set up pins
  servo.attach(servopin);// init servo
}

void loop() {
  // put your main code here, to run repeatedly:
  if(Serial.available()>0){ //enter if there is stuff in the buffer
    in = Serial.parseInt(); // read an int from the serial buffer
    Serial.parseInt(); // read an extra int, for some reason theres an extra 0
    Serial.println(in); // print it for debugging
  }
  servo.write(in); // command the servo
}
