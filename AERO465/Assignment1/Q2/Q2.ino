
// assign the numeric value of the pins we will use for input and output
int inpin1 = 4;
int outpin1 = 10;

int inpin2 = 5;
int outpin2 = 2;

// setup the pins in their correct modes for I/O
void setup() {
  // put your setup code here, to run once:
  pinMode(inpin1,INPUT);
  pinMode(outpin1,OUTPUT);
  pinMode(inpin2,INPUT_PULLUP);
  pinMode(outpin2,OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  if(digitalRead(inpin1)){//read the value of input pin1
    digitalWrite(outpin1,1); //if it is high, set the value of output pin 1 to high
  }else{
    digitalWrite(outpin1,0);//else set it low
  }

  if(digitalRead(inpin2)){ // read the value of input pin2
    digitalWrite(outpin2,1);// if it is high set output pin 2 high
  }else{
    digitalWrite(outpin2,0);//else set it low
  }
}
