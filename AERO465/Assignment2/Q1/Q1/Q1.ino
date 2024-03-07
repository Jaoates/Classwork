

const int LED1 = 6;

int counter = 0;
bool dir = 1;
void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(6,OUTPUT);

}

void loop() {
  // put your main code here, to run repeatedly:
  if(dir){
    counter++;
    }
  else{
    counter--;
    }
  if(counter == 0){
    dir = 1;
  }
  if(counter == 255){
    dir = 0;
  }
  analogWrite(LED1,counter);
  Serial.println(counter);

}
