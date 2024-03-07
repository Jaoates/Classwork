int n = 0; // set up some variables we will use later
float a =0;
float b =1.1;
double result;
int stopVal = 20;

// init serial in the setup funciton
void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
}


void loop() {
  // put your main code here, to run repeatedly:
  result = pow((a+b),n); // compute (a+b)^n
  Serial.println(result);// print this
  Serial.println(n); // print n
  n++; // n = n+1
  if (n>stopVal){ // check if we should kill the program
    Serial.println("done"); 
    delay(100); //needed so that the serial library can complete the transaction
    exit(0);
  }
}
