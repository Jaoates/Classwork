
// this is the conversion factor from degrees to radians
const double deg2rad = 0.0174533;

// Definition of the new factorial function
long unsigned int myFactorial(unsigned int n){
  long unsigned int out = 1; 
  for (int i = n; i > 0; i = i-1){ //loop from n to 1 incrementing by one
    out = out*i; //multiply the previous value by n
  }
  return out;
}
// Defintion of the the new sin function
double mySin(double x, int n){ // x is the value we are taking the sin of, n is the number of terms used in the taylor series expansion
  double out = 0;
  for(int i = 0; i < n; i++){ // loop through all terms in the taylor expansion
    out = out + pow(-1,i)*pow(x,(2*i+1))/myFactorial(2*i+1); // add each term to the previous value
  }
  return out;
}

// Definition of a function to assist with testing and displaying
void testSin(double x,int n){
  Serial.print(x,4); //print the input value to 4 decimal points
  Serial.print(" | ");
  Serial.print(mySin(x,n),4); // print the sin of the input value with n terms in sin to 4 decimal points
  Serial.print(" | ");
  Serial.print(sin(x),4); //same as above but using builtin sin
  Serial.println("");
  delay(100);
}

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  delay(100);
  Serial.println("");

  // test myFactorial()
  // for(int i=0;i<=20;i++){
  //   Serial.println(i);
  //   Serial.println(myFactorial(i));
  // }

  int n = 1; // run test for 1 term
  Serial.print("n = ");
  Serial.println(n);
  double angle = 30*deg2rad; // run test for 30 deg
  testSin(angle,n);
  angle = 60*deg2rad; // run test for 60 deg
  testSin(angle,n);
  angle = 90*deg2rad; // run test for 90 deg
  testSin(angle,n);

  n = 3; //repeat for 3 terms
  Serial.print("n = ");
  Serial.println(n);
  angle = 30*deg2rad;
  testSin(angle,n);
  angle = 60*deg2rad;
  testSin(angle,n);
  angle = 90*deg2rad;
  testSin(angle,n);

  n = 5; //repeat for 5 terms
  Serial.print("n = ");
  Serial.println(n);
  angle = 30*deg2rad;
  testSin(angle,n);
  angle = 60*deg2rad;
  testSin(angle,n);
  angle = 90*deg2rad;
  testSin(angle,n);

}

void loop() { //dont have anything to do in loop
  // put your main code here, to run repeatedly:
  

}
