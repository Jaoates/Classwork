#include <iostream>
#include <vector>
#include <string>

using namespace std;

int myFactorial(int n){
  int out = 1;
  // Serial.println(n);
  for (int i = n; i != 0; i = i-1){
    out = out*i;
    // Serial.println(i);
  }
  return out;
}

double mySin(double x){
    double out = 0;
    for(int i = 0; i<9; i++){
        // Serial.println(i);
        out = out + pow(-1,i)*(pow(x,(2*i+1)))/myFactorial(2*i+1);
        // Serial.println(pow(-1,i)*(pow(x,(2*i+1)))/myFactorial(2*i+1));
    }
    return out;
}

int main(){
    for(double i = 0; i<7;i = i+.01){
        cout << i;
        cout << " | ";
        cout << sin(i);
        cout << '\n';
    }
}