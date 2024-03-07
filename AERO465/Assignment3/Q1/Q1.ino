#include <NewPing.h>

int TRIGGER_PIN = 9;
int ECHO_PIN = 10;
int MAX_DISTANCE = 400;

NewPing sonar(TRIGGER_PIN, ECHO_PIN, MAX_DISTANCE); // NewPing setup of pins and maximum distance.

void setup() {
  Serial.begin(9600);
}

void loop() {
  delay(50);                    // Wait 50ms between pings (about 20 pings/sec). 29ms should be the shortest delay between pings.

  int distance = sonar.ping_cm(); // Send ping, get distance in cm and print result (0 = outside set distance range)

  // Serial.print("Distance: ");
  Serial.println(distance);
  // Serial.println("cm");
}



