#include "Nicla_System.h"

volatile unsigned long pulseStart = 0;
volatile int rcPulse = 1500;  // Default neutral
int ref = 1500;
const int RC_PIN = 1;
const int MOTORPIN = 2;
const int DEADZONE = 50;

void setup() {
  pinMode(RC_PIN, INPUT);
  pinMode(MOTORPIN, OUTPUT);
  attachInterrupt(digitalPinToInterrupt(RC_PIN), rcISR, CHANGE);
  Serial.begin(115200);
  analogWrite(MOTORPIN,0);
}

void loop() {
  Serial.println(digitalRead(RC_PIN));
  int pwmValue;
  if (abs(rcPulse - ref) < DEADZONE) {
    pwmValue = 0;  // Off at neutral
  } else if (rcPulse > ref) {
    pwmValue = map(rcPulse, ref+50, 2000, 1, 255);  // Forward speed
  } else {
    pwmValue = 0;  // Or map to reverse if bidirectional needed (requires H-bridge)
  }
  pwmValue = constrain(pwmValue, 0, 255);
  analogWrite(MOTORPIN, 0);
  delay(20);  // ~50Hz loop
}

void rcISR() {
  if (digitalRead(RC_PIN) == HIGH) {
    pulseStart = micros();
  } else {
    rcPulse = micros() - pulseStart;
  }
}