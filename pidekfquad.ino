#include <Wire.h>
#include <sbus.h>
#include <math.h>
#include <string.h>

#define SBUS_BAUD_RATE 100000
bfs::SbusRx sbus_rx(&Serial5);
#define MICOLINK_MSG_HEAD 0xEF
#define MICOLINK_MAX_PAYLOAD_LEN 64
#define MICOLINK_MAX_LEN (MICOLINK_MAX_PAYLOAD_LEN + 7)
enum {
  MICOLINK_MSG_ID_RANGE_SENSOR = 0x51,
};
typedef struct {
  uint8_t head;
  uint8_t dev_id;
  uint8_t sys_id;
  uint8_t msg_id;
  uint8_t seq;
  uint8_t len;
  uint8_t payload[MICOLINK_MAX_PAYLOAD_LEN];
  uint8_t checksum;
  uint8_t status;
  uint8_t payload_cnt;
} MICOLINK_MSG_t;
#pragma pack(push, 1)
typedef struct {
  uint32_t time_ms;
  uint32_t distance;
  uint8_t  strength;
  uint8_t  precision;
  uint8_t  dis_status;
  uint8_t  reserved1;
  int16_t  flow_vel_x;
  int16_t  flow_vel_y;
  uint8_t  flow_quality;
  uint8_t  flow_status;
  uint16_t reserved2;
} MICOLINK_PAYLOAD_RANGE_SENSOR_t;
#pragma pack(pop)

static const int NX = 12; 
static const int NY = 13;  
static const double MASS   = 1.70;  
static const double GRAV   = 9.81;   
static const double PWM2T = 0.0215;
double ch0 = 1001, ch1 = 1000, ch2 = 172, ch3 = 988, ch4 = 1294, ch5 = 1811;
double phird   = 0, thetard = 0, psird = 0, zrd = 0;
double phir    = 0, thetar  = 0;
double z = 0;
double xd = 0, yd = 0, zd = 0;
double xr = 0, yr = 0, zr = 0;
double izd = 0, ivphi = 0, ivtheta = 0, ivpsi = 0;
double pre_errorzd = 0, pre_errorvphi = 0, pre_errorvtheta = 0, pre_errorvpsi = 0;
double pwm_value1 = 819, pwm_value2 = 819, pwm_value3 = 819, pwm_value4 = 819;
double U1pid = 0, U2pid = 0, U3pid = 0, U4pid = 0;
double vphi = 0, vtheta = 0, vpsi = 0;
double vphicalibration = 0, vthetacalibration = 0, vpsicalibration = 0;
int    RateCalibrationNumber;
double AccX = 0, AccY = 0, AccZ = 0;
double phi   = 0, theta = 0, psi = 0;
uint32_t LoopTimer;
double dt = 0.005;

double phikalman     = 0.0;
double phiunkalman   = 0.07;
double thetakalman   = 0.0;
double thetaunkalman = 0.07;
double dtOptical = 0.02;
unsigned long lastOpticalTime = 0;
double xzd = 0.0;
double yzd = 0.0;

const uint32_t RISING_TIME_MS = 3000;
const double   HOVER_BASE_US  = 1365.0;
bool           rising_active   = false;
uint32_t       rising_start_ms = 0;
uint32_t       rising_end_ms   = 0;
double         pwm_stable      = HOVER_BASE_US;
double         pwm_stable_start = 1000.0f;
double         pwm_stable_end   = HOVER_BASE_US;

inline void start_rising() {
  rising_active    = true;
  rising_start_ms  = millis();
  rising_end_ms    = rising_start_ms + RISING_TIME_MS;
  pwm_stable_start = 1000.0f;
  pwm_stable_end   = HOVER_BASE_US;
}

double x_hat[NX] = {0.0};    
double P[NX][NX];             
double Qm[NX][NX];           
double Rm[NY][NY];            
double Cmat[NY][NX];         
void ekfInit(double phi0, double theta0, double psi0, double z0);
void ekfStep(double dt_local, const double z_meas[NY]);
void ekfGetState(double *x_out);
bool invertMatrix13(const double A[NY][NY], double Ainv[NY][NY]);

void micolink_decode(uint8_t data);
bool micolink_parse_char(MICOLINK_MSG_t* msg, uint8_t data);
void read_receiver();
void gyro_signals(void);
void getoptical();

void pid_controller(double &Uout, double &integral, double &pre_error,
                    double error, double kp, double ki, double kd, double dtc) {
  integral += error * dtc;
  integral = constrain(integral, -1.5, 1.5);
  double derivative = (error - pre_error) / dtc;
  Uout = kp * error + ki * integral + kd * derivative;
  Uout = constrain(Uout, -150.0, 150.0);
  pre_error = error;
}

void kalman_angle(double &angleesti, double &angleun, double rate_gyro, double angle_meas) {
  const double dtKF = 0.005;             
  const double q    = 0.035 * 0.035;     
  const double r    = 0.00275;           
  angleesti = angleesti + dtKF * rate_gyro;
  angleun   = angleun   + dtKF * dtKF * q;
  double kgain = angleun / (angleun + r);
  angleesti = angleesti + kgain * (angle_meas - angleesti);
  angleun   = (1.0 - kgain) * angleun;
}

bool invertMatrix13(const double A[NY][NY], double Ainv[NY][NY])
{
  const int N = NY; 
  double aug[NY][2*NY];
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      aug[i][j]   = A[i][j];
      aug[i][j+N] = (i == j) ? 1.0 : 0.0;
    }
  }
  for (int col = 0; col < N; ++col) {
    int pivot = col;
    double max_val = fabs(aug[col][col]);
    for (int row = col+1; row < N; ++row) {
      double v = fabs(aug[row][col]);
      if (v > max_val) {
        max_val = v;
        pivot = row;
      }
    }
    if (max_val < 1e-12) return false;
    if (pivot != col) {
      for (int j = 0; j < 2*N; ++j) {
        double tmp = aug[col][j];
        aug[col][j] = aug[pivot][j];
        aug[pivot][j] = tmp;
      }
    }
    double diag = aug[col][col];
    for (int j = 0; j < 2*N; ++j) {
      aug[col][j] /= diag;
    }

    for (int row = 0; row < N; ++row) {
      if (row == col) continue;
      double factor = aug[row][col];
      for (int j = 0; j < 2*N; ++j) {
        aug[row][j] -= factor * aug[col][j];
      }
    }
  }
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      Ainv[i][j] = aug[i][j+N];

  return true;
}

void ekfInit(double phi0, double theta0, double psi0, double z0)
{
  for (int i = 0; i < NX; ++i) x_hat[i] = 0.0;
  double z_init = (z0 > 0.025) ? z0 : 0.025;
  x_hat[2] = z_init;  
  x_hat[3] = phi0;    
  x_hat[4] = theta0;   
  x_hat[5] = psi0;    
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NX; ++j) {
      P[i][j]  = (i == j) ? 0.01  : 0.0;
      Qm[i][j] = (i == j) ? 1e-4 : 0.0;  
    }
  }
  for (int i = 0; i < NY; ++i) {
    for (int j = 0; j < NY; ++j) {
      double val = 0.0;
      if (i == j) {
        if (i <= 3)        val = 1e-2;   // z, phi, theta, psi
        else if (i <= 9)   val = 1e-2;   // velocities & rates
        else               val = 1e-1;   // accelerations more noisy
      }
      Rm[i][j] = val;
    }
  }
  for (int i = 0; i < NY; ++i)
    for (int j = 0; j < NX; ++j)
      Cmat[i][j] = 0.0;
  Cmat[0][2] = 1.0;
  Cmat[1][3] = 1.0;
  Cmat[2][4] = 1.0;
  Cmat[3][5] = 1.0;
  Cmat[4][6] = 1.0;
  Cmat[5][7] = 1.0;
  Cmat[6][8] = 1.0;
  Cmat[7][9]  = 1.0;
  Cmat[8][10] = 1.0;
  Cmat[9][11] = 1.0;
}

void ekfStep(double dt_local, const double z_meas[NY])
{
  double A[NX][NX];
  for (int i = 0; i < NX; ++i)
    for (int j = 0; j < NX; ++j)
      A[i][j] = (i == j) ? 1.0 : 0.0;
  A[0][6]  = dt_local;
  A[1][7]  = dt_local; 
  A[2][8]  = dt_local; 
  A[3][9]  = dt_local; 
  A[4][10] = dt_local; 
  A[5][11] = dt_local; 
  double x_pred[NX];
  for (int i = 0; i < NX; ++i) {
    x_pred[i] = 0.0;
    for (int j = 0; j < NX; ++j) {
      x_pred[i] += A[i][j] * x_hat[j];
    }
  }
  double AP[NX][NX];
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NX; ++j) {
      AP[i][j] = 0.0;
      for (int k = 0; k < NX; ++k) {
        AP[i][j] += A[i][k] * P[k][j];
      }
    }
  }
  double P_pred[NX][NX];
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NX; ++j) {
      P_pred[i][j] = 0.0;
      for (int k = 0; k < NX; ++k) {
        P_pred[i][j] += AP[i][k] * A[j][k]; 
      }
      P_pred[i][j] += Qm[i][j];
    }
  }
  double y_pred[NY];
  for (int i = 0; i < NY; ++i) {
    y_pred[i] = 0.0;
    for (int j = 0; j < NX; ++j) {
      y_pred[i] += Cmat[i][j] * x_pred[j];
    }
  }
  double innov[NY];
  for (int i = 0; i < NY; ++i) {
    innov[i] = z_meas[i] - y_pred[i];
  }
  double CP[NY][NX];
  for (int i = 0; i < NY; ++i) {
    for (int j = 0; j < NX; ++j) {
      CP[i][j] = 0.0;
      for (int k = 0; k < NX; ++k) {
        CP[i][j] += Cmat[i][k] * P_pred[k][j];
      }
    }
  }
  double S[NY][NY];
  for (int i = 0; i < NY; ++i) {
    for (int j = 0; j < NY; ++j) {
      S[i][j] = 0.0;
      for (int k = 0; k < NX; ++k) {
        S[i][j] += CP[i][k] * Cmat[j][k];
      }
      S[i][j] += Rm[i][j];
    }
  }
  double S_inv[NY][NY];
  if (!invertMatrix13(S, S_inv)) {
    for (int i = 0; i < NX; ++i)
      x_hat[i] = x_pred[i];
    for (int i = 0; i < NX; ++i)
      for (int j = 0; j < NX; ++j)
        P[i][j] = P_pred[i][j];
    return;
  }
  double PCt[NX][NY];
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NY; ++j) {
      PCt[i][j] = 0.0;
      for (int k = 0; k < NX; ++k) {
        PCt[i][j] += P_pred[i][k] * Cmat[j][k];
      }
    }
  }
  double K[NX][NY];
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NY; ++j) {
      K[i][j] = 0.0;
      for (int k = 0; k < NY; ++k) {
        K[i][j] += PCt[i][k] * S_inv[k][j];
      }
    }
  }
  double x_new[NX];
  for (int i = 0; i < NX; ++i) {
    x_new[i] = x_pred[i];
    for (int j = 0; j < NY; ++j) {
      x_new[i] += K[i][j] * innov[j];
    }
  }
  double I_KC[NX][NX];
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NX; ++j) {
      double kc = 0.0;
      for (int k = 0; k < NY; ++k) {
        kc += K[i][k] * Cmat[k][j];
      }
      I_KC[i][j] = (i == j ? 1.0 : 0.0) - kc;
    }
  }
  double P_new[NX][NX];
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NX; ++j) {
      P_new[i][j] = 0.0;
      for (int k = 0; k < NX; ++k) {
        P_new[i][j] += I_KC[i][k] * P_pred[k][j];
      }
    }
  }
  for (int i = 0; i < NX; ++i)
    x_hat[i] = x_new[i];
  for (int i = 0; i < NX; ++i)
    for (int j = 0; j < NX; ++j)
      P[i][j] = P_new[i][j];
}

void ekfGetState(double *x_out)
{
  for (int i = 0; i < NX; ++i)
    x_out[i] = x_hat[i];
}

void getoptical() {
  while (Serial3.available()) {
    uint8_t sensor_data = Serial3.read();
    micolink_decode(sensor_data);
  }
}

void micolink_decode(uint8_t data) {
  static MICOLINK_MSG_t msg = {};   
  if (!micolink_parse_char(&msg, data)) return;

  if (msg.msg_id == MICOLINK_MSG_ID_RANGE_SENSOR) {
    MICOLINK_PAYLOAD_RANGE_SENSOR_t payload;
    memcpy(&payload, msg.payload, msg.len);
    dtOptical = (millis() - lastOpticalTime) / 1000.0f;
    lastOpticalTime = millis();
    if (payload.distance < 8000 &&
        fabs(payload.flow_vel_x) < 1500 &&
        fabs(payload.flow_vel_y) < 1500) {
      if (dtOptical > 0.0) {
        double range_m = payload.distance / 1000.0;
        double z_meas = range_m * cos(phikalman) * cos(thetakalman);
        zd = (z_meas - z) / dtOptical;
        z  = z_meas;
        if (fabs(zd) > 4.0) zd = 0.0;
        xzd = payload.flow_vel_x * 0.01;   
        yzd = -payload.flow_vel_y * 0.01;  
      }
    }
  }
}
bool micolink_parse_char(MICOLINK_MSG_t* msg, uint8_t data) {
  switch (msg->status) {
    case 0:
      if (data == MICOLINK_MSG_HEAD) {
        msg->head   = data;
        msg->status = 1;
      }
      break;
    case 1: msg->dev_id = data; msg->status = 2; break;
    case 2: msg->sys_id = data; msg->status = 3; break;
    case 3: msg->msg_id = data; msg->status = 4; break;
    case 4: msg->seq    = data; msg->status = 5; break;
    case 5:
      msg->len = data;
      if (msg->len == 0) {
        msg->status = 7;
      } else if (msg->len > MICOLINK_MAX_PAYLOAD_LEN) {
        msg->status      = 0;
        msg->payload_cnt = 0;
      } else {
        msg->payload_cnt = 0;
        msg->status      = 6;
      }
      break;
    case 6:
      msg->payload[msg->payload_cnt++] = data;
      if (msg->payload_cnt == msg->len) {
        msg->status = 7;
      }
      break;
    case 7:
      msg->checksum = data;
      msg->status   = 0;
      return true;
    default:
      msg->status      = 0;
      msg->payload_cnt = 0;
      break;
  }
  return false;
}

void read_receiver() {
  if (sbus_rx.Read()) {
    bfs::SbusData data = sbus_rx.data();
    ch0 = data.ch[0];
    ch1 = data.ch[1];
    ch2 = data.ch[2];
    ch3 = data.ch[3];
    ch4 = data.ch[4];
    ch5 = data.ch[6];
  }
}

void gyro_signals(void) {
  Wire.beginTransmission(0x68);
  Wire.write(0x1A);
  Wire.write(0x05);
  Wire.endTransmission();
  Wire.beginTransmission(0x68);
  Wire.write(0x1C);
  Wire.write(0x10);
  Wire.endTransmission();
  Wire.beginTransmission(0x68);
  Wire.write(0x3B);
  Wire.endTransmission();
  Wire.requestFrom(0x68, 6);
  int16_t AccXLSB = Wire.read() << 8 | Wire.read();
  int16_t AccYLSB = Wire.read() << 8 | Wire.read();
  int16_t AccZLSB = Wire.read() << 8 | Wire.read();
  Wire.beginTransmission(0x68);
  Wire.write(0x1B);
  Wire.write(0x08);
  Wire.endTransmission();
  Wire.beginTransmission(0x68);
  Wire.write(0x43);
  Wire.endTransmission();
  Wire.requestFrom(0x68, 6);
  int16_t GyroX = Wire.read() << 8 | Wire.read();
  int16_t GyroY = Wire.read() << 8 | Wire.read();
  int16_t GyroZ = Wire.read() << 8 | Wire.read();
  vphi   = (double)GyroX / 65.5 * PI / 180.0;
  vtheta = (double)GyroY / 65.5 * PI / 180.0;
  vpsi   = (double)GyroZ / 65.5 * PI / 180.0;
  AccX = (double)AccXLSB / 4096 - 0.039;
  AccY = (double)AccYLSB / 4096 + 0.001;
  AccZ = (double)AccZLSB / 4096 - 0.214;
  phi   = atan2(AccY, sqrt(AccX * AccX + AccZ * AccZ));
  theta = -atan2(AccX, sqrt(AccY * AccY + AccZ * AccZ));
}

void setup() {
  Serial.begin(115200);
  Wire.setClock(400000);
  Wire.begin();
  delay(250);
  Wire.beginTransmission(0x68);
  Wire.write(0x6B);
  Wire.write(0x00);
  Wire.endTransmission();
  for (RateCalibrationNumber = 0; RateCalibrationNumber < 2000; RateCalibrationNumber++) {
    gyro_signals();
    vphicalibration   += vphi;
    vthetacalibration += vtheta;
    vpsicalibration   += vpsi;
    delay(1);
  }
  vphicalibration   /= 2000;
  vthetacalibration /= 2000;
  vpsicalibration   /= 2000;
  Serial5.begin(SBUS_BAUD_RATE, SERIAL_8E2);
  sbus_rx.Begin();
  Serial3.begin(115200);
  analogWriteFrequency(3, 200);
  analogWriteFrequency(4, 200);
  analogWriteFrequency(5, 200);
  analogWriteFrequency(6, 200);
  analogWriteResolution(12);
  analogWrite(3, 819);
  analogWrite(4, 819);
  analogWrite(5, 819);
  analogWrite(6, 819);
  delay(5000);
  gyro_signals();
  double phi0   = phi;
  double theta0 = theta;
  double psi0   = 0.0;
  double z0     = 0.05;   
  phikalman   = phi0;
  thetakalman = theta0;
  ekfInit(phi0, theta0, psi0, z0);
  LoopTimer       = micros();
  lastOpticalTime = millis();
}

void loop() {
  read_receiver();
  gyro_signals();
  getoptical();
  vphi   -= vphicalibration;
  vtheta -= vthetacalibration;
  vpsi   -= vpsicalibration;
  kalman_angle(phikalman,   phiunkalman,   vphi,   phi);
  kalman_angle(thetakalman, thetaunkalman, vtheta, theta);
  double accx = AccX * 9.8065;
  double accy = AccY * 9.8065;
  double accz = (AccZ - 1.0) * 9.8065;
  double vx_meas = xzd * z;   
  double vy_meas = yzd * z;  

  xd  = vx_meas;
  yd  = vy_meas;
  double x_est[12];
  ekfGetState(x_est);
  double xekf      = x_est[0];
  double yekf      = x_est[1];
  double zekf      = x_est[2];
  double phiekf    = x_est[3];
  double thetaekf  = x_est[4];
  double psiekf    = x_est[5];
  double xdekf     = x_est[6];
  double ydekf     = x_est[7];
  double zdekf     = x_est[8];
  double phidekf   = x_est[9];
  double thetadekf = x_est[10];
  double psidekf   = x_est[11];

 
  zr = map(ch2, 180, 1811, 0.0, 2.0);
  xr = xr + (1001 - ch1) / 1200.0 * 0.01;
  yr = yr + (1000 - ch0) / 1200.0 * 0.01;
  zrd    = 0.5 * (zr - zekf) - 0.1 * zdekf;
  thetar = 0.1 * exp(-0.25 * zekf) * (xr - xekf) - 0.2 * exp(-0.5 * zekf) * xdekf;
  phir   = 0.1 * exp(-0.25 * zekf) * (yekf - yr) + 0.2 * exp(-0.5 * zekf) * ydekf;
  psird  = 4e-4 * (988 - ch3);
  phir   = constrain(phir,   -0.3, 0.3);
  thetar = constrain(thetar, -0.3, 0.3);

  phird   = 4.75 * (phir   - phikalman)   - 0.075 * vphi;
  thetard = 4.75 * (thetar - thetakalman) - 0.075 * vtheta;
  phird   = constrain(phird,   -0.75, 0.75);
  thetard = constrain(thetard, -0.75, 0.75);
  psird   = constrain(psird,   -0.75, 0.75);

  double errorzd     = zrd     - zdekf;
  double errorvphi   = phird   - vphi;
  double errorvtheta = thetard - vtheta;
  double errorvpsi   = psird   - vpsi;

  pid_controller(U1pid, izd,     pre_errorzd,     errorzd,     62.75, 122.25, 3.75, dt);
  pid_controller(U2pid, ivphi,   pre_errorvphi,   errorvphi,   43.5,  162.5,  2.25, dt);
  pid_controller(U3pid, ivtheta, pre_errorvtheta, errorvtheta, 43.5,  162.5,  2.25, dt);
  pid_controller(U4pid, ivpsi,   pre_errorvpsi,   errorvpsi,   42.5,  175.5,  1.7,  dt);

  static bool was_auto = false;
  bool is_auto   = (ch4 < 700);
  bool is_manual = (ch4 >= 700 && ch4 < 1050);

  if (is_auto && !was_auto) {
    start_rising();
  }
  was_auto = is_auto;

  if (ch4 < 700) {
    if (rising_active) {
      unsigned long now_ms = millis();
      double alpha = (double)(now_ms - rising_start_ms) / (double)RISING_TIME_MS;
      if (alpha < 0.0f) alpha = 0.0f;
      if (alpha > 1.0f) alpha = 1.0f;
      pwm_stable = pwm_stable_start + alpha * (pwm_stable_end - pwm_stable_start);
      if ((int32_t)(now_ms - rising_end_ms) >= 0) rising_active = false;
    } else {
      pwm_stable = HOVER_BASE_US;
    }

    double m1_us = pwm_stable + U1pid - U2pid - U3pid - U4pid;
    double m2_us = pwm_stable + U1pid + U2pid - U3pid + U4pid;
    double m3_us = pwm_stable + U1pid - U2pid + U3pid + U4pid;
    double m4_us = pwm_stable + U1pid + U2pid + U3pid - U4pid;

    pwm_value1 = 0.819 * m1_us;
    pwm_value2 = 0.819 * m2_us;
    pwm_value3 = 0.819 * m3_us;
    pwm_value4 = 0.819 * m4_us;

    pwm_value1 = constrain(pwm_value1, 819, 1600);
    pwm_value2 = constrain(pwm_value2, 819, 1600);
    pwm_value3 = constrain(pwm_value3, 819, 1600);
    pwm_value4 = constrain(pwm_value4, 819, 1600);
  }
  else if (is_manual) {
    double U1manual = map(ch2, 172, 1804, 1000, 1700);

    double m1_us = U1manual - U2pid - U3pid - U4pid;
    double m2_us = U1manual + U2pid - U3pid + U4pid;
    double m3_us = U1manual - U2pid + U3pid + U4pid;
    double m4_us = U1manual + U2pid + U3pid - U4pid;

    pwm_value1 = 0.819 * m1_us;
    pwm_value2 = 0.819 * m2_us;
    pwm_value3 = 0.819 * m3_us;
    pwm_value4 = 0.819 * m4_us;

    pwm_value1 = constrain(pwm_value1, 819, 1600);
    pwm_value2 = constrain(pwm_value2, 819, 1600);
    pwm_value3 = constrain(pwm_value3, 819, 1600);
    pwm_value4 = constrain(pwm_value4, 819, 1600);

    izd = 0;
  }
  else {
    pwm_value1 = pwm_value2 = pwm_value3 = pwm_value4 = 819;
    izd        = ivphi = ivtheta = ivpsi = 0.0;
    xr         = xekf;
    yr         = yekf;
    zr         = zekf;
  }

  double z_meas[NY];
  z_meas[0]  = z;            
  z_meas[1]  = phikalman;    
  z_meas[2]  = thetakalman; 
  z_meas[3]  = psi;          
  z_meas[4]  = xd;          
  z_meas[5]  = yd;          
  z_meas[6]  = zd;          
  z_meas[7]  = vphi;         
  z_meas[8]  = vtheta;      
  z_meas[9]  = vpsi;         
  z_meas[10] = accx;       
  z_meas[11] = accy;        
  z_meas[12] = accz;        
  ekfStep(dt, z_meas);
  ekfGetState(x_est);
  xekf      = x_est[0];
  yekf      = x_est[1];
  zekf      = x_est[2];
  phiekf    = x_est[3];
  thetaekf  = x_est[4];
  psiekf    = x_est[5];
  xdekf     = x_est[6];
  ydekf     = x_est[7];
  zdekf     = x_est[8];
  phidekf   = x_est[9];
  thetadekf = x_est[10];
  psidekf   = x_est[11];

  Serial.print(xr); Serial.print(",");
  Serial.print(yr); Serial.print(",");
  Serial.print(zr); Serial.print(",");
  Serial.print(phir); Serial.print(",");
  Serial.print(thetar); Serial.print(",");
  Serial.print(xekf); Serial.print(",");
  Serial.print(yekf); Serial.print(",");
  Serial.print(zekf); Serial.print(",");
  Serial.print(phiekf); Serial.print(",");
  Serial.print(thetaekf); Serial.print(",");
  Serial.print(psiekf); Serial.print(",");
  Serial.print(xdekf); Serial.print(",");
  Serial.print(ydekf); Serial.print(",");
  Serial.print(zdekf); Serial.print(",");
  Serial.print(phidekf); Serial.print(",");
  Serial.print(thetadekf); Serial.print(",");
  Serial.print(psidekf); Serial.print(",");
  Serial.print(accx); Serial.print(",");
  Serial.print(accy); Serial.print(",");
  Serial.print(accz); Serial.print(",");
  Serial.print(U1pid); Serial.print(",");
  Serial.print(U2pid); Serial.print(",");
  Serial.print(U3pid); Serial.print(",");
  Serial.print(U4pid); Serial.println(",");
  while (micros() - LoopTimer < 5000) {}
  dt = (micros() - LoopTimer) / 1000000.0;
  LoopTimer = micros();
  psi = psi + dt * vpsi;
  analogWrite(3, (int)pwm_value3);
  analogWrite(4, (int)pwm_value1);
  analogWrite(5, (int)pwm_value4);
  analogWrite(6, (int)pwm_value2);
}
