//
// Created by chris on 9/8/2025.
//

#ifndef PLATFORMIO_NULI_AIRBRAKES_SOFTWARE_APOGEEPREDICTOR_H
#define PLATFORMIO_NULI_AIRBRAKES_SOFTWARE_APOGEEPREDICTOR_H

#include <cstdint>

// Helpers
#define MY_FABS(x) ((x) < 0.0 ? -(x) : (x))
//inline double my_fabs(double x) {
//  return x < 0.0 ? -x : x;
//}


#define APOGEE_BUFFER_SIZE 100

class ApogeePredictor {
public:
  ApogeePredictor();

  void addDataPoint(float time, float altitude);

  bool calculateCoefficients();

  struct ApogeeResult_s {
      float time;
      float altitude;
      bool valid;
  };

  ApogeeResult_s predictApogee();

  float getAltitudeAt(float time);

  float getFitQuality();

  void reset();

  int getDataCount() const;
protected:
private:
  float m_time_buffer[APOGEE_BUFFER_SIZE] = {0};
  float m_height_buffer[APOGEE_BUFFER_SIZE] = {0};

  uint16_t m_write_idx = 0;
  uint16_t m_num_samples = 0;

  // Cached regression coefficients (h = at^2 + bt + c)
  float m_a_coeff;
  float m_b_coeff;
  float m_c_coeff;
  bool m_coeffs_valid;

  // Running sums for efficient regression calculation
  float m_sum_t;
  float m_sum_t2;
  float m_sum_t3;
  float m_sum_t4;
  float m_sum_h;
  float m_sum_th;
  float m_sum_t2h;

};


#endif //PLATFORMIO_NULI_AIRBRAKES_SOFTWARE_APOGEEPREDICTOR_H
