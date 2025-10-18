#include "ApogeePredictor.h"

ApogeePredictor::ApogeePredictor() : m_write_idx(0), m_num_samples(0), m_a_coeff(0.0f), m_b_coeff(0.0f), m_c_coeff(0.0f),
                                     m_coeffs_valid(false), m_sum_t(0.0f), m_sum_t2(0.0f), m_sum_t3(0.0f), m_sum_t4(0.0f),
                                     m_sum_h(0.0f), m_sum_th(0.0f), m_sum_t2h(0.0f) {}

void ApogeePredictor::addDataPoint(float time, float altitude) {
  // if buffer is full
  if (m_num_samples == APOGEE_BUFFER_SIZE) {
    uint16_t oldest_idx = (m_write_idx + 1) % APOGEE_BUFFER_SIZE;
    float old_t = m_time_buffer[oldest_idx];
    float old_h = m_height_buffer[oldest_idx];

    // remove old value from running sums
    float old_t2 = old_t * old_t;
    float old_t3 = old_t2 * old_t;
    float old_t4 = old_t2 * old_t2;

    m_sum_t -= old_t;
    m_sum_t2 -= old_t2;
    m_sum_t3 -= old_t3;
    m_sum_t4 -= old_t4;
    m_sum_h -= old_h;
    m_sum_th -= old_t * old_h;
    m_sum_t2h -= old_t2 * old_h;
  }

  // insert new data point
  m_time_buffer[m_write_idx] = time;
  m_height_buffer[m_write_idx] = altitude;

  // Update running sums
  float t2 = time * time;
  float t3 = t2 * time;
  float t4 = t2 * t2;

  m_sum_t += time;
  m_sum_t2 += t2;
  m_sum_t3 += t3;
  m_sum_t4 += t4;
  m_sum_h += altitude;
  m_sum_th += time * altitude;
  m_sum_t2h += t2 * altitude;

  // Update indices and count
  m_write_idx = (m_write_idx + 1) % APOGEE_BUFFER_SIZE;
  if (m_num_samples < APOGEE_BUFFER_SIZE) {
    m_num_samples++;
  }

  // Invalidate coefficients - need recalculation
  m_coeffs_valid = false;
}

bool ApogeePredictor::calculateCoefficients() {
  // require at least 3 sample for regression model
  if (m_num_samples < 3) {
    return false;
  }

  auto n = static_cast<float>(m_num_samples);

  // Build the normal equations matrix for quadratic regression
  // [sum_t4  sum_t3  sum_t2] [a]   [sum_t2h]
  // [sum_t3  sum_t2  sum_t ] [b] = [sum_th ]
  // [sum_t2  sum_t   n     ] [c]   [sum_h  ]

  // Using Cramer's rule for 3x3 system (no dynamic allocation)
  float det_A = n * (m_sum_t2 * m_sum_t4 - m_sum_t3 * m_sum_t3)
                - m_sum_t * (m_sum_t * m_sum_t4 - m_sum_t2 * m_sum_t3)
                + m_sum_t2 * (m_sum_t * m_sum_t3 - m_sum_t2 * m_sum_t2);

  if (MY_FABS(det_A) < 1e-10) {
    return false; // Matrix is singular
  }

  // Calculate a coefficient
  float det_a = m_sum_t2h * (m_sum_t2 * n - m_sum_t * m_sum_t)
                - m_sum_th * (m_sum_t3 * n - m_sum_t * m_sum_t2)
                + m_sum_h * (m_sum_t3 * m_sum_t - m_sum_t2 * m_sum_t2);

  // Calculate b coefficient
  float det_b = m_sum_t4 * (m_sum_th * n - m_sum_h * m_sum_t)
                - m_sum_t3 * (m_sum_t2h * n - m_sum_h * m_sum_t2)
                + m_sum_t2 * (m_sum_t2h * m_sum_t - m_sum_th * m_sum_t2);

  // Calculate c coefficient
  float det_c = m_sum_t4 * (m_sum_t2 * m_sum_h - m_sum_t * m_sum_th)
                - m_sum_t3 * (m_sum_t3 * m_sum_h - m_sum_t * m_sum_t2h)
                + m_sum_t2 * (m_sum_t3 * m_sum_th - m_sum_t2 * m_sum_t2h);

  m_a_coeff = det_a / det_A;
  m_b_coeff = det_b / det_A;
  m_c_coeff = det_c / det_A;

  m_coeffs_valid = true;
  return true;
}

ApogeePredictor::ApogeeResult_s ApogeePredictor::predictApogee() {
  ApogeeResult_s result = {0, 0, false};

  // Recalculate coefficients if needed
  if (!m_coeffs_valid) {
    if (!calculateCoefficients()) {
      return result;
    }
  }

  // For parabola h = at^2 + bt + c
  // Apogee occurs at t = -b/(2a) if a < 0
  if (m_a_coeff >= 0) {
    // Not a downward parabola - rocket still accelerating
    return result;
  }

  result.time = -m_b_coeff / (2.0f * m_a_coeff);

  // Check if predicted time is reasonable (not in the past)
  if (m_num_samples > 0) {
    float latest_time = m_time_buffer[(m_write_idx + APOGEE_BUFFER_SIZE - 1) % APOGEE_BUFFER_SIZE];
    if (result.time < latest_time) {
      // Apogee already passed or bad fit
      return result;
    }
  }

  // Calculate altitude at apogee
  result.altitude = m_a_coeff * result.time * result.time +
                    m_b_coeff * result.time + m_c_coeff;
  result.valid = true;

  return result;
}

float ApogeePredictor::getAltitudeAt(float time) {
  if (!m_coeffs_valid) {
    if (!calculateCoefficients()) {
      return 0;
    }
  }
  return m_a_coeff * time * time + m_b_coeff * time + m_c_coeff;
}

float ApogeePredictor::getFitQuality() {
  if (!m_coeffs_valid) {
    if (!calculateCoefficients()) {
      return 0;
    }
  }

  float mean_h = m_sum_h / m_num_samples;
  float ss_tot = 0;
  float ss_res = 0;

  // Calculate sum of squares
  uint16_t start_idx = (m_num_samples < APOGEE_BUFFER_SIZE) ? 0 : m_write_idx;
  for (uint16_t i = 0; i < m_num_samples; i++) {
    uint16_t idx = (start_idx + i) % APOGEE_BUFFER_SIZE;
    float h_actual = m_height_buffer[idx];
    float h_pred = getAltitudeAt(m_time_buffer[idx]);

    ss_tot += (h_actual - mean_h) * (h_actual - mean_h);
    ss_res += (h_actual - h_pred) * (h_actual - h_pred);
  }

  if (ss_tot < 1e-10) return 0;
  return 1.0f - (ss_res / ss_tot);
}

void ApogeePredictor::reset() {
  m_write_idx = 0;
  m_num_samples = 0;
  m_coeffs_valid = false;
  m_sum_t = m_sum_t2 = m_sum_t3 = m_sum_t4 = 0;
  m_sum_h = m_sum_th = m_sum_t2h = 0;
}

int ApogeePredictor::getDataCount() const {
  return m_num_samples;
}
