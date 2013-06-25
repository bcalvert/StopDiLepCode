#include <cmath>
float nDigits(float number, int digits) {
  return round(number * std::pow(10.,digits)) / std::pow(10.,digits);
}
