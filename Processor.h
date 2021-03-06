#ifdef __cplusplus
extern "C" {
#endif

#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <stdint.h>
#include "RGB565Image.h"
#include "Config.h"

// Processor Class
typedef struct RGB565Processor {
  RGB565Image *img;

  // Scale image by scale factor > 0
  void          (*Scale)              (struct RGB565Processor *self, float scale, scale_t mode);

  // Rotate image by degree value (deg)
  void          (*Rotate)             (struct RGB565Processor *self, float deg, scale_t mode);

  // Grayscale an image
  void          (*Grayscale)          (struct RGB565Processor *self);

  // Insert a second image at top left (x, y) coordinates using op
  void          (*Insert)             (struct RGB565Processor *self, RGB565Image *second, uint16_t x, uint16_t y, img_op_t op);

  // Simple RGB color inversion
  void          (*Invert)             (struct RGB565Processor *self);

  // Histogram manipulation by power function with exponent (a) applied to colorComponent
  void          (*Point_Curve_Pow)    (struct RGB565Processor *self, float a, uint8_t colorComponent);

  // Histogram manipulation with Lightroom-like point curve, given n points (x, fx), applied to colorComponent (cubic spline interpolated)
  void          (*Point_Curve_Points) (struct RGB565Processor *self, double *x, double *fx, uint8_t n, uint8_t colorComponent);

  // Blur Image. Add a dreamy magical (gaussian) blur. Dreaminess = drunkness (how blurry!)
  void          (*Dreamify)           (struct RGB565Processor *self, uint8_t dreaminess);
} RGB565Processor;

RGB565Processor *RGB565Processor_Init(RGB565Image *img);
void RGB565Processor_Delete(RGB565Processor *self);

#endif

#ifdef __cplusplus
}
#endif
