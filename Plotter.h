#ifdef __cplusplus
extern "C" {
#endif

#ifndef PLOTTER_H
#define PLOTTER_H

#include <stdint.h>
#include "RGB565Image.h"
#include "Config.h"

// Plotter Class
typedef struct RGB565Plotter {
  RGB565Image *img;
  uint16_t      working_x, working_y, working_w, working_h,
                innerOffset,
                plotArea_x, plotArea_y, plotArea_w, plotArea_h;

  void          (*Plot_Single_2D)           (struct RGB565Plotter *self, double *x, double *fx, uint32_t n, uint16_t color);
  void          (*Plot_Multiple_2D)         (struct RGB565Plotter *self, double *x, double **fx, uint32_t n, uint32_t nfx, uint16_t *colors);
  void          (*Add_Title)                (struct RGB565Plotter *self, char *label, font_t fontSize);
  void          (*Add_Label_X)              (struct RGB565Plotter *self, char *label);
  void          (*Add_Label_Y)              (struct RGB565Plotter *self, char *label);
} RGB565Plotter;

RGB565Plotter *RGB565Plotter_Init(RGB565Image *img, uint16_t x, uint16_t y, uint16_t width, uint16_t height);
void RGB565Plotter_Delete(RGB565Plotter *self);

#endif

#ifdef __cplusplus
}
#endif
