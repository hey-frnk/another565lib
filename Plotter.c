#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "DrawHandler.h"
#include "Plotter.h"
#include "Processor.h"

static inline double _d_yminmax(double *arr, uint32_t n, bool max) {
  if(max) {
    double max = arr[0];
    for(uint32_t i = 0; i < n; ++i) if(max < arr[i]) max = arr[i];
    return max;
  } else {
    double min = arr[0];
    for(uint32_t i = 0; i < n; ++i) if(min > arr[i]) min = arr[i];
    return min;
  }
}

// Internal Frame & Grid Drawing
static void __Plot_2D_Create_Axes(struct RGB565Plotter *self) {
  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(self->img);

  for(uint16_t i = 1; i < PLOT_DIVISIONS_X; ++i) _h->DrawLine(_h,
      self->plotArea_x + i * (self->plotArea_w / PLOT_DIVISIONS_X),
      self->plotArea_y + self->plotArea_h - PLOT_MARK_LENGTH - 1,
      self->plotArea_x + i * (self->plotArea_w / PLOT_DIVISIONS_X),
      self->plotArea_y + self->plotArea_h - 1, 0
    );

  for(uint16_t i = 1; i < PLOT_DIVISIONS_Y; ++i) _h->DrawLine(_h,
      self->plotArea_x,
      self->plotArea_y + i * (self->plotArea_h / PLOT_DIVISIONS_Y),
      self->plotArea_x + PLOT_MARK_LENGTH,
      self->plotArea_y + i * (self->plotArea_h / PLOT_DIVISIONS_Y), 0
    );

  // Draw Grid
  for(uint16_t i = 1; i < PLOT_DIVISIONS_X; ++i) _h->DrawVLine(_h,
      self->plotArea_x + i * (self->plotArea_w / PLOT_DIVISIONS_X), // x pos
      self->plotArea_y + 1, // y pos
      self->plotArea_h - PLOT_MARK_LENGTH - 1, 0xC618, LINE_STYLE_DOT
    );
  for(uint16_t i = 1; i < PLOT_DIVISIONS_Y; ++i) _h->DrawHLine(_h,
      self->plotArea_x + PLOT_MARK_LENGTH,
      self->plotArea_y + i * (self->plotArea_h / PLOT_DIVISIONS_Y),
      self->plotArea_w - PLOT_MARK_LENGTH - 1, 0xC618, LINE_STYLE_DOT
    );
  RGB565ImageDrawHandler_Delete(_h);
}

static void __Plot_2D_Create_Scale(struct RGB565Plotter *self, double xMin, double xMax, double yMin, double yMax) {
  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(self->img);

  double xAbs = fabs(xMin) + fabs(xMax), xDiv = xAbs / (double)PLOT_DIVISIONS_X,
         yAbs = fabs(yMin) + fabs(yMax), yDiv = yAbs / (double)PLOT_DIVISIONS_Y;

  for(uint16_t i = 0; i < (PLOT_DIVISIONS_X + 1); ++i) {
    float _xLinSpV = xMin + ((float)i * xDiv);
    uint32_t _prStrSz = snprintf(NULL, 0, "%5.2f", _xLinSpV);
    char *_xLinSpS = (char *)calloc(_prStrSz + 1, sizeof(char));
    sprintf(_xLinSpS, "%5.2f", _xLinSpV);
    _h->DrawString(_h,
      self->plotArea_x + i * (self->plotArea_w / PLOT_DIVISIONS_X),
      self->plotArea_y + self->plotArea_h + PLOT_VALUES_OFFSET,
      _xLinSpS, 0, FONT_5X7, ALIGNMENT_CENTER);
    free(_xLinSpS);
  }
  for(uint16_t i = 0; i < (PLOT_DIVISIONS_Y + 1); ++i) {
    float _yLinSpV = yMax - ((float)i * yDiv);
    uint32_t _prStrSz = snprintf(NULL, 0, "%5.2f", _yLinSpV);
    char *_yLinSpS = (char *)calloc(_prStrSz + 1, sizeof(char));
    sprintf(_yLinSpS, "%5.2f", _yLinSpV);
    _h->DrawString(_h,
      self->plotArea_x - PLOT_VALUES_OFFSET,
      self->plotArea_y + i * (self->plotArea_h / PLOT_DIVISIONS_Y),
      _yLinSpS, 0, FONT_5X7, ALIGNMENT_RIGHT);
    free(_yLinSpS);
  }

  RGB565ImageDrawHandler_Delete(_h);
}

void _Plot_Single_2D(struct RGB565Plotter *self, double *x, double *fx, uint32_t n, uint16_t color) {
  // Draw Framing
  __Plot_2D_Create_Axes(self);

  // Get (x, y) area
  double xMin = _d_yminmax(x, n, false), xMax = _d_yminmax(x, n, true),
         yMin = _d_yminmax(fx, n, false), yMax = _d_yminmax(fx, n, true);

  // Values on Axes
  __Plot_2D_Create_Scale(self, xMin, xMax, yMin, yMax);

  // Scale (x, y) to G(x, y)
  float     _scPAW = self->plotArea_w - 3,
            _scPAH = self->plotArea_h - 2,
            _linAX = _scPAW / (xMax - xMin),
            _linBX = (xMin * _scPAW) / (xMin - xMax),
            _linAY = _scPAH / (yMax - yMin),
            _linBY = (yMin * _scPAH) / (yMin - yMax);

  // Do the actual plot (Linear Interpolation)
  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(self->img);
  for(uint32_t i = 0; i < n - 1; ++i) {
    // Plot area in x dir: self->plotArea_w, y dir: self->plotArea_h
    _h->DrawLine(_h,
      self->plotArea_x + 1 + _linAX * x[i] + _linBX,
      self->plotArea_y - 1 + (self->plotArea_h - (_linAY * fx[i] + _linBY)),
      self->plotArea_x + 1 + _linAX * x[i + 1] + _linBX,
      self->plotArea_y - 1 + (self->plotArea_h - (_linAY * fx[i + 1] + _linBY)),
    color);
  }

  RGB565ImageDrawHandler_Delete(_h);
}

void _Plot_Multiple_2D(struct RGB565Plotter *self, double *x, double **fx, uint32_t n, uint32_t nfx, uint16_t *colors) {
  // Draw Framing
  __Plot_2D_Create_Axes(self);

  double *_yMinArr = (double *)malloc(nfx * sizeof(double));
  double *_yMaxArr = (double *)malloc(nfx * sizeof(double));
  for(uint32_t i = 0; i < nfx; ++i) {
    _yMinArr[i] = _d_yminmax(fx[i], n, false);
    _yMaxArr[i] = _d_yminmax(fx[i], n, true);
  }

  // Get (x, y) area
  double xMin = _d_yminmax(x, n, false), xMax = _d_yminmax(x, n, true),
         yMin = _d_yminmax(_yMinArr, nfx, false), yMax = _d_yminmax(_yMaxArr, nfx, true); // Now take max of both arrays

  free(_yMinArr);
  free(_yMaxArr);

  // Values on Axes
  __Plot_2D_Create_Scale(self, xMin, xMax, yMin, yMax);

  // Scale (x, y) to G(x, y)
  float     _scPAW = self->plotArea_w - 3,
            _scPAH = self->plotArea_h - 2,
            _linAX = _scPAW / (xMax - xMin),
            _linBX = (xMin * _scPAW) / (xMin - xMax),
            _linAY = _scPAH / (yMax - yMin),
            _linBY = (yMin * _scPAH) / (yMin - yMax);

  // Do the actual plot (Linear Interpolation)
  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(self->img);
  for(uint32_t i = 0; i < n - 1; ++i) {
    for(uint32_t j = 0; j < nfx; ++j) {
      // Plot area in x dir: self->plotArea_w, y dir: self->plotArea_h
      _h->DrawLine(_h,
        self->plotArea_x + 1 + _linAX * x[i] + _linBX,
        self->plotArea_y - 1 + (self->plotArea_h - (_linAY * fx[j][i] + _linBY)),
        self->plotArea_x + 1 + _linAX * x[i + 1] + _linBX,
        self->plotArea_y - 1 + (self->plotArea_h - (_linAY * fx[j][i + 1] + _linBY)),
      colors[j]);
    }
  }

  RGB565ImageDrawHandler_Delete(_h);
}

void _Add_Title(struct RGB565Plotter *self, char *label, font_t fontSize) {
  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(self->img);
  _h->DrawString(_h, self->working_x + (self->working_w >> 1), self->working_y + (PLOT_INNER_OFFSET >> 1), label, 0, fontSize, ALIGNMENT_CENTER);
  RGB565ImageDrawHandler_Delete(_h);
}

void _Add_Label_X(struct RGB565Plotter *self, char *label) {
  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(self->img);
  _h->DrawString(_h, self->working_x + (self->working_w >> 1) + (PLOT_INNER_OFFSET >> 1), self->working_y + (self->working_h - (PLOT_INNER_OFFSET >> 1)), label, 0, FONT_5X7, ALIGNMENT_CENTER);
  RGB565ImageDrawHandler_Delete(_h);
}

void _Add_Label_Y(struct RGB565Plotter *self, char *label) {
  const uint8_t _f5x7w = 6, _f5x7h = 10;

  // Create label and write text to it
  RGB565Image *_i = RGB565Image_InitWithColor(_f5x7w * strlen(label) + 3, _f5x7h, 0xFFFF);
  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(_i);
  _h->DrawString(_h, 2, 2, label, 0, FONT_5X7, ALIGNMENT_LEFT);
  RGB565ImageDrawHandler_Delete(_h);

  RGB565Processor *_pi = RGB565Processor_Init(_i);
  _pi->Rotate(_pi, 90);
  RGB565Processor_Delete(_pi);

  RGB565Processor *_p = RGB565Processor_Init(self->img);
  _p->Insert(_p, _i, self->working_x + (PLOT_INNER_OFFSET >> 1), self->working_y + (self->working_h >> 1) - (_i->height >> 1), OP_REPLACE);
  RGB565Processor_Delete(_p);

  RGB565Image_Delete(_i);
}

RGB565Plotter *RGB565Plotter_Init(RGB565Image *img, uint16_t x, uint16_t y, uint16_t width, uint16_t height) {
  RGB565Plotter *self = (RGB565Plotter *)calloc(1, sizeof(RGB565Plotter));
  self->img = img;
  self->working_x = x;
  self->working_y = y;
  self->working_w = width;
  self->working_h = height;

  self->innerOffset = PLOT_INNER_OFFSET;

  self->plotArea_x = x + 2 * self->innerOffset;
  self->plotArea_y = y + self->innerOffset;
  self->plotArea_w = width - 3 * self->innerOffset;
  self->plotArea_h = height - (self->innerOffset << 1);

  RGB565ImageDrawHandler *_h = RGB565ImageDrawHandler_Init(self->img);
  // Working Area
  _h->DrawRectangle(_h, self->working_x, self->working_y, self->working_w, self->working_h, 0, LINE_STYLE_DOT);

  // Coordinate System Area
  _h->DrawRectangle(_h, self->plotArea_x, self->plotArea_y, self->plotArea_w, self->plotArea_h, 0, LINE_STYLE_SOLID);
  RGB565ImageDrawHandler_Delete(_h);

  self->Add_Title = _Add_Title;
  self->Plot_Single_2D = _Plot_Single_2D;
  self->Plot_Multiple_2D = _Plot_Multiple_2D;
  self->Add_Label_X = _Add_Label_X;
  self->Add_Label_Y = _Add_Label_Y;

  return self;
}

void RGB565Plotter_Delete(RGB565Plotter *self) {
  free(self);
}
