#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Processor.h"
#include "ImageException.h"

static const uint16_t __cnstmask[3] = {0xF800, 0x07E0, 0x001F};
static const uint8_t  __cnstshft[3] = {11, 5, 0};
static const float    __cnstfact[3] = {8.225806451612f, 4.047619047619f, 8.225806451612f};

static inline uint16_t __bilinearInterpolation565(double dx, double dy, uint16_t topLeft, uint16_t topRight, uint16_t bottomLeft, uint16_t bottomRight) {
  float t[3], b[3];
  uint16_t bIp[3];
  for(uint8_t i = 0; i < 3; ++i) {
    // Interpolate in x direction
    t[i] = (1.0f - dx) * ((((topLeft) & __cnstmask[i]) >> __cnstshft[i]) * __cnstfact[i]) + dx * ((((topRight) & __cnstmask[i]) >> __cnstshft[i]) * __cnstfact[i]);
    b[i] = (1.0f - dx) * ((((bottomLeft) & __cnstmask[i]) >> __cnstshft[i]) * __cnstfact[i]) + dx * ((((bottomRight) & __cnstmask[i]) >> __cnstshft[i]) * __cnstfact[i]);
    bIp[i] = (int16_t)roundf((1.0f - dy) * t[i] + dy * b[i]); // In y direction
    if(bIp[i] < 0) bIp[i] = 0;
    else if(bIp[i] > 255) bIp[i] = 255; // Clip
  }
  return ((bIp[0] & 0b11111000) << 8) | ((bIp[1] & 0b11111100) << 3) | (bIp[2] >> 3);
}

void _RGB565Processor_Rotate(struct RGB565Processor *self, float deg) {
  float   _radV = ((float)deg * M_PI) / 180.0f,
          _sinPhi = sinf(_radV), _cosPhi = cosf(_radV);

  uint16_t _NewHeight   = (uint16_t)roundf(fabsf(self->img->width * _sinPhi) + fabsf(self->img->height * _cosPhi)),
           _NewWidth    = (uint16_t)roundf(fabsf(self->img->width * _cosPhi) + fabsf(self->img->height * _sinPhi));

  uint16_t _cX = (self->img->width >> 1), _cY = (self->img->height >> 1);

  uint16_t **_rTarget = (uint16_t **)calloc(_NewHeight, sizeof(uint16_t *));
  for(uint16_t i = 0; i < _NewHeight; ++i) {
    _rTarget[i] = (uint16_t *)calloc(_NewWidth, sizeof(uint16_t));
    memset(_rTarget[i], 0xFFFF, _NewWidth * sizeof(uint16_t));
  }

  for(uint16_t i = 0; i < _NewHeight; ++i) {
    for(uint16_t j = 0; j < _NewWidth; ++j) {
      int16_t   x = j - (_NewWidth >> 1),
                y = (_NewHeight >> 1) - i;
      float     fAbs = hypotf((double)x, (double)y),
                fArg = 0.0f;

      if (x == 0) {
        if (y == 0) {
          _rTarget[i][j] = self->img->bitmap[_cY][_cX];
          continue;
        }
        else if (y < 0)   fArg = 1.5f * M_PI;
        else              fArg = 0.5f * M_PI;
      }
      else fArg = atan2f((double)y, (double)x);

      // Rotate
      fArg -= _radV;

      float _fx = fAbs * cosf(fArg);
      float _fy = fAbs * sinf(fArg);
      // convert Cartesian to raster
      _fx = _fx + _cX;
      _fy = _cY - _fy;
      int16_t _iffx = (int16_t)floorf(_fx),  _iffy = (int16_t)floorf(_fy),
              _icfx = (int16_t)ceilf(_fx),   _icfy = (int16_t)ceilf(_fy);

      // check bounds
      if(_iffx < 0 || _icfx < 0 || _iffx >= self->img->width || _icfx >= self->img->width \
        || _iffy < 0 || _icfy < 0 || _iffy >= self->img->height || _icfy >= self->img->height)
          continue;

      float _dx = _fx - (float)_iffx,
            _dy = _fy - (float)_iffy;

      _rTarget[i][j] = __bilinearInterpolation565(_dx, _dy,
        self->img->bitmap[_iffy][_iffx], self->img->bitmap[_iffy][_icfx],
        self->img->bitmap[_icfy][_iffx], self->img->bitmap[_icfy][_icfx]
      );
    }
  }

  // As the old image is single reference, realloc!
  for(uint16_t i = 0; i < self->img->height; ++i) free(self->img->bitmap[i]);
  free(self->img->bitmap);

  self->img->bitmap = _rTarget;
  self->img->width = _NewWidth;
  self->img->height = _NewHeight;
}

void _RGB565Processor_Insert(struct RGB565Processor *self, RGB565Image *second, uint16_t x, uint16_t y, img_op_t op) {
  for(uint16_t i = 0; i < second->height; ++i) {
    for(uint16_t j = 0; j < second->width; ++j) {
      if(i + y >= self->img->height || j + x >= self->img->width) ThrowImageException(RGB565_IMAGE_EXCEPTION_OUT_OF_BOUNDS);

      if(op == OP_REPLACE) self->img->bitmap[i + y][j + x] = second->bitmap[i][j];
      else {
        float cvs[3] = {0};
        int16_t cvt[3] = {0};
        for(uint8_t k = 0; k < 3; ++k) {
          cvs[k] = (((self->img->bitmap[i + y][j + x]) & __cnstmask[k]) >> __cnstshft[k]) * __cnstfact[k];
          if(op == OP_ADD) cvt[k] = (int16_t)roundf(cvs[k] + ((((second->bitmap[i][j]) & __cnstmask[k]) >> __cnstshft[k]) * __cnstfact[k]));
          else if(op == OP_SUBTRACT) cvt[k] = (int16_t)roundf(cvs[k] - ((((second->bitmap[i][j]) & __cnstmask[k]) >> __cnstshft[k]) * __cnstfact[k]));
          else if(op == OP_MULTIPLY) cvt[k] = (int16_t)roundf(cvs[k] * ((((second->bitmap[i][j]) & __cnstmask[k]) >> __cnstshft[k]) * __cnstfact[k]));
          else if(op == OP_DIVIDE) cvt[k] = (int16_t)roundf(cvs[k] / ((((second->bitmap[i][j]) & __cnstmask[k]) >> __cnstshft[k]) * __cnstfact[k]));
          if(cvt[k] < 0) cvt[k] = 0;
          else if(cvt[k] > 255) cvt[k] = 255; // Clip
        }
        self->img->bitmap[i + y][j + x] = ((cvt[0] & 0b11111000) << 8) | ((cvt[1] & 0b11111100) << 3) | (cvt[2] >> 3);
      }
    }
  }
}

void _RGB565Processor_Invert(struct RGB565Processor *self) {
  for(uint16_t i = 0; i < self->img->height; ++i)
    for(uint16_t j = 0; j < self->img->width; ++j)
      self->img->bitmap[i][j] = ~self->img->bitmap[i][j]; // Bit flip (lol)
}

void _RGB565Processor_Point_Curve_Pow(struct RGB565Processor *self, float a, uint8_t colorComponent) {
  for(uint16_t i = 0; i < self->img->height; ++i) {
    for(uint16_t j = 0; j < self->img->width; ++j) {
      uint16_t _resClr[3];
      for(uint8_t k = 0; k < 3; ++k) {
        double _clrExtract = ((((self->img->bitmap[i][j]) & __cnstmask[k]) >> __cnstshft[k]) * __cnstfact[k]);
        if((colorComponent >> k) & 0x01) _resClr[k] = (uint16_t)round(255.0f * pow(_clrExtract / 255.0, (double)a));
      }
      self->img->bitmap[i][j] = ((_resClr[0] & 0b11111000) << 8) | ((_resClr[1] & 0b11111100) << 3) | (_resClr[2] >> 3);
    }
  }
}

// Given n points (x_i, f(x_i)), return n - 1 cubic spline parameters a...d
// Numerical Analysis 9th ed - Burden, Faires (Ch. 3 Natural Cubic Spline, Pg. 149)
// Source: https://gist.github.com/svdamani/1015c5c4b673c3297309
static double **__nCubicInterpolation(double *x, double *a, uint8_t n) {
  int32_t i, j;
  n--;
  double h[n], A[n], l[n + 1], u[n + 1], z[n + 1],
        *c = malloc((n + 1) * sizeof(double)), *b = malloc(n * sizeof(double)), *d = malloc(sizeof(double));

  for (i = 0; i <= n - 1; ++i) h[i] = x[i + 1] - x[i];
  for (i = 1; i <= n - 1; ++i)
      A[i] = 3.0f * (a[i + 1] - a[i]) / h[i] - 3.0f * (a[i] - a[i - 1]) / h[i - 1];
  l[0] = 1.0f;
  u[0] = 0.0f;
  z[0] = 0.0f;
  for (i = 1; i <= n - 1; ++i) {
      l[i] = 2.0f * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
      u[i] = h[i] / l[i];
      z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
  }
  l[n] = 1;
  z[n] = 0;
  c[n] = 0;
  for (j = n - 1; j >= 0; --j) {
      c[j] = z[j] - u[j] * c[j + 1];
      b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0f * c[j]) / 3.0f;
      d[j] = (c[j + 1] - c[j]) / (3.0f * h[j]);
  }

  double **res = calloc(4, sizeof(double *));
  res[0] = a; res[1] = b; res[2] = c; res[3] = d;
  return res;
}

void _RGB565Processor_Point_Curve_Points(struct RGB565Processor *self, double *x, double *fx, uint8_t n, uint8_t colorComponent) {
  double **s = __nCubicInterpolation(x, fx, n);

  for(uint16_t i = 0; i < self->img->height; ++i) {
    for(uint16_t j = 0; j < self->img->width; ++j) {
      uint16_t _resClr[3];

      for(uint8_t k = 0; k < 3; ++k) {
        // Get pixel value
        double _clrExtract = ((((self->img->bitmap[i][j]) & __cnstmask[k]) >> __cnstshft[k]) * __cnstfact[k]);
        // Get spline position
        register volatile uint8_t l;
        for(l = 0; (l < n - 1) && (_clrExtract >= x[l]); ++l);
        --l;

        if((colorComponent >> k) & 0x01) {
          // a_l + b_l * (x - x_i) + c_l * (x - x_i)^2 + d_l * (x - x_i)^3
          double _x_xi = _clrExtract - x[l],
                 _clrIp = round(s[0][l] + s[1][l] * (_x_xi) + s[2][l] * pow(_x_xi, 2) + s[3][l] * pow(_x_xi, 3));
          _resClr[k] = (_clrIp > 255.0f) ? 255 : ((_clrIp < 0.0f) ? 0 : (uint16_t)_clrIp);
        }
      }
      self->img->bitmap[i][j] = ((_resClr[0] & 0b11111000) << 8) | ((_resClr[1] & 0b11111100) << 3) | (_resClr[2] >> 3);
    }
  }

  free(s[1]); free(s[2]); free(s[3]);
  free(s);
}

RGB565Processor *RGB565Processor_Init(RGB565Image *img) {
  RGB565Processor *self = (RGB565Processor *)calloc(1, sizeof(RGB565Processor));
  self->img = img;
  self->Rotate = _RGB565Processor_Rotate;
  self->Insert = _RGB565Processor_Insert;
  self->Invert = _RGB565Processor_Invert;
  self->Point_Curve_Pow = _RGB565Processor_Point_Curve_Pow;
  self->Point_Curve_Points = _RGB565Processor_Point_Curve_Points;
  return self;
}

void RGB565Processor_Delete(RGB565Processor *self) {
  free(self);
}
