#include <stdlib.h>
#include <string.h>
#include "DrawHandler.h"
#include "ImageException.h"

// Methods of Image Draw Handler
void _DrawPixel(struct RGB565ImageDrawHandler *self, uint16_t x, uint16_t y, uint16_t color) {
  if(x < 0 || y < 0 || x >= self->img->width || y >= self->img->height) ThrowImageException(RGB565_IMAGE_EXCEPTION_OUT_OF_BOUNDS);
  self->img->bitmap[y][x] = color;
}

void _DrawHLine(struct RGB565ImageDrawHandler *self, uint16_t x, uint16_t y, uint16_t l, uint16_t color, line_style_t s) {
  switch(s) {
    case LINE_STYLE_SOLID: {
      for(uint16_t i = 0; i < l; ++i) self->DrawPixel(self, x + i, y, color);
      break;
    }
    case LINE_STYLE_DOT: {
      for(uint16_t i = 0; i < l; ++i) if(i % 2) self->DrawPixel(self, x + i, y, color);
      break;
    }
    case LINE_STYLE_DASH: {
      for(uint16_t i = 0; i < l; ++i)
        if(i % (LINE_STYLE_LENGTH << 1) < LINE_STYLE_LENGTH)
          self->DrawPixel(self, x + i, y, color);
      break;
    }
    case LINE_STYLE_DASHDOT: {
      int16_t plottedPx = 0, i = 0;
      while(i < l) {
        if(i % (LINE_STYLE_LENGTH << 1) < LINE_STYLE_LENGTH && plottedPx++ < LINE_STYLE_LENGTH) {
          self->DrawPixel(self, x + i++, y, color);
        } else {
          ++i;
          self->DrawPixel(self, x + i++, y, color);
          plottedPx = 0;
        }
      }
      break;
    }
  }
}

void _DrawVLine(struct RGB565ImageDrawHandler *self, uint16_t x, uint16_t y, uint16_t l, uint16_t color, line_style_t s) {
  switch(s) {
    case LINE_STYLE_SOLID: {
      for(uint16_t i = 0; i < l; ++i) self->DrawPixel(self, x, y + i, color);
      break;
    }
    case LINE_STYLE_DOT: {
      for(uint16_t i = 0; i < l; ++i) if(i % 2) self->DrawPixel(self, x, y + i, color);
      break;
    }
    case LINE_STYLE_DASH: {
      for(uint16_t i = 0; i < l; ++i)
        if(i % (LINE_STYLE_LENGTH << 1) < LINE_STYLE_LENGTH)
          self->DrawPixel(self, x, y + i, color);
      break;
    }
    case LINE_STYLE_DASHDOT: {
      int16_t plottedPx = 0, i = 0;
      while(i < l) {
        if(i % (LINE_STYLE_LENGTH << 1) < LINE_STYLE_LENGTH && plottedPx++ < LINE_STYLE_LENGTH) {
          self->DrawPixel(self, x, y + i++, color);
        } else {
          ++i;
          self->DrawPixel(self, x, y + i++, color);
          plottedPx = 0;
        }
      }
      break;
    }
  }
}

void _DrawLine(struct RGB565ImageDrawHandler *self, uint16_t x0, uint16_t y0, uint16_t x1, uint16_t y1, uint16_t color) {
  // According to Bresenham algorithm
  int16_t dx =  abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int16_t dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int16_t err = dx + dy, e2;

  while(1) {
    self->DrawPixel(self, x0, y0, color);
    if(x0 == x1 && y0 == y1) break;
    e2 = err << 1;
    if(e2 > dy) { err += dy; x0 += sx; } /* e_xy + e_x > 0 */
    if(e2 < dx) { err += dx; y0 += sy; } /* e_xy + e_y < 0 */
  }
}

void _DrawRectangle(struct RGB565ImageDrawHandler *self, uint16_t x, uint16_t y, uint16_t width, uint16_t height, uint16_t color, line_style_t s) {
  self->DrawHLine(self, x, y, width, color, s);
  self->DrawHLine(self, x, y + height - 1, width, color, s);
	self->DrawVLine(self, x, y, height, color, s);
	self->DrawVLine(self, x + width - 1, y, height, color, s);
}

void _DrawFilledRectangle(struct RGB565ImageDrawHandler *self, uint16_t x, uint16_t y, uint16_t width, uint16_t height, uint16_t color) {
  if(width > height) for(uint8_t i = 0; i < height; ++i) self->DrawHLine(self, x, y + i, width, color, LINE_STYLE_SOLID);
	else for(uint8_t i = 0; i < width; ++i) self->DrawVLine(self, x + i, y, height, color, LINE_STYLE_SOLID);
}

static inline void _DrawChar_retrieveWidthHeight(font_t fontSize, uint8_t *fWidth, uint8_t *fHeight, uint8_t *xspPd, unsigned char **_cptr){
	switch(fontSize){
		case FONT_5X7: 		*fWidth = 5; 		*fHeight = 7; 	*xspPd = 1;		if(_cptr) *_cptr = (unsigned char *)Font5x7;		break;
		case FONT_8X12: 	*fWidth = 8; 		*fHeight = 12; 	*xspPd = 0;		if(_cptr) *_cptr = (unsigned char *)Font8x12; 	break;
		case FONT_8X14: 	*fWidth = 8; 		*fHeight = 14; 	*xspPd = 1;		if(_cptr) *_cptr = (unsigned char *)Font8x14; 	break;
		case FONT_12X16: 	*fWidth = 12;		*fHeight = 16; 	*xspPd = 1;		if(_cptr) *_cptr = (unsigned char *)Font12x16; 	break;
		// case 3: *fWidth = 16; *fHeight = 26; if(_cptr) *_cptr = (unsigned char *)Font16x26; break;
	}
}

void _DrawChar(struct RGB565ImageDrawHandler *self, uint16_t x, uint16_t y, uint8_t chr, uint16_t color, font_t fontSize) {
  uint8_t fWidth = 0, fHeight = 0, xspPd;
	unsigned char *_cptr;
	_DrawChar_retrieveWidthHeight(fontSize, &fWidth, &fHeight, &xspPd, &_cptr);

	// Aaaah, the perks of copying from different libraries!
	if(fontSize == 0){
		uint16_t _ypj;
		uint8_t *buffer = &_cptr[(chr - 0x20) * fWidth];
		for(uint8_t j = 0; j < fHeight; ++j){
			_ypj = y + j;
			for(uint8_t i = 0; i < fWidth; ++i) if((buffer[i] >> j) & 0x01) self->DrawPixel(self, x + i, _ypj, color);
		}
	} else {
		uint8_t _ishn;
    uint16_t _ypi, _xpfw = x + fWidth, _fwmj;
    uint8_t _fwecp = fWidth > 8,
            // Extra sort of alignment?
            n           = _fwecp ? 1 : 0,
            _rshAdj     = _fwecp ? 8 : fWidth,
            _fHeight    = fHeight << n,
            *buffer     = &_cptr[(chr - 0x20) * _fHeight];

		for(uint8_t i = 0; i < fHeight; ++i){
			_ishn = i << n;
			_ypi = i + y;
			for(uint8_t j = 0; j < _rshAdj; ++j){
				_fwmj = _xpfw - j;
				for(uint8_t k = 0; k <= n; ++k) if((buffer[_ishn + k] >> j) & 0x01) self->DrawPixel(self, _fwmj - (k << 3), _ypi, color);
			}
		}
	}
}

void _DrawString(struct RGB565ImageDrawHandler *self, uint16_t x, uint16_t y, const char *str, uint16_t color, font_t fontSize, alignment_t alignment) {
  // cL tells the length of str
  uint8_t fWidth = 0, fHeight = 0, xspPd = 0;
  _DrawChar_retrieveWidthHeight(fontSize, &fWidth, &fHeight, &xspPd, NULL);

  // Length of string
  uint16_t cL = strlen(str);

  switch(alignment){
    case ALIGNMENT_CENTER: 	x = x - (((fWidth + xspPd) * cL) >> 1); break;
    case ALIGNMENT_RIGHT:		x = x - ((fWidth + xspPd) * cL); 				break;
    default: break;
  }

  // _x, _y are local backups of X and Y
  uint16_t _x = x; // _y = y;

  while (*str) {
    self->DrawChar(self, x, y, *str++, color, fontSize);
    // Drawing inside
    if(x < self->img->width - ((fWidth + xspPd) << 1)) x += (fWidth + xspPd);
    // Word wrap
    else if (y < self->img->height - ((fHeight + 1) << 1)) {x =_x; y += (fHeight + 1);}
    // Reset otherwise
    else {x =_x;} //y =_y;}
  }
}


RGB565ImageDrawHandler *RGB565ImageDrawHandler_Init(RGB565Image *img) {
  RGB565ImageDrawHandler *self = (RGB565ImageDrawHandler *)calloc(1, sizeof(RGB565ImageDrawHandler));
  self->img = img;
  self->DrawPixel = _DrawPixel;
  self->DrawHLine = _DrawHLine;
  self->DrawVLine = _DrawVLine;
  self->DrawLine = _DrawLine;
  self->DrawRectangle = _DrawRectangle;
  self->DrawFilledRectangle = _DrawFilledRectangle;
  self->DrawChar = _DrawChar;
  self->DrawString = _DrawString;
  return self;
}

void RGB565ImageDrawHandler_Delete(RGB565ImageDrawHandler *self) {
  free(self);
}
