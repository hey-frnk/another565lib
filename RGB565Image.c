#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include "RGB565Image.h"
#include "ImageException.h"

#define     RGB565_HEADER_SIZE 70

void _RGB565Image_Render(struct RGB565Image *self, const char *path) {
  // Attribute access overhead reduction
  uint16_t  _width    = self->width,
            _height   = self->height;

  // 1 padding byte fill if width is odd
  uint8_t   _paddingBytes = (_width & 0x01) ? 1 : 0;
  // BMP data size: 2 Bytes * (Width + Padding) * Height
  uint32_t  _bmpDataSize = ((_width + _paddingBytes) * _height) << 1,
            _bmpFileSize = RGB565_HEADER_SIZE + _bmpDataSize;

  uint8_t _header[RGB565_HEADER_SIZE] = {
    0x42, 0x4D,                         // 0x00 BMP Magic number
    // 0x02 File Size (Bytes)
    (uint8_t)(_bmpFileSize & 0xFF), (uint8_t)((_bmpFileSize >> 8) & 0xFF), (uint8_t)((_bmpFileSize >> 16) & 0xFF), (uint8_t)((_bmpFileSize >> 24) & 0xFF),
    0x00, 0x00, 0x00, 0x00,             // 0x06 Reserved
    0x46, 0x00, 0x00, 0x00,             // 0x0A 70 Bytes of Header

    0x38, 0x00, 0x00, 0x00,             // 0x0E 56 Bytes of Info
    // 0x12 Width in pixels (little endian)
    (uint8_t)(_width & 0xFF), (uint8_t)((_width >> 8) & 0xFF), 0x00, 0x00,
    // 0x16 Height in pixels (little endian)
    (uint8_t)(_height & 0xFF), (uint8_t)((_height >> 8) & 0xFF), 0x00, 0x00,
    0x01, 0x00,                         // 0x1A 1 plane
    0x10, 0x00,                         // 0x1C 16 bit RGB 565 image
    0x03, 0x00, 0x00, 0x00,             // 0x1E BI_BITFIELDS for 16 bit
    // 0x22 Image Data Size (Bytes)
    (uint8_t)(_bmpDataSize & 0xFF), (uint8_t)((_bmpDataSize >> 8) & 0xFF), (uint8_t)((_bmpDataSize >> 16) & 0xFF), (uint8_t)((_bmpDataSize >> 24) & 0xFF),
    0x12, 0x0B, 0x00, 0x00,             // 0x26 xPelsPerMeter
    0x12, 0x0B, 0x00, 0x00,             // 0x2A yPelsPerMeter
    0x00, 0x00, 0x00, 0x00,             // 0x2E 0 biClrUsed (all bits)
    0x00, 0x00, 0x00, 0x00,             // 0x32 0 biClrImportant (all colors)

    0x00, 0xF8, 0x00, 0x00,             // Red Mask
    0xE0, 0x07, 0x00, 0x00,             // Green Mask
    0x1F, 0x00, 0x00, 0x00,             // Blue Mask
    0x00, 0x00, 0x00, 0x00              // biClrUsed == 0 -> no table
  };

  uint8_t *_outputImage = (uint8_t *)malloc(_bmpFileSize * sizeof(uint8_t));
  uint32_t i;

  // Input header
  for(i = 0; i < RGB565_HEADER_SIZE; ++i) _outputImage[i] = _header[i];

  // Input Bitmap Data
  for(uint16_t y = 0; y < _height; ++y) {
    for(uint16_t x = 0; x < _width; ++x) {
      _outputImage[i++] = (uint8_t)(self->bitmap[_height - y - 1][x] & 0xFF); // Bottom left corner first
      _outputImage[i++] = (uint8_t)((self->bitmap[_height - y - 1][x] >> 8) & 0xFF);
    }
    // Only works for padding == 0 or == 1
    if(_paddingBytes) {
      _outputImage[i++] = 0;
      _outputImage[i++] = 0;
    }
  }

  FILE *_imageWriter = fopen(path, "w");
  fwrite(_outputImage, 1, _bmpFileSize, _imageWriter);
  fclose(_imageWriter);

  free(_outputImage);
}

// Methods of Image
RGB565Image *RGB565Image_Init(uint16_t width, uint16_t height) {
  RGB565Image *self = (RGB565Image *)calloc(1, sizeof(RGB565Image));
  if(!self) ThrowImageException(RGB565_IMAGE_EXCEPTION_MEM_ERR);

  self->width = width;
  self->height = height;

  self->bitmap = (uint16_t **)calloc(height, sizeof(uint16_t *));
  if(!self->bitmap) ThrowImageException(RGB565_IMAGE_EXCEPTION_MEM_ERR);

  for(uint16_t i = 0; i < height; ++i) {
    self->bitmap[i] = (uint16_t *)calloc(width, sizeof(uint16_t));
    if(!self->bitmap[i]) ThrowImageException(RGB565_IMAGE_EXCEPTION_MEM_ERR);
  }

  self->Render = _RGB565Image_Render;

  return self;
}

RGB565Image *RGB565Image_InitWithColor(uint16_t width, uint16_t height, uint16_t color) {
  RGB565Image *self = RGB565Image_Init(width, height);
  for(uint16_t i = 0; i < height; ++i) memset(self->bitmap[i], color, width * sizeof(uint16_t));
  return self;
}

RGB565Image *RGB565Image_InitWithFile(const char *path) {
  RGB565Image *self = (RGB565Image *)calloc(1, sizeof(RGB565Image));
  if(!self) ThrowImageException(RGB565_IMAGE_EXCEPTION_MEM_ERR);

  // Load image size and open
  FILE *_imageReader;
  _imageReader = fopen(path, "r");
  fseek(_imageReader, 0, SEEK_END);
  uint32_t _bmpFileSize = ftell(_imageReader);
  fseek(_imageReader, 0, SEEK_SET);

  uint8_t *_inputImage = malloc(_bmpFileSize * sizeof(uint8_t)), *_inputImagePtr;
  if(!_inputImage) ThrowImageException(RGB565_IMAGE_EXCEPTION_MEM_ERR);

  _inputImagePtr = _inputImage;
  fread(_inputImage, 1, _bmpFileSize, _imageReader);

  // Format check
  if(_inputImagePtr[0] != 'B' && _inputImagePtr[1] != 'M') ThrowImageException(RGB565_IMAGE_EXCEPTION_FILE_TYPE_ERR);
  _inputImagePtr += 10; // Skip reserved, file size

  // Format check 2:
  bool _bBMPchecker = _inputImagePtr[0] == 70 && _inputImagePtr[4] == 56;
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 8;

  self->width = (uint16_t)_inputImagePtr[0] | (((uint16_t)_inputImagePtr[1]) << 8);
  _inputImagePtr += 4;
  self->height = (uint16_t)_inputImagePtr[0] | (((uint16_t)_inputImagePtr[1]) << 8);
  _inputImagePtr += 4;

  _bBMPchecker &= _inputImagePtr[0] == 1 && !_inputImagePtr[1];
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 2;

  _bBMPchecker &= _inputImagePtr[0] == 16 && !_inputImagePtr[1];
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 2;

  _bBMPchecker &= _inputImagePtr[0] == 3 && !_inputImagePtr[1] && !_inputImagePtr[2] && !_inputImagePtr[3];
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 8; // skip data size

  // X & Y Pixels per meter check
  _bBMPchecker &= _inputImagePtr[0] == 0x12 && _inputImagePtr[1] == 0x0B;
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 4;
  _bBMPchecker &= _inputImagePtr[0] == 0x12 && _inputImagePtr[1] == 0x0B;
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 12;

  // Bit mask check
  _bBMPchecker &= _inputImagePtr[1] == 0xF8;
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 4;
  _bBMPchecker &= _inputImagePtr[0] == 0xE0 && _inputImagePtr[1] == 0x07;
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 4;
  _bBMPchecker &= _inputImagePtr[0] == 0x1F && _inputImagePtr[4] == 0x00;
  if(!_bBMPchecker) ThrowImageException(RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR);
  _inputImagePtr += 8;

  // Allocate space for new image
  self->bitmap = (uint16_t **)calloc(self->height, sizeof(uint16_t *));
  if(!self->bitmap) ThrowImageException(RGB565_IMAGE_EXCEPTION_MEM_ERR);

  uint32_t i = 0;
  uint8_t _paddingBytes = self->width & 0x01;
  for(uint16_t y = 0; y < self->height; ++y) {
    uint16_t _altY = self->height - y - 1;
    self->bitmap[_altY] = (uint16_t *)calloc(self->width, sizeof(uint16_t));
    if(!self->bitmap[_altY]) ThrowImageException(RGB565_IMAGE_EXCEPTION_MEM_ERR);

    // Subsequently read in image
    for(uint16_t x = 0; x < self->width; ++x, i += 2)
      self->bitmap[_altY][x] = (uint16_t)_inputImagePtr[i] | (((uint16_t)_inputImagePtr[i + 1]) << 8);
    if(_paddingBytes) i += 2;
  }

  fclose(_imageReader);
  free(_inputImage);

  // Init methods, finish
  self->Render = _RGB565Image_Render;
  return self;
}

void RGB565Image_Delete(RGB565Image *self) {
  for(uint16_t i = 0; i < self->height; ++i) free(self->bitmap[i]);
  free(self->bitmap);
  free(self);
}
