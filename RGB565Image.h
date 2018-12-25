#ifdef __cplusplus
extern "C" {
#endif

#ifndef RGB565IMAGE_H
#define RGB565IMAGE_H

#include <stdint.h>

// RGB 565 Image Class
typedef struct RGB565Image {
  uint16_t width, height;
  uint16_t **bitmap;

  void            (*Render)(struct RGB565Image *self, const char *path);
} RGB565Image;

// Init black image of size (width x height)
RGB565Image *RGB565Image_Init(uint16_t width, uint16_t height);

// Init image of size (width x height) with background color 'color'
RGB565Image *RGB565Image_InitWithColor(uint16_t width, uint16_t height, uint16_t color);

// Init image by loading from a saved file
RGB565Image *RGB565Image_InitWithFile(const char *path);

// Destructor
void RGB565Image_Delete(RGB565Image *self);

#endif

#ifdef __cplusplus
}
#endif
