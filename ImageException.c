#include <stdio.h>
#include <stdlib.h>
#include "ImageException.h"

void ThrowImageException(RGB565_IMAGE_EXCEPTION e) {
  switch(e) {
    case RGB565_IMAGE_EXCEPTION_OUT_OF_BOUNDS:
      printf("RGB 565 Image Exception: Pixel Out of Bounds!\n"); break;
    case RGB565_IMAGE_EXCEPTION_MEM_ERR:
      printf("RGB 565 Image Exception: Memory Allocation Error!\n"); break;
    case RGB565_IMAGE_EXCEPTION_FILE_TYPE_ERR:
      printf("RGB 565 Image Exception: Image Format Incorrect!\n"); break;
    case RGB565_IMAGE_EXCEPTION_BMP_TYPE_ERR:
      printf("RGB 565 Image Exception: Bitmap is not of type 16 bit RGB 565!\n"); break;
    default:
      printf("RGB 565 Image Exception: Undefined Error!\n"); break;
  }
  printf("Terminating Program with Error Code %d\n", (int)e);
  exit(e);
}
