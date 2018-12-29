#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "RGB565Image.h"
#include "Processor.h"
#include "DrawHandler.h"
#include "Plotter.h"
#include "DrawHandler.h"
#include "ImageException.h"

/*  cd documents/github/another565lib
    gcc RGB565Image.c Processor.c DrawHandler.c Plotter.c ImageException.c Main.c -O3 -std=c17 -Wextra -o Main
    gcc RGB565Image.c Processor.c DrawHandler.c Plotter.c ImageException.c Main.c -O0 -g3 -std=c17 -Wextra -o Main
    valgrind ./Main instagrammable2.jpg
*/

// Random ability demo

int main(int argc, char **argv) {
  // New Image
  int w = 833, h = 639;
  RGB565Image *img = RGB565Image_InitWithColor(w, h, 0xFFFF);

  // Generate Data
  const int samples = 512;
  double timeSample = 1.0f;
  double f = 1.4;
  double signalArray[samples];
  double outputArray[samples];
  double output2[samples];
  double *outputArrayArray[2] = {outputArray, output2};
  uint16_t colors[2] = {0xF800, 0x001F};
  for(int i = 0; i < samples; i++) {
    signalArray[i] = (i / (double)samples) * timeSample;
    outputArray[i] = sin(2 * M_PI * f * signalArray[i]) + 0.5 * sin(2 * M_PI * 4 * f * signalArray[i]) + 3 * sin(2 * M_PI * 16 * f * signalArray[i]) + 0.5 * sin(2 * M_PI * 47 * f * signalArray[i] + 0.84) + 2.23 * sin(2 * M_PI * 78 * f * signalArray[i] - 0.30) + 2.08 * cos(2 * M_PI * 101 * f * signalArray[i] + 0.9) + 1.87;
    output2[i] = sin(6.44 * M_PI * f * signalArray[i]);
  }

  // Plot data to image
  RGB565Plotter *p = RGB565Plotter_Init(img, 2, 2, w - 4, 200);
  RGB565Plotter *p2 = RGB565Plotter_Init(img, 2, 2 + 210, w - 4, 200);
  p->Plot_Single_2D(p, signalArray, outputArray, samples, 0x5458);
  p->Add_Title(p, (char *)"This is an amazing function", FONT_8X12);
  p->Add_Label_X(p, (char *)"t [s]");
  p->Add_Label_Y(p, (char *)"Output Voltage [V]");
  p2->Plot_Multiple_2D(p2, signalArray, outputArrayArray, samples, 2, colors);
  p2->Add_Title(p2, (char *)"This is an even more amazing function", FONT_8X12);
  p2->Add_Label_X(p2, (char *)"f [Hz]");
  p2->Add_Label_Y(p2, (char *)"Output amazingness [kV]");
  RGB565Plotter_Delete(p);
  RGB565Plotter_Delete(p2);

  RGB565Processor *processor = RGB565Processor_Init(img);
  processor->Rotate(processor, 10);
  // processor->Rotate(processor, -10);
  RGB565Processor_Delete(processor);

  img->Render(img, (char *)"Eyeyetwo.bmp");
  RGB565Image_Delete(img);


  if(argc <= 1) {
    printf("Oops, use the second parameter for the image!\n");
    ThrowImageException(RGB565_IMAGE_NO_ERROR);
  }
  RGB565Image *img2 = RGB565Image_InitWithFile(argv[1]);

  RGB565Plotter *k = RGB565Plotter_Init(img2, 2, 2, img2->width - 4, 200);
  k->Plot_Multiple_2D(k, signalArray, outputArrayArray, samples, 2, colors);
  k->Add_Title(k, (char *)"How awesome is this. You can overlay a graph, wohooooo!", FONT_8X12);
  k->Add_Label_X(k, (char *)"f [MHz] #reasontoroam #TheVFDCollective");
  k->Add_Label_Y(k, (char *)"Output amazingness [GB]");
  RGB565Plotter_Delete(k);

  RGB565Processor *pro2 = RGB565Processor_Init(img2);
  pro2->Point_Curve_Pow(pro2, 0.9f, OP_RGB);
  double px[5] = {0, 100, 180, 220, 255};
  double py[5] = {55, 100, 180, 210, 220};

  pro2->Point_Curve_Points(pro2, px, py, 5, OP_RGB);

  pro2->Rotate(pro2, sqrt(3) + log(2.81)); // lol
  pro2->Rotate(pro2, -(sqrt(3) + log(2.81))); // lolwut
  pro2->Dreamify(pro2, 10);

  RGB565Processor_Delete(pro2);

  img2->Render(img2, "output.bmp");
  RGB565Image_Delete(img2);
}
