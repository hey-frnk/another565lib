# Another RGB 565 Image Library
This is yet another RGB 565 image library in plain C (plain C OOP) using C standard libraries. It has a whole lot of fun features, more or less efficient, let me explain.

 **RGB565Image.h/c:**  
- This is the basic image 'class'. It can read and write 16 bit RGB565 BMP files or create new black or custom initialized images
- Depends on: `stdint.h`, `stdlib.h`, `string.h`, `stdbool.h`, `ImageException.h`

**Processor.h/c:**
 - The processor can do basic image manipulations
 - Functions `Rotate` and `Insert` are self explaining
 - Use `Point_Curve_Points` to make a picture instagrammable. Define a custom point curve 
 - It depends on `RGB565Image.h`, `stdio.h`, `stdlib.h`, `string.h`, `math.h`, `ImageException.h`

**DrawHandler.h/c:**
 - Use the draw handler to create simple objects
 - Basic geometric functions: Draw a single pixel, a line (horizontal, vertical or point to point) and rectangles
 - Text functions: Print single characters or custom aligned texts using built in fonts
 - Depends on: `RGB565Image.h`, `stdlib.h`, `string.h`, `ImageException.h`

**Plotter.h/c:**

 - Plotter is a very specific tool to create 2D graphs demonstrating the ability of  `DrawHandler` and `Processor`
 - Use Single `(x, f(x))` or Multiple `(x, f1(x) ... fn(x))` to graph and the title and label functions to add caption
 - Depends on: `RGB565Image.h`, `stdio.h`, `stdlib.h`, `string.h`, `math.h`, `stdbool.h`, `DrawHandler.h`, `Processor.h`, `ImageException.h`

More on the way. `Main.c` is a small preview of what's possible (or not).
