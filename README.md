# Another RGB 565 Image Library
This is yet another RGB 565 image library. It has a whole lot of fun features, let me explain.

 **RGB565Image.h/c:**  
- This is the basic image 'class'. It can read and write 16 bit RGB565 BMP files or create new black or custom initialized images
- Depends on: `stdint.h`, `stdlib.h`, `string.h`, `stdbool.h`, `ImageException.h`

**Processor.h/c:**
 - The processor can do basic image manipulations
 - It depends on `RGB565Image.h`, `stdio.h`, `stdlib.h`, `string.h`, `math.h`, `ImageException.h`
