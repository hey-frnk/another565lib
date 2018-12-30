#ifdef __cplusplus
extern "C" {
#endif

#ifndef _CONFIG_H
#define _CONFIG_H

// Basic Colors
#define   COLOR_WHITE        0xFFFF
#define   COLOR_BLACK        0x0000
#define   COLOR_BLUE         0x001F
#define   COLOR_BRED         0xF81F
#define   COLOR_GRED         0xFFE0
#define   COLOR_GBLUE        0x07FF
#define   COLOR_RED          0xF800
#define   COLOR_MAGENTA      0xF81F
#define   COLOR_GREEN        0x07E0
#define   COLOR_CYAN         0x7FFF
#define   COLOR_YELLOW       0xFFE0
#define   COLOR_BROWN        0xBC40
#define   COLOR_BRRED        0xFC07
#define   COLOR_GRAY         0x8430
#define   COLOR_DARKBLUE     0x01CF
#define   COLOR_LIGHTBLUE    0x7D7C
#define   COLOR_GRAYBLUE     0x5458
#define   COLOR_LIGHTGREEN   0x841F
#define   COLOR_LGRAY        0xC618
#define   COLOR_LGRAYBLUE    0xA651
#define   COLOR_LBBLUE       0x2B12

// ###################################### Processor
// Operations, color selection
typedef enum { OP_REPLACE, OP_ADD, OP_SUBTRACT, OP_MULTIPLY, OP_DIVIDE, OP_ABSDIFF } img_op_t;
typedef enum { SCALE_PLAIN, SCALE_BILINEAR, SCALE_BICUBIC } scale_t;

#define       OP_RED      0x01    // 0001
#define       OP_GREEN    0x02    // 0010
#define       OP_BLUE     0x04    // 0100
#define       OP_RGB      0x07    // 0111

// ###################################### Draw Handler
// Font, font alignment, line style
typedef enum { FONT_5X7, FONT_8X12, FONT_8X14, FONT_12X16 } font_t;
typedef enum { ALIGNMENT_LEFT, ALIGNMENT_CENTER, ALIGNMENT_RIGHT } alignment_t;
typedef enum { LINE_STYLE_SOLID, LINE_STYLE_DOT, LINE_STYLE_DASH, LINE_STYLE_DASHDOT } line_style_t;

#define LINE_STYLE_LENGTH 4

// ###################################### Plotter
// Aesthetic parameters
#define PLOT_INNER_OFFSET   45
#define PLOT_DIVISIONS_X    10
#define PLOT_DIVISIONS_Y    5
#define PLOT_MARK_LENGTH    10
#define PLOT_VALUES_OFFSET  8



#endif

#ifdef __cplusplus
}
#endif
