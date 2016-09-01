/*****************************************

 GL_Display.h

 General purpose frames per second counter for OpenGL/GLUT GNU/Linux
 programs. Displays "Frames per second: N" at an arbitrary position in
 the window. Saves and restores the app's modelview and projection matrices,
 colour, and lighting.

 Author: Toby Howard. toby@cs.man.ac.uk.
 Version 2.1, 1 June 2006

 ====================

 Usage: to add an on-screen frames per second counter to your program, save
 this file alongside your app as "frames.h", and add:

    #include "GL_Display.h"

 immediately after all your app's other #includes; then bracket all the
 code in your display() function, before swapping buffers, with

   frameStart();

 and

   frameEnd(void *font, GLclampf r, GLclampf g, GLclampf b,
            float x, float y);

     font:    font to use, e.g., GLUT_BITMAP_HELVETICA_10
     r, g, b: RGB text colour
     x, y:    text position in window: range [0,0] (bottom left of window)
              to [1,1] (top right of window).

 ====================

 Example:

    void display(void) {
      glClear(GL_COLOR_BUFFER_BIT);

      frameStart();

      // all the graphics code

      frameEnd(GLUT_BITMAP_HELVETICA_10, 1.0, 1.0, 1.0, 0.05, 0.95);

      glutSwapBuffers();
    }
*****************************************/
#pragma once

#include <ctime>
#include <cmath>
#include <cstdio>

namespace SF {

  static const unsigned int DISPLAY_COUNT = 199;

  static clock_t _frameStart;
  static float _deltaTime = 0;
  static float _fps;
  static unsigned int _counter = 0;

  static char _fpsStr[64];
  static const char *_ch;
  static GLint _matrixMode;

  inline void displayOnScreen (void *font, GLclampf r, GLclampf g, GLclampf b, GLfloat x, GLfloat y, const char *str)
  {
    glDisable (GL_LIGHTING);

    glGetIntegerv (GL_MATRIX_MODE, &_matrixMode);

    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glLoadIdentity ();
    gluOrtho2D (0., 1., 0., 1.);

    glMatrixMode (GL_MODELVIEW);
    glPushMatrix ();
    glLoadIdentity ();
    glPushAttrib (GL_COLOR_BUFFER_BIT);
    glColor3f (r, g, b);
    glRasterPos3f (x, y, 0.);

    for (_ch = str; *_ch; ++_ch){
      glutBitmapCharacter (font, static_cast <int> (*_ch));
    }

    glPopAttrib ();
    glPopMatrix ();
    glMatrixMode (GL_PROJECTION);
    glPopMatrix ();
    glMatrixMode (_matrixMode);

    glEnable (GL_LIGHTING);
  }

  // function to start frame
  void
  frameStart ()
  {
    _frameStart = clock ();
  }

  // function to end frame
  void
  frameEnd (void *font, GLclampf r, GLclampf g, GLclampf b, GLfloat x, GLfloat y)
  {
    _deltaTime += static_cast <float> (clock () - _frameStart)/ CLOCKS_PER_SEC;
    ++_counter;

    if (_counter > DISPLAY_COUNT){
      if (_deltaTime < EPSILON){
        _fps = __builtin_inf();
      } else {
        _fps = 200./ _deltaTime;
        _counter = 0;
        _deltaTime = 0;
      }
    }

    sprintf (_fpsStr, "FPS %g", _fps);
    displayOnScreen (font, r, g, b, x, y, _fpsStr);
  }
}
