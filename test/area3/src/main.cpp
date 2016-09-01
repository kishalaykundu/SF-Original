/**
 * @file main.cpp
 * @author Kishalay Kundu <kishalay.kundu@gmail.com>
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 */

#include <cassert>
#include <cstdio>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>

#include "Preprocess.h"

extern "C" {
#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glxext.h>
#include <GL/glut.h>

#include <cuda.h>
#include <cuda_gl_interop.h>
#include <cudaGL.h>
}

#include <boost/thread/thread.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

#include "vec4.h"
#include "/home/kish1/Projects/SF/common/GL/common.h"
#include "/home/kish1/Projects/SF/common/CUDA/common.h"

static int GLX_ATTRIBUTE_LIST [] = {GLX_RGBA, None};

using namespace std;
using namespace boost;
using namespace SF;

extern "C" {
#include <cuda_runtime_api.h>
}

// Data-related parameters
typedef struct mesh_t{
  vector <vec4> verts0;
  vector <vec4> verts1;
  vector <unsigned int> inds;
} Mesh;

static Mesh mesh;
boost::thread myThread;

boost::interprocess::interprocess_semaphore *mutexes [2];

// OpenGL-related parameters
bool bFlag = false;
GLuint windowWidth = 0, windowHeight = 0;
GLuint vId1 = 0, vId2 = 0, vaId1 = 0, vaId2 = 0, iId = 0;
GLuint pId = 0;

string vs_prog = string
("#version 410\nin vec4 pos;\nvoid main ()\n{\nvec4 position = vec4 ((pos.x + 150.)/100., (pos.y - 110.)/50., (pos.z + 50.)/150., 1.);\ngl_Position = position;\n}\n");
string fs_prog = string ("#version 410\nout vec4 fragColor;\nvoid main ()\n{\nfragColor = vec4 (.5, .5, .5, 1.);\n}\n");

GLXContext gl_context, cugl_context;
Display *gl_display;
GLXDrawable gl_drawable;

// CUDA-related parameters
CUcontext cu_context;
CUfunction cu_func;

CUdeviceptr devBufferPtr [2];
CUgraphicsResource cuBufferId [2];

// CUDA-thread related functions
static void run ()
{
  size_t numBytes = 0;

  XVisualInfo *visualInfo = glXChooseVisual (gl_display, 0, GLX_ATTRIBUTE_LIST);
  GLenum error;
  checkGLError (error);

  cugl_context = glXCreateContext (gl_display, visualInfo, gl_context, GL_TRUE);
  checkGLError (error);
  if (!glXMakeCurrent (gl_display, gl_drawable, cugl_context)){
    checkGLError (error);
    PRINT ("glXMakeCurrent error\n");
  }

  CUresult status = cuInit (0);
  printCUResult (status);

  // check if any device is CUDA capable
  int numCudaDevices = 0;
  status = cuDeviceGetCount (&numCudaDevices);
  printCUResult (status);
  if (!numCudaDevices){
    PRINT ("Error: Could not find CUDA devices\n");
  }

  // get handle for device 0
  CUdevice cudaDevice;
  status = cuDeviceGet (&cudaDevice, 0);
  printCUResult (status);

  // create context
  status = cuGLCtxCreate (&cu_context, CU_CTX_SCHED_BLOCKING_SYNC, cudaDevice);
  printCUResult (status);

  status = cuCtxPushCurrent (cu_context);
  printCUResult (status);

  status = cuCtxSynchronize ();
  printCUResult (status);

  CUmodule cu_mod;
  status = cuModuleLoad (&cu_mod, "/home/kish1/Projects/bin/CudaXFEM_xfem.cu.ptx");
  printCUResult (status);
  cuModuleGetFunction (&cu_func, cu_mod, "conjugate_gradient");
  cuGraphicsGLRegisterBuffer (&cuBufferId [0], vId1, CU_GRAPHICS_MAP_RESOURCE_FLAGS_NONE);
  cuGraphicsGLRegisterBuffer (&cuBufferId [1], vId2, CU_GRAPHICS_MAP_RESOURCE_FLAGS_NONE);

  cuCtxPopCurrent (&cu_context);

  int threadsPerBlock = 32;
  unsigned int problemsize = mesh.verts0.size ();
  void *kernelArgs1 [] = { &(devBufferPtr [0]), &(devBufferPtr [1]), &problemsize};
  void *kernelArgs2 [] = { &(devBufferPtr [1]), &(devBufferPtr [0]), &problemsize};

  unsigned int nPrints = 0;

  while (true){
    mutexes [1]->wait ();
    cuCtxPushCurrent (cu_context);

    cuGraphicsMapResources (2, cuBufferId, 0);
    cuGraphicsResourceGetMappedPointer (&(devBufferPtr [0]), &numBytes, cuBufferId [0]);
    if (nPrints < 3){
      fprintf (stdout, "%lu ", numBytes); fflush (stdout);
    }
    cuGraphicsResourceGetMappedPointer (&(devBufferPtr [1]), &numBytes, cuBufferId [1]);
    if (nPrints < 3){
      fprintf (stdout, "%lu\n", numBytes); fflush (stdout);
      ++nPrints;
    }

    // launch CUDA kernel
    if (bFlag){
      cuLaunchKernel (cu_func, problemsize/ threadsPerBlock, 1, 1, threadsPerBlock, 1, 1, 0, 0, kernelArgs1, 0);
    } else {
      cuLaunchKernel (cu_func, problemsize/ threadsPerBlock, 1, 1, threadsPerBlock, 1, 1, 0, 0, kernelArgs2, 0);
    }

    // unmap GL buffers so that OpenGL can use them
    cuGraphicsUnmapResources (2, cuBufferId, 0);

    cuCtxPopCurrent (&cu_context);
    mutexes [0]->post ();
  }
}

// GLUT-related functions
static void display ()
{
  mutexes [0]->wait ();

  glUseProgram (pId);

  glClearColor (0., 0., 0., 0.);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (bFlag){
    glBindVertexArray (vaId1);
  } else {
    glBindVertexArray (vaId2);
  }
  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, iId);
  glDrawElements (GL_TRIANGLES, mesh.inds.size (), GL_UNSIGNED_INT, 0);
  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
  glBindVertexArray (0);

  glutSwapBuffers ();
  glUseProgram (0);

  bFlag = !bFlag;

  mutexes [1]->post ();
}

static void idle ()
{
  display ();
}

static void keys (unsigned char k, int a, int b)
{
  switch (k) {
    case 'q': case 'Q':
      exit (EXIT_SUCCESS);
  }
}

static void resize (int w, int h)
{
  windowWidth = w;
  windowHeight = h;

  char wStr [16];
  sprintf (wStr, "%d x %d", w, h);

  string title ("Test Area3 System - ");
  title.append (wStr);

  glutSetWindowTitle (title.c_str ());

  display ();
}

static void initGL ()
{
  glGenBuffers (1, &vId1);
  glBindBuffer (GL_ARRAY_BUFFER, vId1);
  glBufferData (GL_ARRAY_BUFFER, 4*sizeof (real)*mesh.verts0.size (), &(mesh.verts0 [0]._v[0]), GL_DYNAMIC_DRAW);
  glBindBuffer (GL_ARRAY_BUFFER, 0);

  glGenBuffers (1, &vId2);
  glBindBuffer (GL_ARRAY_BUFFER, vId2);
  glBufferData (GL_ARRAY_BUFFER, 4*sizeof (real)*mesh.verts1.size (), &(mesh.verts1 [0]._v[0]), GL_DYNAMIC_DRAW);
  glBindBuffer (GL_ARRAY_BUFFER, 0);

  glGenBuffers (1, &iId);
  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, iId);
  glBufferData (GL_ELEMENT_ARRAY_BUFFER, 3*sizeof (unsigned int)*mesh.inds.size (), &(mesh.inds [0]), GL_DYNAMIC_DRAW);
  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, 0);
}

static void initGPUProgs ()
{
  pId = glCreateProgram ();

  GLuint shader = 0;
  shader = glCreateShader (GL_VERTEX_SHADER);

  const char *vSource = vs_prog.c_str ();
  glShaderSource (shader, 1, &vSource, NULL);
  glCompileShader (shader);
  glAttachShader (pId, shader);

  GLint compileStatus;
  glGetShaderiv (shader, GL_COMPILE_STATUS, &compileStatus);
  if (compileStatus == GL_FALSE){

    char log [1024];
    GLsizei loglength;
    glGetShaderInfoLog (shader, 1024, &loglength, log);
    PRINT ("[1]Shader Info Log\n%s\n", log);
  }

  shader = glCreateShader (GL_FRAGMENT_SHADER);

  const char *fSource = fs_prog.c_str ();
  glShaderSource (shader, 1, &fSource, NULL);
  glCompileShader (shader);
  glAttachShader (pId, shader);

  glGetShaderiv (shader, GL_COMPILE_STATUS, &compileStatus);
  if (compileStatus == GL_FALSE){

    char log [1024];
    GLsizei loglength;
    glGetShaderInfoLog (shader, 1024, &loglength, log);
    PRINT ("[2]Shader Info Log\n%s\n", log);
  }

  glLinkProgram (pId);

  glUseProgram (pId);
  GLint posLoc = glGetAttribLocation (pId, "pos");
  glBindFragDataLocation (pId, 0, "fragColor");

  glGenVertexArrays (1, &vaId1);
  glBindVertexArray (vaId1);
  glBindBuffer (GL_ARRAY_BUFFER, vId1);
  glVertexAttribPointer (posLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray (posLoc);

  glGenVertexArrays (1, &vaId2);
  glBindVertexArray (vaId2);
  glBindBuffer (GL_ARRAY_BUFFER, vId2);
  glVertexAttribPointer (posLoc, 4, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray (posLoc);

  glBindBuffer (GL_ARRAY_BUFFER, 0);
  glBindVertexArray (0);
  glUseProgram (0);
}

// Data-related functions
static void readCube ()
{
  string prefix ("/home/kish1/Data/Cube/fem/0/cube.");

  string filename (prefix);
  filename.append ("node");

  FILE *fp = fopen (filename.c_str (), "r");

  unsigned int nelems = 0;
  int status = fscanf (fp, "%u\n", &nelems);
  assert (nelems);
  mesh.verts0.resize (nelems);

  float tmpf [3];
  for (unsigned int i = 0; i < nelems; ++i){
    status = fscanf (fp, "%f %f %f\n", &(tmpf [0]), &(tmpf [1]), &(tmpf [2]));
    mesh.verts0 [i] = vec4 (tmpf);
  }
  fclose (fp);

  mesh.verts1 = mesh.verts0;

  filename = prefix;
  filename.append ("0.trio.ele");

  fp = fopen (filename.c_str (), "r");
  nelems = 0;
  status = fscanf (fp, "%u\n", &nelems);
  assert (nelems);

  mesh.inds.resize (3*nelems);

  unsigned int tmpu [3];
  for (unsigned int i = 0; i < nelems; ++i){
    status = fscanf (fp, "%u %u %u\n", &(tmpu [0]), &(tmpu [1]), &(tmpu [2]));
    for (unsigned int j = 0; j < 3; ++j){
      mesh.inds [3*i + j] = tmpu [j];
    }
  }

  fclose (fp);
}

int main (int argc, char **argv)
{
  readCube ();

  mutexes [0] = new boost::interprocess::interprocess_semaphore (0);
  mutexes [1] = new boost::interprocess::interprocess_semaphore (1);

  windowWidth = 512;
  windowHeight = 512;

  // initialize OpenGL window
  glutInit (&argc, argv);
  glutInitWindowPosition (0, 0);
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitDisplayMode (GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
  glutCreateWindow ("Simulate OpenGL GL_Window");

  // enable GL states
  glShadeModel (GL_SMOOTH);
  glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

  glEnable (GL_DEPTH_TEST);
  glDepthFunc (GL_LESS);

  glEnable (GL_NORMALIZE);
  glEnable (GL_POLYGON_SMOOTH);

  glEnable (GL_CULL_FACE);
  glCullFace (GL_BACK);

  initGL ();
  initGPUProgs ();

  gl_context = glXGetCurrentContext ();
  gl_display = glXGetCurrentDisplay ();
  gl_drawable = glXGetCurrentDrawable ();

  // initialize the CUDA-responsible thread
  boost::thread myt (&run);
  myThread = boost::move (myt);

  // GLUT function binding
  glutDisplayFunc (display);
  glutIdleFunc (idle);
  glutReshapeFunc (resize);
  glutKeyboardFunc (keys);

  glutMainLoop ();
}
