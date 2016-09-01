in vec2 vertex;
in vec4 texCoord;

out vec4 color;

void main ()
{
  color = texCoord;
  gl_Position = vec4 (vertex * 2. - 1., 0., 1.);
}
