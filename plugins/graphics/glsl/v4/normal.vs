in vec4 vertex;
in vec2 inTexCoord;

out vec2 texCoord;

void main ()
{
  gl_Position = vertex;
  texCoord = inTexCoord;
}
