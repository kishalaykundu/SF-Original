in vec4 vertex;
in vec4 restVertex;
in vec2 inTexCoord;

out vec4 restPosition;
out vec2 texCoord;

void main ()
{
  gl_Position = vertex;
  restPosition = restVertex;
  texCoord = inTexCoord;
}
