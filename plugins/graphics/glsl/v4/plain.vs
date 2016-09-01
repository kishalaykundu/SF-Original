uniform mat4 modelview;
uniform mat4 projection;

in vec4 vertex;

in vec2 normalTexCoord;
uniform sampler2D normalTexture;

out vec4 position;
out vec4 eye;
out vec4 normal;

void main ()
{
  position = modelview * vertex;

  vec4 nrm = normalize (texture (normalTexture, normalTexCoord));
  normal = modelview * vec4 (nrm.xyz, 0.);

  eye = inverse (projection) * vec4 (0., 0., -1., 0.);

  gl_Position = projection * position;
}
