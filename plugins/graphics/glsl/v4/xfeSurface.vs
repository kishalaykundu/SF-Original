uniform mat4 modelview;
uniform mat4 projection;

in vec4 vertex;
in vec2 inColorTexCoord;

in vec2 normalTexCoord;
uniform sampler2D normalTexture;

out vec4 position;
out vec4 normal;
out vec4 eye;
out vec2 colorTexCoord;

void main ()
{
  position = modelview * vertex;

  vec4 tmpNormal = normalize (texture (normalTexture, normalTexCoord));
  normal = modelview * vec4(tmpNormal.xyz, 0.);

  eye = inverse (projection) * vec4 (0., 0., -1., 0.);

  colorTexCoord = inColorTexCoord;

  gl_Position = projection * position;
}
