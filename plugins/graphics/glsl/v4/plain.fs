uniform mat4 projection;
uniform vec3 color;

in vec4 position;
in vec4 eye;
in vec4 normal;

uniform samplerCube envTexture;

uniform vec3 lightDir1;
uniform vec3 lightAmbient1;
uniform vec3 lightDiffuse1;
uniform float lightSpecular1;
uniform float lightExp1;

const vec2 RF = vec2 (.01, 1.);

out vec4 fragColor;

void main ()
{
  // get contribution from light
  vec3 L = normalize (-lightDir1);
  vec3 N = normalize (normal.xyz);
  vec3 V = normalize (eye.xyz*position.w - position.xyz*eye.w);
  vec3 H = normalize (L + V);

  float diff = max (0., dot (N, L));
  float spec = pow (max (0., dot (N, H)), lightExp1);

  vec3 lightSpec1 = vec3 (lightSpecular1);
  vec3 sColor3 = lightAmbient1 + color.rgb + diff*lightDiffuse1 + spec*lightSpec1;
  vec4 surfaceColor = vec4 (sColor3, 1.);

  // get contribution from environment
  vec3 R = normalize (reflect (-V, N));
  vec4 envTexColor = texture (envTexture, R);
  float scale = length (envTexColor.xyz) > 2. ? 2./ length (envTexColor.xyz) : 1.;
  vec4 envColor = vec4 (scale*envTexColor.xyz, 1.);

  // calculate the fresnel term
  float fresnel = mix (RF.x, RF.y, pow (1. - dot (N, V), 5.));

  // mix the two colors
  fragColor = mix (surfaceColor, envColor, fresnel);
}
