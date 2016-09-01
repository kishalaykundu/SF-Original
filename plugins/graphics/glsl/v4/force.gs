layout (lines) in;
layout (points, max_vertices = 2) out;

in vec4 restPosition [2];
in vec2 texCoord [2];

out vec4 force;

void main ()
{
  int i;
  vec3 currentLength = gl_in [1].gl_Position.xyz - gl_in [0].gl_Position.xyz;
  vec3 restLength = restPosition [1].xyz - restPosition [0].xyz;

  vec3 difference = restLength - currentLength;

  for (i = 0; i < 2; ++i){
    gl_Position = vec4 (texCoord [i].x * 2. - 1., texCoord [i].y * 2. - 1., 0., 1.);
    force = vec4 (difference, 1.);

    EmitVertex ();
  }

  EndPrimitive ();
}
