layout (triangles) in;
layout (points, max_vertices = 3) out;

in vec2 texCoord [3];
out vec4 normal;

void main ()
{
  int i;
  vec3 edge1 = gl_in [1].gl_Position.xyz - gl_in [0].gl_Position.xyz;
  vec3 edge2 = gl_in [2].gl_Position.xyz - gl_in [0].gl_Position.xyz;

  vec3 product = cross (edge1, edge2);

  for (i = 0; i < 3; ++i){
    gl_Position = vec4 (texCoord [i].x * 2. - 1., texCoord [i].y * 2. - 1., 0., 1.);
    normal = vec4 (product, 1.);

    EmitVertex ();
  }

  EndPrimitive ();
}
