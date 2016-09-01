uniform mat4 modelview;
uniform mat4 projection;

in vec4 vertex;
in vec2 inTexCoordCoord;
in float inSurfaceFlag;
in vec3 inTexCoord;

out vec2 texCoordCoord;
out float surfaceFlag;
out vec3 texCoord;
out vec4 ineye;

void main( )
{
  texCoordCoord = inTexCoordCoord;
  surfaceFlag = inSurfaceFlag;
  texCoord = inTexCoord;
	ineye = inverse( projection ) * vec4( 0., 0., -1., 0. );
	gl_Position = modelview * vertex;
}
