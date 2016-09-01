layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

uniform mat4 projection;
uniform mat4 modelview;

in vec2 texCoordCoord[ 3 ];
in float surfaceFlag[ 3 ];
in vec3 texCoord[ 3 ];
in vec4 ineye[ 3 ];

out vec4 position;
out vec3 normal;
out vec4 eye;
out vec4 trigVerts[ 3 ];
out vec2 trigTexCoordCoords[ 3 ];
out float trigSurfaceFlags[ 3 ];
out vec3 trigTexCoords[ 3 ];
out vec2 outTexCoordCoord;
out vec3 outTexCoord;
out vec2 uvCoords;

void main( )
{
	int i;
	for( i = 0; i < 3; ++i ){
		
		trigVerts[ i ] = gl_in[ i ].gl_Position;
		trigTexCoordCoords[ i ] = texCoordCoord[ i ];
		trigSurfaceFlags[ i ] = surfaceFlag[ i ];
		trigTexCoords[ i ] = texCoord[ i ];
	}

	vec3 edge1 = gl_in[ 1 ].gl_Position.xyz - gl_in[ 0 ].gl_Position.xyz;
	vec3 edge2 = gl_in[ 2 ].gl_Position.xyz - gl_in[ 0 ].gl_Position.xyz;
	vec3 nrm = normalize( cross( edge1, edge2 ));
	vec4 tnormal = modelview * vec4( nrm, 0. );
	normal = tnormal.xyz;

	for( i = 0; i < 3; ++i ){
		
		if( i == 0 ){
			uvCoords = vec2( 1., 0. );
		} else if( i == 1 ){
			uvCoords = vec2( 0., 1. );
		} else {
			uvCoords = vec2( 0., 0. );
		}

		position = gl_in[ i ].gl_Position;
		outTexCoord = texCoord[ i ];
		outTexCoordCoord = texCoordCoord[ i ];
		eye = ineye[ i ];
		gl_Position = projection * gl_in[ i ].gl_Position;
		EmitVertex( );
	}
	
	EndPrimitive( );
}
