uniform sampler3D colorTexture;
uniform sampler2D texCoordTexture;

in vec4 position;
in vec3 normal;
in vec4 eye;

in vec4 trigVerts[ 3 ];
in vec2 trigTexCoordCoords[ 3 ];
in float trigSurfaceFlags[ 3 ];
in vec3 trigTexCoords[ 3 ];

in vec2 outTexCoordCoord;
in vec3 outTexCoord;
in vec2 uvCoords;

// lighting parameters -- could pass in as uniform
const vec3 ldir = vec3( .4, -.4, -1. );
const vec3 lambient = vec3( 0., .1, 0. );
const vec3 ldiffuse = vec3( 0., .2, 0. );
const vec3 lspecular = vec3( .5 );
const float lexp = 50.;

out vec4 fragColor;

void main()
{
	int i;
	int counter = 0;
	for( i = 0; i < 3; ++i ){
		if( trigSurfaceFlags[ i ] > .5 ){
			++counter;
		}
	}

	vec2 tcc;
	vec3 tc, tmp;

	/******** ALL SURFACE VERTICES ********/
	if( counter == 3 ){
		tc = texture( texCoordTexture, outTexCoordCoord ).xyz;
		vec4 texColor = texture( colorTexture, tc );

		// get contribution from light
		vec3 L = normalize( -ldir );
		vec3 N = normalize( normal*0.5 );
		vec3 V = normalize( eye.xyz*position.w - position.xyz*eye.w );
		vec3 H = normalize( L + V );
	
		float diff = max( 0., dot( N, L ));
		float spec = pow( max( 0., dot( N, H )), lexp );

		vec3 col = lambient + texColor.rgb/texColor.a + ldiffuse*diff + lspecular*spec;

		fragColor = vec4( col , 1. ); return;
	}
	/******** TWO SURFACE VERTICES ********/
	else if( counter == 2 ){
		int index1, index2, index3;

		if( trigSurfaceFlags[ 0 ] > .5 ){
			index1 = 0;
			if( trigSurfaceFlags[ 1 ] > .5 ){
				index2 = 1;
				index3 = 2;
			} else {
				index2 = 2;
				index3 = 1;
			}
		}
		else {
			index1 = 1;
			index2 = 2;
			index3 = 0;
		}

		vec3 edge = position.xyz - trigVerts[ index3 ].xyz;
		vec3 V = normalize( cross( normal, edge ));

		vec3 edge1 = trigVerts[ index3 ].xyz - trigVerts[ index1 ].xyz;
		vec3 edge2 = trigVerts[ index2 ].xyz - trigVerts[ index1 ].xyz;
		float ratio = dot( V, edge1 )/ dot( V, edge2 );
	
		tcc = mix( trigTexCoordCoords[ index1 ], trigTexCoordCoords[ index2 ], ratio );
		vec3 pointTexCoord = texture( texCoordTexture, tcc ).xyz;

		vec4 point = mix( trigVerts[ index1 ], trigVerts[ index2 ], ratio );
		vec3 edge3 = point.xyz - position.xyz;

		ratio = length( edge )/( length( edge ) + length( edge3 ));

		tc = mix( trigTexCoords[ index3 ], pointTexCoord, ratio );
	}
	/******** ONE SURFACE VERTEX ********/
	else if( counter == 1 ){

		if( trigSurfaceFlags[ 0 ] > 0.5 ){
			vec4 texColr = texture( texCoordTexture, trigTexCoordCoords[ 0 ]);
			tc = texColr.xyz*uvCoords.x + 
				trigTexCoords[ 1 ]*uvCoords.y + trigTexCoords[ 2 ]*( 1. - uvCoords.x - uvCoords.y );
		}
		else if( trigSurfaceFlags[ 1 ] > 0.5 ){
			vec4 texColr = texture( texCoordTexture, trigTexCoordCoords[ 1 ]);
			tc = trigTexCoords[ 0 ]*uvCoords.x + texColr.xyz*uvCoords.y +
				trigTexCoords[ 2 ]*( 1. - uvCoords.x - uvCoords.y );
		}
		else {
			vec4 texColr = texture( texCoordTexture, trigTexCoordCoords[ 2 ]);
			tc = trigTexCoords[ 0 ]*uvCoords.x + trigTexCoords[ 1 ]*uvCoords.y +
				texColr.xyz*( 1. - uvCoords.x - uvCoords.y );
		}
	}
	/******** ALL INTERNAL VERTICES ********/
	else {
		tc = outTexCoord;
	}

	tc = outTexCoord;
	vec4 col = texture( colorTexture, tc );
	fragColor = vec4( col.rgb, 1. );
}
