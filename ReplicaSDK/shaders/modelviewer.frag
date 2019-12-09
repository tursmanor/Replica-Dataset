#version 130
#expect SHOW_COLOR
#expect SHOW_NORMAL
#expect SHOW_TEXTURE
#expect SHOW_MATCAP
#expect SHOW_UV
#expect SHOW_DEPTH
#expect SHOW_BINARY

#if SHOW_COLOR
    varying vec4 vColor;
#elif SHOW_NORMAL
    varying vec3 vNormal;
#elif SHOW_TEXTURE
    varying vec2 vUV;
    uniform sampler2D texture_0;
#elif SHOW_MATCAP
    varying vec3 vNormalCam;
    uniform sampler2D matcap;
#elif SHOW_UV
    varying vec2 vUV;
#elif SHOW_DEPTH
    varying float depth;
#else
    varying vec3 vP;
#endif

void main() {
#if SHOW_COLOR
    gl_FragColor = vColor;
#elif SHOW_NORMAL
    gl_FragColor = vec4((vNormal + vec3(1.0,1.0,1.0)) / 2.0, 1.0);
#elif SHOW_TEXTURE
    gl_FragColor = texture2D(texture_0, vUV);
#elif SHOW_MATCAP
    vec2 uv = 0.5 * vNormalCam.xy + vec2(0.5, 0.5);
    gl_FragColor = texture2D(matcap, uv);
#elif SHOW_UV
    gl_FragColor = vec4(vUV,1.0-vUV.x,1.0);
#elif SHOW_DEPTH
    //float mDepth = depth * 1.0f/16.0f;
    float mDepth = depth;
    gl_FragColor = vec4(mDepth, mDepth, mDepth, 1.0f);
#elif SHOW_BINARY
    gl_FragColor = vec4(1.0f,1.0f,1.0f,1.0f); // Make the mesh/foreground white
#else
    gl_FragColor = vec4(vP / 100.0,1.0);
#endif
}